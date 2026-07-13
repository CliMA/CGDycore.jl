using KernelAbstractions
using StaticArrays
using LinearAlgebra

# Typendefinitionen für die 3x3 Blöcke im Thread-Register
const BlockMat = SMatrix{3, 3, Float64, 9}
const BlockVec = SVector{3, Float64}

"""
    mul_A_structure(A, C)
Spezifischer Multiplikationskern: Nutzt aus, dass die erste Spalte von A null ist.
Reduziert FLOPs direkt im Register des GPU/CPU-Threads.
"""
@inline function mul_A_structure(A::BlockMat, C::BlockMat)
    return BlockMat(
        A[1,2]*C[2,1] + A[1,3]*C[3,1], A[1,2]*C[2,2] + A[1,3]*C[3,2], A[1,2]*C[2,3] + A[1,3]*C[3,3],
        A[2,2]*C[2,1] + A[2,3]*C[3,1], A[2,2]*C[2,2] + A[2,3]*C[3,2], A[2,2]*C[2,3] + A[2,3]*C[3,3],
        A[3,2]*C[2,1] + A[3,3]*C[3,1], A[3,2]*C[2,2] + A[3,3]*C[3,2], A[3,2]*C[2,3] + A[3,3]*C[3,3]
    )
end

@inline function mul_A_structure(A::BlockMat, v::BlockVec)
    return BlockVec(
        A[1,2]*v[2] + A[1,3]*v[3],
        A[2,2]*v[2] + A[2,3]*v[3],
        A[3,2]*v[2] + A[3,3]*v[3]
    )
end

# Definition des plattformunabhängigen Rechenkerns
@kernel function block_thomas_3d_kernel!(X, @Const(A), @Const(B), @Const(C), @Const(D), n)
    # JEDER Thread kümmert sich um genau eine (x, y) Spalte des 3D-Gitters
    x, y = @index(Global, NTuple)
    
    # Lokaler Thread-Speicher (wird vom Compiler in Registern gehalten, da Größe n zur Kompilierzeit feststeht)
    # Da n < 25, passen diese Arrays vollständig in den schnellsten Speicher des Streaming-Multiprozessors (SM)
    C_prime = MVector{24, BlockMat}(undef)
    D_prime = MVector{24, BlockVec}(undef)

    # --- 1. Vorwärtselimination ---
    C_prime[1] = B[1, x, y] \ C[1, x, y]
    D_prime[1] = B[1, x, y] \ D[1, x, y]

    for z in 2:n
        Gamma = B[z, x, y] - mul_A_structure(A[z, x, y], C_prime[z-1])
        
        if z < n
            C_prime[z] = Gamma \ C[z, x, y]
        end
        D_prime[z] = Gamma \ (D[z, x, y] - mul_A_structure(A[z, x, y], D_prime[z-1]))
    end

    # --- 2. Rückwärtssubstitution ---
    X_local = D_prime[n]
    X[n, x, y] = X_local

    for z in (n-1):-1:1
        X_local = D_prime[z] - C_prime[z] * X_local
        X[z, x, y] = X_local
    end
end

# ==============================================================================
# TREIBER-FUNKTION (Wechselt transparent zwischen CPU und GPU)
# ==============================================================================
function solve_3d_system(backend, NX, NY, NZ)
    n = NZ # Anzahl der Blöcke (z.B. 20)
    
    # 1. Speicher auf dem gewählten Backend allokieren (CPU, CUDA, AMD, etc.)
    # In KA werden Arrays passend zum Backend erzeugt
    A_device = KernelAbstractions.allocate(backend, BlockMat, n, NX, NY)
    B_device = KernelAbstractions.allocate(backend, BlockMat, n, NX, NY)
    C_device = KernelAbstractions.allocate(backend, BlockMat, n, NX, NY)
    D_device = KernelAbstractions.allocate(backend, BlockVec, n, NX, NY)
    X_device = KernelAbstractions.allocate(backend, BlockVec, n, NX, NY)

    # (Hier würden Sie Ihre echten Simulationsdaten auf das Device kopieren oder dort initialisieren)
    # Aus Veranschaulichungsgründen füllen wir die Daten hier fiktiv:
    KernelAbstractions.fill!(A_device, rand(BlockMat))
    KernelAbstractions.fill!(B_device, rand(BlockMat) + 4I)
    KernelAbstractions.fill!(C_device, rand(BlockMat))
    KernelAbstractions.fill!(D_device, rand(BlockVec))

    # 2. Kernel für das spezifische Backend kompilieren
    kernel! = block_thomas_3d_kernel!(backend)

    # 3. Kernel asynchron aufrufen
    # Die ND-Range entspricht der (X, Y) Ausdehnung unseres Gitters
    # Work-Group-Size (z. B. 16x16 Threads pro Block auf der GPU)
    work_groups = (16, 16)
    nd_range    = (NX, NY)
    
    # Kernel ausführen
    kernel!(X_device, A_device, B_device, C_device, D_device, n, ndrange=nd_range, workgroupsize=work_groups)
    
    # 4. Auf Fertigstellung der GPU/CPU warten
    KernelAbstractions.synchronize(backend)
    
    return X_device
end

# ==============================================================================
# AUFRUF-BEISPIEL
# ==============================================================================

# Für CPU-Ausführung (Multi-Threading):
using CPUBackend
backend_cpu = CPU()
X_cpu = solve_3d_system(backend_cpu, 128, 128, 20)
println("CPU Kernel erfolgreich berechnet.")

# Für NVIDIA-GPU-Ausführung (einfach einkommentieren, falls CUDA verfügbar ist):
# using CUDA
# backend_gpu = CUDA.CUDABackend()
# X_gpu = solve_3d_system(backend_gpu, 128, 128, 20)

