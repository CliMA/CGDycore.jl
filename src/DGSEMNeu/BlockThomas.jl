using StaticArrays
using LinearAlgebra

"""
    block_thomas_optimized!(X, A, B, C, D, C_prime, D_prime)

Löst das block-tridiagonale System in-place ohne Speicherallokationen.
Nutzt StaticArrays für maximale CPU-Registernutzung und Loop-Unrolling.
"""
function block_thomas_optimized!(X, A, B, C, D, C_prime, D_prime)
    n = length(B)
    
    # --- 1. Vorwärtselimination ---
    
    # Erster Block (i = 1)
    # B[1] \ C[1] löst das System exakt via Cramerscher Regel/unrolled LU auf Registern
    C_prime[1] = B[1] \ C[1]
    D_prime[1] = B[1] \ D[1]
    
    # Folgende Blöcke (i = 2 bis n)
    @inbounds for i in 2:n
        # Keine Allokation: Operationen auf SMatrix/SVector erzeugen keinen Heap-Speicher
        Gamma = B[i] - A[i] * C_prime[i-1]
        
        if i < n
            C_prime[i] = Gamma \ C[i]
        end
        
        D_prime[i] = Gamma \ (D[i] - A[i] * D_prime[i-1])
    end

    # --- 2. Rückwärtssubstitution ---
    
    # Letzter Block (i = n)
    X[n] = D_prime[n]
    
    # Von n-1 rückwärts bis 1
    @inbounds for i in (n-1):-1:1
        X[i] = D_prime[i] - C_prime[i] * X[i+1]
    end
    
    return nothing
end

# ==============================================================================
# BENCHMARK & PERFORMANCE-TEST
# ==============================================================================
using BenchmarkTools

n = 1000 # Anzahl der 3x3 Blöcke

# Definition der statischen Typen für 3x3 Matrizen und 3x1 Vektoren
const BlockMat = SMatrix{3, 3, Float64, 9}
const BlockVec = SVector{3, Float64}

# Eingabedaten initialisieren (Daten liegen kontinuierlich im Speicher)
A = [rand(BlockMat) for _ in 1:n]
B = [rand(BlockMat) + 4I for _ in 1:n] 
C = [rand(BlockMat) for _ in 1:n]
D = [rand(BlockVec) for _ in 1:n]

# Vorallokation der Arbeits- und Ergebnis-Container (nur einmalig vor der Schleife)
X_res   = Vector{BlockVec}(undef, n)
C_tmp   = Vector{BlockMat}(undef, n)
D_tmp   = Vector{BlockVec}(undef, n)

# Benchmark ausführen
println("--- Performance-Ergebnis ---")
@btime block_thomas_optimized!($X_res, $A, $B, $C, $D, $C_tmp, $D_tmp)
