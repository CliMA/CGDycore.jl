using AMDGPU
using KernelAbstractions
@kernel inbounds = true function MomentumKernelAtomic!(F, U, D, dXdxI, MRho, M, Glob)
    I, J, iz = @index(Local, NTuple)
    _, _, Iz, IF = @index(Global, NTuple)
    ColumnTilesDim = @uniform @groupsize()[3]
    N = @uniform @groupsize()[1]
    Nz = @uniform @ndrange()[3]
    ID = I + (J - 1) * N
    ind = Glob[ID, IF]
    RhoCol = @localmem eltype(F) (N, N, ColumnTilesDim)
    uCol = @localmem eltype(F) (N, N, ColumnTilesDim)
    vCol = @localmem eltype(F) (N, N, ColumnTilesDim)
    wCol = @localmem eltype(F) (N, N, ColumnTilesDim + 1)
    if Iz <= Nz
        RhoCol[I, J, iz] = U[Iz, ind, 1]
        uCol[I, J, iz] = U[Iz, ind, 2]
        vCol[I, J, iz] = U[Iz, ind, 3]
        wCol[I, J, iz+1] = U[Iz, ind, 4]
        if Iz == 1
            wCol[I, J, 1] = -(dXdxI[3, 1, 1, ID, 1, IF] * U[Iz, ind, 2] +
                              dXdxI[3, 2, 1, ID, 1, IF] * U[Iz, ind, 3]) / dXdxI[3, 3, 1, ID, 1, IF]
        elseif iz == 1
            wCol[I, J, 1] = U[Iz-1, ind, 4]
        end
    end
    @synchronize
    ID = I + (J - 1) * N
    ind = Glob[ID, IF]
    if Iz <= Nz
        ind = Glob[ID, IF]
        uCon1 = -RhoCol[I, J, iz] * (dXdxI[1, 1, 1, ID, Iz, IF] * uCol[I, J, iz] +
                                     dXdxI[1, 2, 1, ID, Iz, IF] * vCol[I, J, iz] + dXdxI[1, 3, 1, ID, Iz, IF] * wCol[I, J, iz])
        uCon2 = -RhoCol[I, J, iz] * (dXdxI[1, 1, 2, ID, Iz, IF] * uCol[I, J, iz] +
                                     dXdxI[1, 2, 2, ID, Iz, IF] * vCol[I, J, iz] + dXdxI[1, 3, 2, ID, Iz, IF] * wCol[I, J, iz+1])
        vCon1 = -RhoCol[I, J, iz] * (dXdxI[2, 1, 1, ID, Iz, IF] * uCol[I, J, iz] +
                                     dXdxI[2, 2, 1, ID, Iz, IF] * vCol[I, J, iz] + dXdxI[2, 3, 1, ID, Iz, IF] * wCol[I, J, iz])
        vCon2 = -RhoCol[I, J, iz] * (dXdxI[2, 1, 2, ID, Iz, IF] * uCol[I, J, iz] +
                                     dXdxI[2, 2, 2, ID, Iz, IF] * vCol[I, J, iz] + dXdxI[2, 3, 2, ID, Iz, IF] * wCol[I, J, iz+1])
        wCon1 = -RhoCol[I, J, iz] * (dXdxI[3, 1, 1, ID, Iz, IF] * uCol[I, J, iz] +
                                     dXdxI[3, 2, 1, ID, Iz, IF] * vCol[I, J, iz] + dXdxI[3, 3, 1, ID, Iz, IF] * wCol[I, J, iz])
        wCon2 = -RhoCol[I, J, iz] * (dXdxI[3, 1, 2, ID, Iz, IF] * uCol[I, J, iz] +
                                     dXdxI[3, 2, 2, ID, Iz, IF] * vCol[I, J, iz] + dXdxI[3, 3, 2, ID, Iz, IF] * wCol[I, J, iz+1])
        Dxu = D[I, 1] * uCol[1, J, iz]
        Dyu = D[J, 1] * uCol[I, 1, iz]
        Dxv = D[I, 1] * vCol[1, J, iz]
        Dyv = D[J, 1] * vCol[I, 1, iz]
        Dxw1 = D[I, 1] * wCol[1, J, iz]
        Dyw1 = D[J, 1] * wCol[I, 1, iz]
        Dxw2 = D[I, 1] * wCol[1, J, iz+1]
        Dyw2 = D[J, 1] * wCol[I, 1, iz+1]
        Izp = min(Iz + 1, Nz)
        Izm = max(Iz - 1, 1)
        ind = Glob[ID, IF]
        Dzu2 = eltype(F)(0.5) * (U[Izp, ind, 2] - uCol[I, J, iz])
        Dzv2 = eltype(F)(0.5) * (U[Izp, ind, 3] - vCol[I, J, iz])
        Dzu1 = eltype(F)(0.5) * (uCol[I, J, iz] - U[Izm, ind, 2])
        Dzv1 = eltype(F)(0.5) * (vCol[I, J, iz] - U[Izm, ind, 3])
        Dzw = eltype(F)(0.5) * (wCol[I, J, iz+1] - wCol[I, J, iz])
        for k = 2:N
            Dxu += D[I, k] * uCol[k, J, iz]
            Dyu += D[J, k] * uCol[I, k, iz]
            Dxv += D[I, k] * vCol[k, J, iz]
            Dyv += D[J, k] * vCol[I, k, iz]
            Dxw1 += D[I, k] * wCol[k, J, iz]
            Dyw1 += D[J, k] * wCol[I, k, iz]
            Dxw2 += D[I, k] * wCol[k, J, iz+1]
            Dyw2 += D[J, k] * wCol[I, k, iz+1]
        end
        KernelAbstractions.@atomic F[Iz, ind, 2] += ((uCon1 + uCon2) * Dxu + (vCon1 + vCon2) * Dyu +
                                  wCon1 * Dzu1 + wCon2 * Dzu2) / M[Iz, ind] / RhoCol[I, J, iz]
        KernelAbstractions.@atomic F[Iz, ind, 3] += ((uCon1 + uCon2) * Dxv + (vCon1 + vCon2) * Dyv +
                                  wCon1 * Dzv1 + wCon2 * Dzv2) / M[Iz, ind] / RhoCol[I, J, iz]
    end
    if Iz > 1
        KernelAbstractions.@atomic F[Iz-1, ind, 4] += (uCon1 * Dxw1 + vCon1 * Dyw1 + wCon1 * Dzw) / MRho[Iz-1, ind]
    end
    if Iz < Nz
        KernelAbstractions.@atomic F[Iz, ind, 4] += (uCon2 * Dxw2 + vCon2 * Dyw2 + wCon2 * Dzw) / MRho[Iz, ind]
    end
end

@kernel inbounds = true function MomentumKernelAtomicMonotone!(F, U, D, dXdxI, MRho, M, Glob)
    I, J, iz = @index(Local, NTuple)
    _, _, Iz, IF = @index(Global, NTuple)
    ColumnTilesDim = @uniform @groupsize()[3]
    N = @uniform @groupsize()[1]
    Nz = @uniform @ndrange()[3]
    ID = I + (J - 1) * N
    ind = Glob[ID, IF]
    RhoCol = @localmem eltype(F) (N, N, ColumnTilesDim)
    uCol = @localmem eltype(F) (N, N, ColumnTilesDim)
    vCol = @localmem eltype(F) (N, N, ColumnTilesDim)
    wCol = @localmem eltype(F) (N, N, ColumnTilesDim + 1)
    if Iz <= Nz
        RhoCol[I, J, iz] = U[Iz, ind, 1]
        uCol[I, J, iz] = U[Iz, ind, 2]
        vCol[I, J, iz] = U[Iz, ind, 3]
        wCol[I, J, iz+1] = U[Iz, ind, 4]
        if Iz == 1
            wCol[I, J, 1] = -(dXdxI[3, 1, 1, ID, 1, IF] * U[Iz, ind, 2] +
                              dXdxI[3, 2, 1, ID, 1, IF] * U[Iz, ind, 3]) / dXdxI[3, 3, 1, ID, 1, IF]
        elseif iz == 1
            wCol[I, J, 1] = U[Iz-1, ind, 4]
        end
    end
    @synchronize
    ID = I + (J - 1) * N
    ind = Glob[ID, IF]
    if Iz <= Nz
        ind = Glob[ID, IF]
        uCon1 = -RhoCol[I, J, iz] * (dXdxI[1, 1, 1, ID, Iz, IF] * uCol[I, J, iz] +
                                     dXdxI[1, 2, 1, ID, Iz, IF] * vCol[I, J, iz] + dXdxI[1, 3, 1, ID, Iz, IF] * wCol[I, J, iz])
        uCon2 = -RhoCol[I, J, iz] * (dXdxI[1, 1, 2, ID, Iz, IF] * uCol[I, J, iz] +
                                     dXdxI[1, 2, 2, ID, Iz, IF] * vCol[I, J, iz] + dXdxI[1, 3, 2, ID, Iz, IF] * wCol[I, J, iz+1])
        vCon1 = -RhoCol[I, J, iz] * (dXdxI[2, 1, 1, ID, Iz, IF] * uCol[I, J, iz] +
                                     dXdxI[2, 2, 1, ID, Iz, IF] * vCol[I, J, iz] + dXdxI[2, 3, 1, ID, Iz, IF] * wCol[I, J, iz])
        vCon2 = -RhoCol[I, J, iz] * (dXdxI[2, 1, 2, ID, Iz, IF] * uCol[I, J, iz] +
                                     dXdxI[2, 2, 2, ID, Iz, IF] * vCol[I, J, iz] + dXdxI[2, 3, 2, ID, Iz, IF] * wCol[I, J, iz+1])
        wCon1 = -RhoCol[I, J, iz] * (dXdxI[3, 1, 1, ID, Iz, IF] * uCol[I, J, iz] +
                                     dXdxI[3, 2, 1, ID, Iz, IF] * vCol[I, J, iz] + dXdxI[3, 3, 1, ID, Iz, IF] * wCol[I, J, iz])
        wCon2 = -RhoCol[I, J, iz] * (dXdxI[3, 1, 2, ID, Iz, IF] * uCol[I, J, iz] +
                                     dXdxI[3, 2, 2, ID, Iz, IF] * vCol[I, J, iz] + dXdxI[3, 3, 2, ID, Iz, IF] * wCol[I, J, iz+1])
        Dxu = D[I, 1] * uCol[1, J, iz]
        Dyu = D[J, 1] * uCol[I, 1, iz]
        Dxv = D[I, 1] * vCol[1, J, iz]
        Dyv = D[J, 1] * vCol[I, 1, iz]
        Dxw1 = D[I, 1] * wCol[1, J, iz]
        Dyw1 = D[J, 1] * wCol[I, 1, iz]
        Dxw2 = D[I, 1] * wCol[1, J, iz+1]
        Dyw2 = D[J, 1] * wCol[I, 1, iz+1]
        Izp = min(Iz + 1, Nz)
        Izm = max(Iz - 1, 1)
        ind = Glob[ID, IF]
        Dzu2 = eltype(F)(0.5) * (U[Izp, ind, 2] - uCol[I, J, iz])
        Dzv2 = eltype(F)(0.5) * (U[Izp, ind, 3] - vCol[I, J, iz])
        Dzu1 = eltype(F)(0.5) * (uCol[I, J, iz] - U[Izm, ind, 2])
        Dzv1 = eltype(F)(0.5) * (vCol[I, J, iz] - U[Izm, ind, 3])
        Dzw = eltype(F)(0.5) * (wCol[I, J, iz+1] - wCol[I, J, iz])
        for k = 2:N
            Dxu += D[I, k] * uCol[k, J, iz]
            Dyu += D[J, k] * uCol[I, k, iz]
            Dxv += D[I, k] * vCol[k, J, iz]
            Dyv += D[J, k] * vCol[I, k, iz]
            Dxw1 += D[I, k] * wCol[k, J, iz]
            Dyw1 += D[J, k] * wCol[I, k, iz]
            Dxw2 += D[I, k] * wCol[k, J, iz+1]
            Dyw2 += D[J, k] * wCol[I, k, iz+1]
        end
        KernelAbstractions.@atomic :monotonic F[Iz, ind, 2] += ((uCon1 + uCon2) * Dxu + (vCon1 + vCon2) * Dyu +
                                  wCon1 * Dzu1 + wCon2 * Dzu2) / M[Iz, ind] / RhoCol[I, J, iz]
        KernelAbstractions.@atomic :monotonic F[Iz, ind, 3] += ((uCon1 + uCon2) * Dxv + (vCon1 + vCon2) * Dyv +
                                  wCon1 * Dzv1 + wCon2 * Dzv2) / M[Iz, ind] / RhoCol[I, J, iz]
    end
    if Iz > 1
        KernelAbstractions.@atomic :monotonic F[Iz-1, ind, 4] += (uCon1 * Dxw1 + vCon1 * Dyw1 + wCon1 * Dzw) / MRho[Iz-1, ind]
    end
    if Iz < Nz
        KernelAbstractions.@atomic :monotonic F[Iz, ind, 4] += (uCon2 * Dxw2 + vCon2 * Dyw2 + wCon2 * Dzw) / MRho[Iz, ind]
    end
end

@kernel inbounds = true function MomentumKernel!(F, U, D, dXdxI, MRho, M, Glob)
    I, J, iz = @index(Local, NTuple)
    _, _, Iz, IF = @index(Global, NTuple)
    ColumnTilesDim = @uniform @groupsize()[3]
    N = @uniform @groupsize()[1]
    Nz = @uniform @ndrange()[3]
    NF = @uniform @ndrange()[4]
    ID = I + (J - 1) * N
    ind = Glob[ID, IF]
    RhoCol = @localmem eltype(F) (N, N, ColumnTilesDim)
    uCol = @localmem eltype(F) (N, N, ColumnTilesDim)
    vCol = @localmem eltype(F) (N, N, ColumnTilesDim)
    wCol = @localmem eltype(F) (N, N, ColumnTilesDim + 1)
    if Iz <= Nz
        RhoCol[I, J, iz] = U[Iz, ind, 1]
        uCol[I, J, iz] = U[Iz, ind, 2]
        vCol[I, J, iz] = U[Iz, ind, 3]
        wCol[I, J, iz+1] = U[Iz, ind, 4]
        if Iz == 1
            wCol[I, J, 1] = -(dXdxI[3, 1, 1, ID, 1, IF] * U[Iz, ind, 2] +
                              dXdxI[3, 2, 1, ID, 1, IF] * U[Iz, ind, 3]) / dXdxI[3, 3, 1, ID, 1, IF]
        elseif iz == 1
            wCol[I, J, 1] = U[Iz-1, ind, 4]
        end
    end
    @synchronize
    ID = I + (J - 1) * N
    ind = Glob[ID, IF]
    if Iz <= Nz
        ind = Glob[ID, IF]
        uCon1 = -RhoCol[I, J, iz] * (dXdxI[1, 1, 1, ID, Iz, IF] * uCol[I, J, iz] +
                                     dXdxI[1, 2, 1, ID, Iz, IF] * vCol[I, J, iz] + dXdxI[1, 3, 1, ID, Iz, IF] * wCol[I, J, iz])
        uCon2 = -RhoCol[I, J, iz] * (dXdxI[1, 1, 2, ID, Iz, IF] * uCol[I, J, iz] +
                                     dXdxI[1, 2, 2, ID, Iz, IF] * vCol[I, J, iz] + dXdxI[1, 3, 2, ID, Iz, IF] * wCol[I, J, iz+1])
        vCon1 = -RhoCol[I, J, iz] * (dXdxI[2, 1, 1, ID, Iz, IF] * uCol[I, J, iz] +
                                     dXdxI[2, 2, 1, ID, Iz, IF] * vCol[I, J, iz] + dXdxI[2, 3, 1, ID, Iz, IF] * wCol[I, J, iz])
        vCon2 = -RhoCol[I, J, iz] * (dXdxI[2, 1, 2, ID, Iz, IF] * uCol[I, J, iz] +
                                     dXdxI[2, 2, 2, ID, Iz, IF] * vCol[I, J, iz] + dXdxI[2, 3, 2, ID, Iz, IF] * wCol[I, J, iz+1])
        wCon1 = -RhoCol[I, J, iz] * (dXdxI[3, 1, 1, ID, Iz, IF] * uCol[I, J, iz] +
                                     dXdxI[3, 2, 1, ID, Iz, IF] * vCol[I, J, iz] + dXdxI[3, 3, 1, ID, Iz, IF] * wCol[I, J, iz])
        wCon2 = -RhoCol[I, J, iz] * (dXdxI[3, 1, 2, ID, Iz, IF] * uCol[I, J, iz] +
                                     dXdxI[3, 2, 2, ID, Iz, IF] * vCol[I, J, iz] + dXdxI[3, 3, 2, ID, Iz, IF] * wCol[I, J, iz+1])
        Dxu = D[I, 1] * uCol[1, J, iz]
        Dyu = D[J, 1] * uCol[I, 1, iz]
        Dxv = D[I, 1] * vCol[1, J, iz]
        Dyv = D[J, 1] * vCol[I, 1, iz]
        Dxw1 = D[I, 1] * wCol[1, J, iz]
        Dyw1 = D[J, 1] * wCol[I, 1, iz]
        Dxw2 = D[I, 1] * wCol[1, J, iz+1]
        Dyw2 = D[J, 1] * wCol[I, 1, iz+1]
        Izp = min(Iz + 1, Nz)
        Izm = max(Iz - 1, 1)
        ind = Glob[ID, IF]
        Dzu2 = eltype(F)(0.5) * (U[Izp, ind, 2] - uCol[I, J, iz])
        Dzv2 = eltype(F)(0.5) * (U[Izp, ind, 3] - vCol[I, J, iz])
        Dzu1 = eltype(F)(0.5) * (uCol[I, J, iz] - U[Izm, ind, 2])
        Dzv1 = eltype(F)(0.5) * (vCol[I, J, iz] - U[Izm, ind, 3])
        Dzw = eltype(F)(0.5) * (wCol[I, J, iz+1] - wCol[I, J, iz])
        for k = 2:N
            Dxu += D[I, k] * uCol[k, J, iz]
            Dyu += D[J, k] * uCol[I, k, iz]
            Dxv += D[I, k] * vCol[k, J, iz]
            Dyv += D[J, k] * vCol[I, k, iz]
            Dxw1 += D[I, k] * wCol[k, J, iz]
            Dyw1 += D[J, k] * wCol[I, k, iz]
            Dxw2 += D[I, k] * wCol[k, J, iz+1]
            Dyw2 += D[J, k] * wCol[I, k, iz+1]
        end
        F[Iz, ind, 2] += ((uCon1 + uCon2) * Dxu + (vCon1 + vCon2) * Dyu +
                          wCon1 * Dzu1 + wCon2 * Dzu2) / M[Iz, ind] / RhoCol[I, J, iz]
        F[Iz, ind, 3] += ((uCon1 + uCon2) * Dxv + (vCon1 + vCon2) * Dyv +
                          wCon1 * Dzv1 + wCon2 * Dzv2) / M[Iz, ind] / RhoCol[I, J, iz]
    end
    if Iz > 1
        F[Iz-1, ind, 4] += (uCon1 * Dxw1 + vCon1 * Dyw1 + wCon1 * Dzw) / MRho[Iz-1, ind]
    end
    if Iz < Nz
        F[Iz, ind, 4] += (uCon2 * Dxw2 + vCon2 * Dyw2 + wCon2 * Dzw) / MRho[Iz, ind]
    end
end
NF = 5400 # 30*30*6
NumG = 48602
Ord = 4
DoF = Ord * Ord
NumV = 5
NumTr = 2
Nz = 64
NumberThreadGPU = 1024
GlobCPU = zeros(Int, DoF, NF)
read!("./TestAMD/GlobInd", GlobCPU)
#backend = ROCBackend()
backend = CPU()
FT = Float32
NzG = min(div(NumberThreadGPU, DoF), Nz)
group = (Ord, Ord, NzG, 1)
ndrange = (Ord, Ord, Nz, NF)
F = KernelAbstractions.ones(backend, FT, Nz, NumG, NumV + NumTr)
U = KernelAbstractions.ones(backend, FT, Nz, NumG, NumV + NumTr)
D = KernelAbstractions.ones(backend, FT, 4, 4)
dXdxI = KernelAbstractions.ones(backend, FT, 3, 3, 2, DoF, Nz, NF)
MRho = KernelAbstractions.ones(backend, FT, Nz - 1, NumG)
M = KernelAbstractions.ones(backend, FT, Nz, NumG)
Glob = KernelAbstractions.zeros(backend, Int, DoF, NF)
copyto!(Glob, GlobCPU)
KMomentumKernelAtomic! = MomentumKernelAtomic!(backend, group)
KMomentumKernelAtomic!(F, U, D, dXdxI, MRho, M, Glob, ndrange=ndrange)
KernelAbstractions.synchronize(backend)
@show "Atomic"
@time for iter = 1:100
    KMomentumKernelAtomic!(F, U, D, dXdxI, MRho, M, Glob, ndrange=ndrange)
    KernelAbstractions.synchronize(backend)
end
KMomentumKernelAtomicMonotone! = MomentumKernelAtomicMonotone!(backend, group)
KMomentumKernelAtomicMonotone!(F, U, D, dXdxI, MRho, M, Glob, ndrange=ndrange)
KernelAbstractions.synchronize(backend)
@show "AtomicMonotone"
@time for iter = 1:100
    KMomentumKernelAtomicMonotone!(F, U, D, dXdxI, MRho, M, Glob, ndrange=ndrange)
    KernelAbstractions.synchronize(backend)
end
KMomentumKernel! = MomentumKernel!(backend, group)
KMomentumKernel!(F, U, D, dXdxI, MRho, M, Glob, ndrange=ndrange)
KernelAbstractions.synchronize(backend)
@show "Nothing"
@time for iter = 1:100
    KMomentumKernel!(F, U, D, dXdxI, MRho, M, Glob, ndrange=ndrange)
    KernelAbstractions.synchronize(backend)
end





