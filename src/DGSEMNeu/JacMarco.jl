# =====================================================================
# JacobianCacheMarcoSplitNS  —  nested-Schur HEVI Jacobian
#
# Faithful copy of the nested-Schur JacobianCacheMarco, with:
#   * schur_m promoted to an ARRAY Schurm[iz,ID]  (= D21/fac per element),
#     because the split gravity makes D21 = G[M,M] vary per element.
#     For gravity=:pointwise, D21 = g  ->  Schurm = g/fac (constant), so this
#     reproduces Marco exactly.
#   * phi_cache gathered in the fused front (needed by the split gravity).
#
# =====================================================================

using MuladdMacro
using Trixi: get_node_vars, set_node_vars!, get_node_coords, @threaded

mutable struct JacobianCacheMarcoSplitNS{NZ, M, RealT, AType2 <: AbstractArray,
                             AType3 <: AbstractArray, AType4 <: AbstractArray,
                             zColType <: AbstractArray, JacobianAuxType,
                             RotationMatrixType, IDMapType, IZMapType, UVertType}
    M::Int
    nz::Int
    nx::Int
    fac::RealT
    schur_dinv :: RealT      # = inv(fac) = γΔt  (D11 = fac always; scalar)
    gravity::Symbol          # :pointwise or :split
    A12_11::AType2
    A22_diag::AType3
    A33_diag::AType3
    FacGrav::RealT
    DWS::AType2
    Th_cache::AType3
    dpdRhoTh_cache::AType3
    phi_cache::AType3        # (M, nz, NumG) geopotential (frozen) — for split gravity
    Gblk::AType4             # (M, M, nz, NumG) gravity block A31 = G[i,k] (constant)
    grav_done::Base.RefValue{Bool}   # lazy precompute flag for Gblk
    SA::AType4
    SchurD  :: AType4
    SchurL  :: AType4
    SchurU  :: AType4
    Schurm      :: AType2    # (nz, NumG) per-element ρ-elimination multiplier D21/fac
    SchurLhat   :: AType4
    SchurDinv   :: AType4
    SchurUhat   :: AType4
    rs      :: AType2
    NumG::Int
    dz::zColType
    cS::RealT
    invwB::RealT
    jacobian_aux::JacobianAuxType
    wPos::Int
    ThPos::Int
    rotation_matrix::RotationMatrixType
    ID_map::IDMapType
    iz_map::IZMapType
    u_vert::UVertType
end

@inline nvelem(::JacobianCacheMarcoSplitNS{NZ, M}) where {NZ, M} = NZ
@inline nvnodes(::JacobianCacheMarcoSplitNS{NZ, M}) where {NZ, M} = M

struct JacobianAuxMarcoNS{R2Type, ThType, SALType}
    r1::R2Type
    r2::R2Type
    r3::R2Type
    Th::ThType
    dpdRhoTh::ThType
    SAL::SALType
end
function JacobianAuxMarcoNS(M, FT)
    r1  = zeros(FT, M - 1, 3)   # 3 driver columns (ϱ_M, ϱw_M, ϱθ_M)
    r2  = zeros(FT, M - 1, 3)
    r3  = zeros(FT, M - 1, 3)
    Th  = zeros(FT, M)
    dpdRhoTh = zeros(FT, M)
    SAL = zeros(FT, M - 1, M - 1)
    return JacobianAuxMarcoNS{typeof(r2), typeof(Th), typeof(SAL)}(r1, r2, r3, Th, dpdRhoTh, SAL)
end

function JacobianCacheMarcoSplitNS(M, nz, nx, NumG, semi, cS, dt,
                             equations::Union{CompressibleEulerPotentialTemperatureEquationsWithGravity2D,
                             CompressibleEulerPotentialTemperatureEquationsWithGravity3D};
                             gravity::Symbol = :split)
    FT = Float64
    polydeg = M - 1
    fac = inv(dt)
    FacGrav = semi.equations.g
    schur_dinv = inv(fac)
    A12_11 = zeros(FT, nz, NumG)
    A22_diag = zeros(FT, polydeg, nz, NumG)
    A33_diag = zeros(FT, polydeg, nz, NumG)
    A22_diag .= inv(fac)
    A33_diag .= fac
    DWS             = copy(semi.solver.basis.derivative_hat)
    Th_cache        = zeros(FT, M, nz, NumG)
    dpdRhoTh_cache  = zeros(FT, M, nz, NumG)
    phi_cache       = zeros(FT, M, nz, NumG)
    Gblk            = zeros(FT, M, M, nz, NumG)
    grav_done       = Ref(false)
    SA = zeros(FT, polydeg, polydeg, nz, NumG)
    SchurD = zeros(FT, 3, 3, nz, NumG)
    # Full 3×3 block-tridiagonal (the split gravity makes ϱ_M couple densely,
    # so the pointwise 2×2 reduction is invalid — see memory note).
    SchurL = zeros(FT, 3, 3, nz - 1, NumG)
    SchurU = zeros(FT, 3, 3, nz - 1, NumG)
    Schurm = zeros(FT, nz, NumG)             # unused by the 3×3 Thomas (kept for struct compat)
    Schurm .= FacGrav / fac
    SchurLhat = zeros(FT, 3, 3, nz - 1, NumG)
    SchurDinv = zeros(FT, 3, 3, nz, NumG)
    SchurUhat = zeros(FT, 3, 3, nz - 1, NumG)
    rs = zeros(FT, 3 * nz, NumG)

    invwB = inv(semi.solver.basis.weights[1])
    zCol = compute_vertical_step_size(NumG, polydeg, semi, nx, nz, typeof(semi.mesh))
    wPos, ThPos = compute_variable_position(equations)
    jacobian_aux = [JacobianAuxMarcoNS(M, FT) for _ in 1:Threads.nthreads()]
    rotation_matrix = pre_compute_rotation_matrix(semi, equations)
    ID_map, iz_map = build_vertical_mapping(semi, nx, nz, equations)
    u_vert = zeros(FT, M, nz, NumG, nvariables(semi))

    return JacobianCacheMarcoSplitNS{nz, M, FT, typeof(A12_11), typeof(Th_cache),
                         typeof(SA), typeof(zCol), typeof(jacobian_aux),
                         typeof(rotation_matrix), typeof(ID_map), typeof(iz_map),
                         typeof(u_vert)}(M, nz, nx, fac, schur_dinv, gravity,
                                         A12_11, A22_diag, A33_diag, FacGrav,
                                         DWS, Th_cache, dpdRhoTh_cache, phi_cache,
                                         Gblk, grav_done,
                                         SA, SchurD, SchurL, SchurU,
                                         Schurm, SchurLhat, SchurDinv, SchurUhat,
                                         rs, NumG, zCol, cS, invwB,
                                         jacobian_aux, wPos, ThPos,
                                         rotation_matrix, ID_map, iz_map, u_vert)
end

# Lazy precompute of the (state-independent) gravity block G[i,k] per element.
#   :pointwise  -> g·I
#   :split      -> G_ik = inv2dz[ 0.5 D_ik (φ_k-φ_i) + δ_ik 0.5 Σ_m D_im (φ_m-φ_i) ]
# Requires phi_cache (filled by rotate_wrap_to_vertical!).  Run once.
function precompute_gravity!(cache::JacobianCacheMarcoSplitNS)
    cache.grav_done[] && return
    (; M, DWS, phi_cache, Gblk, dz, FacGrav, gravity, NumG) = cache
    nz = nvelem(cache)
    @threaded for ID in 1:NumG
        @inbounds for iz in 1:nz
            inv2dz = 2 / dz[iz, ID]
           # if gravity === :pointwise
           #     for k in 1:M, i in 1:M
           #         Gblk[i, k, iz, ID] = (i == k) ? FacGrav : zero(eltype(Gblk))
           #     end
           # else
                for i in 1:M
                    acc = zero(eltype(Gblk))
                    for k in 1:M
                        dphi = phi_cache[k, iz, ID] - phi_cache[i, iz, ID]
                        Gblk[i, k, iz, ID] = inv2dz * 0.5 * DWS[i, k] * dphi
                        acc += DWS[i, k] * dphi
                    end
                    Gblk[i, i, iz, ID] += inv2dz * 0.5 * acc
                end
            #end
        end
    end
    cache.grav_done[] = true
    return nothing
end

@muladd begin
function FillJacDGVertKernel!(cache_jacobian::JacobianCacheMarcoSplitNS, u, semi, equations)
    (; A12_11, A22_diag, A33_diag, SA, SchurD, SchurL, SchurU, Schurm, M, jacobian_aux,
       fac, FacGrav, DWS, Th_cache, dpdRhoTh_cache, Gblk, wPos, ThPos, invwB) = cache_jacobian
    (; NumG, dz, cS) = cache_jacobian
    kappa   = equations.R / equations.c_p
    invcS   = 1 / cS
    kexp    = kappa / (eltype(u)(1) - kappa)
    kfac    = eltype(u)(1) / (eltype(u)(1) - kappa) * equations.R
    pre_fac = kfac * (equations.R / equations.p_0)^kexp
    RhoPos  = 1
    invfac  = inv(fac)
    nz      = nvelem(cache_jacobian)

    @threaded for ID in 1:NumG
        aux      = jacobian_aux[Threads.threadid()]
        Th       = aux.Th
        dpdRhoTh = aux.dpdRhoTh
        SAL      = aux.SAL

        @inbounds for iz in 1:nz
            inv2dz = eltype(u)(2) / dz[iz, ID]
            invdz  = eltype(u)(1) / dz[iz, ID]

            @inbounds for i in 1:M
                Th[i]       = u[i, iz, ID, ThPos] / u[i, iz, ID, RhoPos]
                dpdRhoTh[i] = pre_fac * u[i, iz, ID, ThPos]^kexp
            end
            @inbounds for i in 1:M
                Th_cache[i, iz, ID]       = Th[i]
                dpdRhoTh_cache[i, iz, ID] = dpdRhoTh[i]
            end

            if iz == 1
                A33_diag[1, iz, ID] = fac + inv2dz * cS * invwB
            else
                A12_11[iz, ID]      = dpdRhoTh[1] * invcS * invdz * invwB
                A22_diag[1, iz, ID] = inv(fac + Th[1] * dpdRhoTh[1] * invcS * invdz * invwB)
                A33_diag[1, iz, ID] = fac + cS * invdz * invwB
            end

            # ── SA = A33 - invfac·(G_II·A13) - A32·A22⁻¹·A23 + invfac·(G_II·A12)·A22⁻¹·A23
            # G is the (constant) gravity block; pointwise G=g·I reduces to Marco.
            is_iz1 = (iz == 1)
            dws11  = DWS[1, 1]

            @inbounds for j in 1:(M - 1)
                thj_inv2dz = Th[j] * inv2dz
                @simd for k in 1:(M - 1)
                    SAL[k, j] = DWS[k, j] * thj_inv2dz * A22_diag[k, iz, ID]
                end
            end
            if !is_iz1
                SAL[1, 1] = zero(eltype(u))
            end

            # SA gravity term: -invfac·(G_II·A13),  A13[l,j]=inv2dz·DWS[l,j], A13[1,1]=0 (iz>1)
            @inbounds for j in 1:(M - 1)
                @simd for i in 1:(M - 1)
                    SA[i, j, iz, ID] = zero(eltype(u))
                end
            end
            @inbounds for j in 1:(M - 1)
                for l in 1:(M - 1)
                    a13lj = inv2dz * DWS[l, j] * invfac
                    @simd for i in 1:(M - 1)
                        SA[i, j, iz, ID] -= Gblk[i, l, iz, ID] * a13lj
                    end
                end
                SA[j, j, iz, ID] += A33_diag[j, iz, ID]
            end
            if !is_iz1
                # A13[1,1]=0 for iz>1: remove the (l=1,j=1) contribution G[i,1]·A13[1,1]
                corr13 = inv2dz * dws11 * invfac
                @simd for i in 1:(M - 1)
                    SA[i, 1, iz, ID] += Gblk[i, 1, iz, ID] * corr13
                end
            end

            # SA -= A32·A22⁻¹·A23
            @inbounds for j in 1:(M - 1)
                @inbounds for k in 1:(M - 1)
                    sal_kj      = SAL[k, j]
                    dpdk_inv2dz = dpdRhoTh[k] * inv2dz
                    @simd for i in 1:(M - 1)
                        SA[i, j, iz, ID] -= DWS[i, k] * dpdk_inv2dz * sal_kj
                    end
                end
            end
            _corr32 = (is_iz1 ? 2 : 1) * dws11 * dpdRhoTh[1] * inv2dz
            @simd for j in 1:(M - 1)
                SA[1, j, iz, ID] += _corr32 * SAL[1, j]
            end
            # A12 cross-term: +invfac·A12_11·G_II[:,1]·SAL[1,:]  (column G[:,1], all rows)
            a12scale = invfac * A12_11[iz, ID]
            @inbounds for j in 1:(M - 1)
                sal1j = SAL[1, j]
                @simd for i in 1:(M - 1)
                    SA[i, j, iz, ID] += a12scale * Gblk[i, 1, iz, ID] * sal1j
                end
            end
            LUFull!(iz, ID, SA, Val(nvnodes(cache_jacobian) - 1))

            # ── D block ──────────────────────────────────────────────────
            D11 = fac
            D21 = Gblk[M, M, iz, ID]      # boundary gravity ρw_M ← ρ_M (=g pointwise)
            D31 = eltype(u)(0)
            if iz < nz
                D12 = eltype(u)(0)
                D22 = fac + cS * invdz * invwB
                D32 = eltype(u)(0)
                D13 = dpdRhoTh[M] * invcS * invdz * invwB
                D23 = eltype(u)(0)
                D33 = fac + Th[M] * dpdRhoTh[M] * invcS * invdz * invwB
            else
                D12 = eltype(u)(2) * DWS[M, M] / dz[iz, ID]
                D22 = fac + inv2dz * cS * invwB
                D32 = inv2dz * DWS[M, M] * Th[M]
                D13 = eltype(u)(0)
                D23 = -inv2dz * DWS[M, M] * dpdRhoTh[M]
                D33 = fac
            end
            SchurD[1, 1, iz, ID] = D11; SchurD[2, 1, iz, ID] = D21; SchurD[3, 1, iz, ID] = D31
            SchurD[1, 2, iz, ID] = D12; SchurD[2, 2, iz, ID] = D22; SchurD[3, 2, iz, ID] = D32
            SchurD[1, 3, iz, ID] = D13; SchurD[2, 3, iz, ID] = D23; SchurD[3, 3, iz, ID] = D33
            Schurm[iz, ID] = D21 * invfac            # per-element ρ-elim multiplier

            if iz < nz
                f0 = eltype(u)(0)
                @simd for c in 1:3
                    SchurL[1, c, iz, ID] = f0; SchurL[2, c, iz, ID] = f0; SchurL[3, c, iz, ID] = f0
                    SchurU[1, c, iz, ID] = f0; SchurU[2, c, iz, ID] = f0; SchurU[3, c, iz, ID] = f0
                end
            end
        end
    end
end

# Generic single-column local interior solve for SchurBoundary.
#   On entry:  R2 = bθ (interior ϱθ rhs), R3 = bw (interior ϱw rhs), BR = bρ (interior ϱ rhs).
#   On exit :  R1 = ϱ response, R2 = ϱθ response, R3 = ϱw response  (= A_ll⁻¹ · B-column).
# This mirrors ldivVerticalFKernel!'s interior solve (gravity matvec + Bgrav-free,
# A32, A12 cross, SA solve, ϱθ recovery, ϱ recovery).  Pointwise (G=g·I) it reduces
# to the validated scalar-gravity form.
# Phase A (one column): transform the momentum rhs R3 *before* the SA solve.
#   On entry per column:  R1 = bρ, R2 = bθ, R3 = bw.  Leaves R1,R2 untouched.
@inline function schur_presa_col!(R1, R2, R3, iz, ID, A22_diag, DWS, dpdRhoTh_cache, Gblk,
                                  a12, inv2dz, invfac, is_iz1, polydeg)
    @inbounds begin
        # R3 -= invfac · G_II · bρ   (matvec; pointwise = -fgi·bρ)
        for l in 1:polydeg
            brl = R1[l] * invfac
            @simd for i in 1:polydeg
                R3[i] -= Gblk[i, l, iz, ID] * brl
            end
        end
        # R3 -= A32 · A22⁻¹ · bθ
        for j in 1:polydeg
            r2js = R2[j] * A22_diag[j, iz, ID]
            dj   = dpdRhoTh_cache[j, iz, ID] * inv2dz
            @simd for i in 1:polydeg
                R3[i] -= DWS[i, j] * dj * r2js
            end
        end
        R3[1] += (is_iz1 ? 2 : 1) * DWS[1, 1] * dpdRhoTh_cache[1, iz, ID] * inv2dz *
                  R2[1] * A22_diag[1, iz, ID]
        # A12 cross:  R3 += invfac·a12·G[:,1]·(A22⁻¹ bθ)_1
        a12c = invfac * a12 * R2[1] * A22_diag[1, iz, ID]
        @simd for i in 1:polydeg
            R3[i] += a12c * Gblk[i, 1, iz, ID]
        end
    end
    return nothing
end

# Phase C (one column): recover ϱθ (R2) and ϱ (R1) after the SA solve has overwritten R3.
#   On entry:  R1 = bρ, R2 = bθ, R3 = SA-solved momentum.  On exit R1,R2 = responses.
@inline function schur_postsa_col!(R1, R2, R3, iz, ID, A22_diag, DWS, Th_cache,
                                   a12, inv2dz, invfac, is_iz1, polydeg)
    @inbounds begin
        # ϱθ recovery:  R2 = A22⁻¹·(bθ - A23·R3)
        for j in 1:polydeg
            r3j = R3[j];  tj = Th_cache[j, iz, ID] * inv2dz
            @simd for i in 1:polydeg
                R2[i] -= DWS[i, j] * tj * r3j
            end
        end
        if !is_iz1
            R2[1] += DWS[1, 1] * Th_cache[1, iz, ID] * inv2dz * R3[1]
        end
        @simd for i in 1:polydeg
            R2[i] *= A22_diag[i, iz, ID]
        end
        # ϱ recovery:  R1 = invfac·(bρ - a12·R2[1]·e1 - A13·R3)   (R1 still holds bρ)
        for j in 1:polydeg
            r3j = R3[j]
            @simd for i in 1:polydeg
                R1[i] -= inv2dz * DWS[i, j] * r3j
            end
        end
        R1[1] -= a12 * R2[1]
        if !is_iz1
            R1[1] += inv2dz * DWS[1, 1] * R3[1]
        end
        @simd for i in 1:polydeg
            R1[i] *= invfac
        end
    end
    return nothing
end

# Batched LU triangular solve of the SA factor over the 3 driver columns of r3.
@inline function ldivFull3!(iz, ID, SA, r3, ::Val{n}) where {n}
    @inbounds begin
        for k in 1:(n - 1)
            for i in (k + 1):n
                sa = SA[i, k, iz, ID]
                r3[i, 1] -= sa * r3[k, 1];  r3[i, 2] -= sa * r3[k, 2];  r3[i, 3] -= sa * r3[k, 3]
            end
        end
        for k in n:-1:1
            sak = inv(SA[k, k, iz, ID])
            r3[k, 1] *= sak;  r3[k, 2] *= sak;  r3[k, 3] *= sak
            for i in 1:(k - 1)
                sa = SA[i, k, iz, ID]
                r3[i, 1] -= sa * r3[k, 1];  r3[i, 2] -= sa * r3[k, 2];  r3[i, 3] -= sa * r3[k, 3]
            end
        end
    end
    return nothing
end

# Same, over the first 2 columns of r3 (PART 2 has only 2 nonzero drivers).
@inline function ldivFull2!(iz, ID, SA, r3, ::Val{n}) where {n}
    @inbounds begin
        for k in 1:(n - 1)
            for i in (k + 1):n
                sa = SA[i, k, iz, ID]
                r3[i, 1] -= sa * r3[k, 1];  r3[i, 2] -= sa * r3[k, 2]
            end
        end
        for k in n:-1:1
            sak = inv(SA[k, k, iz, ID])
            r3[k, 1] *= sak;  r3[k, 2] *= sak
            for i in 1:(k - 1)
                sa = SA[i, k, iz, ID]
                r3[i, 1] -= sa * r3[k, 1];  r3[i, 2] -= sa * r3[k, 2]
            end
        end
    end
    return nothing
end

# Full 3×3 boundary Schur complement.  Three drivers per element (ϱ_M, ϱw_M, ϱθ_M
# → columns 1,2,3), with the split-gravity couplings Bgrav (ϱ_M → interior ϱw) and
# Cgrav (boundary ϱw_M ← interior ϱ).  Pointwise (G=g·I) the ϱ_M column vanishes and
# Cgrav=0, reproducing the previous 2×2-reducible result.
function SchurBoundaryKernel!(cache::JacobianCacheMarcoSplitNS, equations)
    (; A12_11, A22_diag, DWS, Th_cache, dpdRhoTh_cache, dz, Gblk,
       SA, SchurD, SchurL, SchurU, M, NumG, jacobian_aux, fac, cS, invwB) = cache
    FT      = eltype(SA)
    invfac  = inv(fac)
    invcS   = 1 / cS
    nz      = cache.nz
    polydeg = M - 1
    valpd   = Val(nvnodes(cache) - 1)

    @threaded for ID in 1:NumG
        aux = jacobian_aux[Threads.threadid()]
        r1 = aux.r1;  r2 = aux.r2;  r3 = aux.r3      # (polydeg, 3) — one column per driver
        @inbounds for iz in 1:nz
            a12    = A12_11[iz, ID]
            inv2dz = 2 / dz[iz, ID]
            is_iz1 = (iz == 1)
            ThM    = Th_cache[M, iz, ID]
            dpdM   = dpdRhoTh_cache[M, iz, ID]

            # interface (iz-1/iz) C_p coefficients — shared by PART 1 (→U[iz-1]) and PART 2
            if !is_iz1
                invdzm1_invwB = invwB / dz[iz-1, ID]
                dpd1     = dpdRhoTh_cache[1, iz, ID]
                th1      = Th_cache[1, iz, ID]
                thM_prev = Th_cache[M, iz-1, ID]
                c2p1 = -dpd1 * invcS * invdzm1_invwB
                c2p2 =  dpd1 * invdzm1_invwB
                c2p3 = -thM_prev * dpd1 * invcS * invdzm1_invwB
                c3p1 =  invdzm1_invwB
                c3p2 = -cS * invdzm1_invwB
                c3p3 =  th1 * invdzm1_invwB
            end

            # ===== PART 1: drivers = boundary(iz) → D[iz] (self) + U[iz-1] (C_p) =====
            # Build the 3 B-columns:  bρ→r1, bθ→r2, bw→r3.
            @simd for i in 1:polydeg
                r1[i,1]=zero(FT); r1[i,2]=zero(FT); r1[i,3]=zero(FT)
                r2[i,1]=zero(FT); r2[i,2]=zero(FT); r2[i,3]=zero(FT)
                r3[i,1]=zero(FT); r3[i,2]=zero(FT); r3[i,3]=zero(FT)
            end
            @simd for i in 1:polydeg
                bdwsM   = inv2dz * DWS[i, M]
                r3[i,1] = Gblk[i, M, iz, ID]   # col1 ϱ_M  : Bgrav (0 for pointwise interior)
                r1[i,2] = bdwsM                # col2 ϱw_M : A13[:,M]
                r2[i,2] = bdwsM * ThM          # col2 ϱw_M : A23[:,M]
                r3[i,3] = bdwsM * dpdM         # col3 ϱθ_M : A32[:,M]
            end
            # Only col 2 (ϱw_M) needs the pre-SA transform: cols 1,3 have bρ=bθ=0,
            # so schur_presa_col! would leave R3 unchanged (it only touches R3 via bρ,bθ).
            schur_presa_col!(view(r1,:,2), view(r2,:,2), view(r3,:,2), iz, ID,
                             A22_diag, DWS, dpdRhoTh_cache, Gblk, a12, inv2dz, invfac, is_iz1, polydeg)
            ldivFull3!(iz, ID, SA, r3, valpd)
            for col in 1:3
                schur_postsa_col!(view(r1,:,col), view(r2,:,col), view(r3,:,col), iz, ID,
                                  A22_diag, DWS, Th_cache, a12, inv2dz, invfac, is_iz1, polydeg)
                # C_self → SchurD[:,col,iz]
                p1 = zero(FT);  p2 = zero(FT);  p3 = zero(FT);  pg = zero(FT)
                @simd for j in 1:polydeg
                    c1j = inv2dz * DWS[M, j]
                    p1 += c1j * r3[j,col]
                    p3 += c1j * Th_cache[j, iz, ID] * r3[j,col]
                    p2 += c1j * dpdRhoTh_cache[j, iz, ID] * r2[j,col]
                    pg += Gblk[M, j, iz, ID] * r1[j,col]
                end
                SchurD[1, col, iz, ID] -= p1
                SchurD[2, col, iz, ID] -= p2 + pg          # Cgrav (ϱw_M ← interior ϱ)
                SchurD[3, col, iz, ID] -= p3
                # C_p → SchurU[:,col,iz-1]  (node-1 response)
                if !is_iz1
                    r2_1 = r2[1,col];  r3_1 = r3[1,col]
                    SchurU[1, col, iz-1, ID] -= c2p1*r2_1 + c3p1*r3_1
                    SchurU[2, col, iz-1, ID] -= c2p2*r2_1 + c3p2*r3_1
                    SchurU[3, col, iz-1, ID] -= c2p3*r2_1 + c3p3*r3_1
                end
            end

            # ===== PART 2 (iz>1): drivers = boundary(iz-1) via B_m → D[iz-1] (C_p) + L[iz-1] =====
            # B_m is interface-flux only (ϱ-independent), so the ϱ_M(iz-1) driver (col 1) is a
            # zero column → no SchurL col-1, no D[iz-1] col-1 correction.
            if !is_iz1
                invdz_cur = inv2dz / 2
                dpdM_prev = dpdRhoTh_cache[M, iz-1, ID]
                b1m1 = -invwB * invdz_cur
                b1m2 = -dpdM_prev * invcS * invwB * invdz_cur
                b2m1 = -thM_prev * invdz_cur * invwB
                b2m2 = -th1 * dpdM_prev * invcS * invdz_cur * invwB
                b3m1 = -cS * invdz_cur * invwB
                b3m2 = -dpdM_prev * invdz_cur * invwB
                # Only 2 nonzero drivers: store ϱw_M(iz-1) in slot 1, ϱθ_M(iz-1) in slot 2.
                @simd for i in 1:polydeg
                    r1[i,1]=zero(FT); r1[i,2]=zero(FT)
                    r2[i,1]=zero(FT); r2[i,2]=zero(FT)
                    r3[i,1]=zero(FT); r3[i,2]=zero(FT)
                end
                r1[1,1]=b1m1; r2[1,1]=b2m1; r3[1,1]=b3m1   # slot1 = ϱw_M(iz-1) (suffix 1)
                r1[1,2]=b1m2; r2[1,2]=b2m2; r3[1,2]=b3m2   # slot2 = ϱθ_M(iz-1) (suffix 2)
                schur_presa_col!(view(r1,:,1), view(r2,:,1), view(r3,:,1), iz, ID,
                                 A22_diag, DWS, dpdRhoTh_cache, Gblk, a12, inv2dz, invfac, is_iz1, polydeg)
                schur_presa_col!(view(r1,:,2), view(r2,:,2), view(r3,:,2), iz, ID,
                                 A22_diag, DWS, dpdRhoTh_cache, Gblk, a12, inv2dz, invfac, is_iz1, polydeg)
                ldivFull2!(iz, ID, SA, r3, valpd)
                for (slot, col) in ((1, 2), (2, 3))   # slot1→col2 (ϱw_M), slot2→col3 (ϱθ_M)
                    schur_postsa_col!(view(r1,:,slot), view(r2,:,slot), view(r3,:,slot), iz, ID,
                                      A22_diag, DWS, Th_cache, a12, inv2dz, invfac, is_iz1, polydeg)
                    # C_p → SchurD[:,col,iz-1] self correction (node-1 response)
                    r2_1 = r2[1,slot];  r3_1 = r3[1,slot]
                    SchurD[1, col, iz-1, ID] -= c2p1*r2_1 + c3p1*r3_1
                    SchurD[2, col, iz-1, ID] -= c2p2*r2_1 + c3p2*r3_1
                    SchurD[3, col, iz-1, ID] -= c2p3*r2_1 + c3p3*r3_1
                    # C_self → SchurL[:,col,iz-1]  (full interior(iz) response + Cgrav)
                    p1 = zero(FT);  p2 = zero(FT);  p3 = zero(FT);  pg = zero(FT)
                    @simd for j in 1:polydeg
                        c1j = inv2dz * DWS[M, j]
                        p1 += c1j * r3[j,slot]
                        p3 += c1j * Th_cache[j, iz, ID] * r3[j,slot]
                        p2 += c1j * dpdRhoTh_cache[j, iz, ID] * r2[j,slot]
                        pg += Gblk[M, j, iz, ID] * r1[j,slot]
                    end
                    SchurL[1, col, iz-1, ID] -= p1
                    SchurL[2, col, iz-1, ID] -= p2 + pg     # Cgrav
                    SchurL[3, col, iz-1, ID] -= p3
                end
            end
        end
    end
end
end # @muladd

# In-place 3×3 inverse (cofactors) written into Ainv[:,:,i,ID].
@inline function invert_3x3!(Ainv, a11,a12,a13, a21,a22,a23, a31,a32,a33, i, ID)
    c11 =  (a22*a33 - a23*a32)
    c12 = -(a21*a33 - a23*a31)
    c13 =  (a21*a32 - a22*a31)
    det = a11*c11 + a12*c12 + a13*c13
    invdet = inv(det)
    @inbounds begin
        Ainv[1,1,i,ID] =  c11 * invdet
        Ainv[2,1,i,ID] =  c12 * invdet
        Ainv[3,1,i,ID] =  c13 * invdet
        Ainv[1,2,i,ID] = -(a12*a33 - a13*a32) * invdet
        Ainv[2,2,i,ID] =  (a11*a33 - a13*a31) * invdet
        Ainv[3,2,i,ID] = -(a11*a32 - a12*a31) * invdet
        Ainv[1,3,i,ID] =  (a12*a23 - a13*a22) * invdet
        Ainv[2,3,i,ID] = -(a11*a23 - a13*a21) * invdet
        Ainv[3,3,i,ID] =  (a11*a22 - a12*a21) * invdet
    end
    return nothing
end

# Full 3×3 block-tridiagonal LU (Thomas).  Diagonal D, super-diag U, sub-diag L
# are stored full 3×3 (rows = eqs ϱ_M,ϱw_M,ϱθ_M; cols = vars in the same order).
# Lhat[i-1] = L[i-1]·Dinv[i-1];  D[i] -= Lhat[i-1]·U[i-1];  Dinv[i] = inv(D[i]).
function block_lu!(F::JacobianCacheMarcoSplitNS)
    nz   = nvelem(F)
    NumG = F.NumG
    D = F.SchurD;  U = F.SchurU;  L = F.SchurL;  Lh = F.SchurLhat;  Di = F.SchurDinv
    @threaded for ID in 1:NumG
        @inbounds begin
            invert_3x3!(Di, D[1,1,1,ID],D[1,2,1,ID],D[1,3,1,ID],
                            D[2,1,1,ID],D[2,2,1,ID],D[2,3,1,ID],
                            D[3,1,1,ID],D[3,2,1,ID],D[3,3,1,ID], 1, ID)
            for i in 2:nz
                # Lhat = L[i-1] · Dinv[i-1]
                for c in 1:3
                    a1 = Di[1,c,i-1,ID];  a2 = Di[2,c,i-1,ID];  a3 = Di[3,c,i-1,ID]
                    Lh[1,c,i-1,ID] = L[1,1,i-1,ID]*a1 + L[1,2,i-1,ID]*a2 + L[1,3,i-1,ID]*a3
                    Lh[2,c,i-1,ID] = L[2,1,i-1,ID]*a1 + L[2,2,i-1,ID]*a2 + L[2,3,i-1,ID]*a3
                    Lh[3,c,i-1,ID] = L[3,1,i-1,ID]*a1 + L[3,2,i-1,ID]*a2 + L[3,3,i-1,ID]*a3
                end
                # D[i] -= Lhat · U[i-1]
                for c in 1:3
                    u1 = U[1,c,i-1,ID];  u2 = U[2,c,i-1,ID];  u3 = U[3,c,i-1,ID]
                    D[1,c,i,ID] -= Lh[1,1,i-1,ID]*u1 + Lh[1,2,i-1,ID]*u2 + Lh[1,3,i-1,ID]*u3
                    D[2,c,i,ID] -= Lh[2,1,i-1,ID]*u1 + Lh[2,2,i-1,ID]*u2 + Lh[2,3,i-1,ID]*u3
                    D[3,c,i,ID] -= Lh[3,1,i-1,ID]*u1 + Lh[3,2,i-1,ID]*u2 + Lh[3,3,i-1,ID]*u3
                end
                invert_3x3!(Di, D[1,1,i,ID],D[1,2,i,ID],D[1,3,i,ID],
                                D[2,1,i,ID],D[2,2,i,ID],D[2,3,i,ID],
                                D[3,1,i,ID],D[3,2,i,ID],D[3,3,i,ID], i, ID)
            end
        end
    end
    return F
end

function solve!(F::JacobianCacheMarcoSplitNS, b::AbstractMatrix{T}) where T
    nz   = nvelem(F)
    NumG = F.NumG
    U = F.SchurU;  Lh = F.SchurLhat;  Di = F.SchurDinv
    @threaded for ID in 1:NumG
        @inbounds begin
            # forward:  y[i] = rhs[i] - Lhat[i-1]·y[i-1]
            for i in 2:nz
                r = 3*(i-1);  rm = 3*(i-2)
                y1 = b[rm+1,ID];  y2 = b[rm+2,ID];  y3 = b[rm+3,ID]
                b[r+1,ID] -= Lh[1,1,i-1,ID]*y1 + Lh[1,2,i-1,ID]*y2 + Lh[1,3,i-1,ID]*y3
                b[r+2,ID] -= Lh[2,1,i-1,ID]*y1 + Lh[2,2,i-1,ID]*y2 + Lh[2,3,i-1,ID]*y3
                b[r+3,ID] -= Lh[3,1,i-1,ID]*y1 + Lh[3,2,i-1,ID]*y2 + Lh[3,3,i-1,ID]*y3
            end
            # back:  x[nz] = Dinv[nz]·y[nz]
            r = 3*(nz-1)
            y1 = b[r+1,ID];  y2 = b[r+2,ID];  y3 = b[r+3,ID]
            b[r+1,ID] = Di[1,1,nz,ID]*y1 + Di[1,2,nz,ID]*y2 + Di[1,3,nz,ID]*y3
            b[r+2,ID] = Di[2,1,nz,ID]*y1 + Di[2,2,nz,ID]*y2 + Di[2,3,nz,ID]*y3
            b[r+3,ID] = Di[3,1,nz,ID]*y1 + Di[3,2,nz,ID]*y2 + Di[3,3,nz,ID]*y3
            for i in nz-1:-1:1
                r = 3*(i-1);  rp = 3*i
                x1 = b[rp+1,ID];  x2 = b[rp+2,ID];  x3 = b[rp+3,ID]
                t1 = b[r+1,ID] - (U[1,1,i,ID]*x1 + U[1,2,i,ID]*x2 + U[1,3,i,ID]*x3)
                t2 = b[r+2,ID] - (U[2,1,i,ID]*x1 + U[2,2,i,ID]*x2 + U[2,3,i,ID]*x3)
                t3 = b[r+3,ID] - (U[3,1,i,ID]*x1 + U[3,2,i,ID]*x2 + U[3,3,i,ID]*x3)
                b[r+1,ID] = Di[1,1,i,ID]*t1 + Di[1,2,i,ID]*t2 + Di[1,3,i,ID]*t3
                b[r+2,ID] = Di[2,1,i,ID]*t1 + Di[2,2,i,ID]*t2 + Di[2,3,i,ID]*t3
                b[r+3,ID] = Di[3,1,i,ID]*t1 + Di[3,2,i,ID]*t2 + Di[3,3,i,ID]*t3
            end
        end
    end
    return b
end

luBandkernel!(cache::JacobianCacheMarcoSplitNS) = block_lu!(cache)

function ldivVerticalFKernel!(cache_jacobian::JacobianCacheMarcoSplitNS, b, equations)
    (; M, A12_11, A22_diag, DWS, Th_cache, dpdRhoTh_cache, Gblk,
       rs, SA, jacobian_aux, NumG, ThPos, wPos, dz, cS, invwB) = cache_jacobian
    FT      = eltype(SA)
    invfac  = inv(cache_jacobian.fac)
    invcS   = 1 / cS
    FacGrav = equations.g
    nz      = nvelem(cache_jacobian)
    polydeg = M - 1
    RhoPos  = 1
    @threaded for ID in 1:NumG
        aux = jacobian_aux[Threads.threadid()]
        r2  = @view aux.r2[:, 1]
        r3  = @view aux.r3[:, 1]
        @inbounds for iz in 1:nz
            sh     = (iz - 1) * 3
            inv2dz = 2 / dz[iz, ID]
            is_iz1 = (iz == 1)
            rs[sh+1, ID] = b[M, iz, ID, RhoPos]
            rs[sh+2, ID] = b[M, iz, ID, wPos]
            rs[sh+3, ID] = b[M, iz, ID, ThPos]
            @simd for i in 1:polydeg
                r2[i] = b[i, iz, ID, ThPos]
                r3[i] = b[i, iz, ID, wPos]
            end
            # r3 -= invfac · G_II · b_ρ   (matvec; pointwise = -fgi·b_ρ)
            for l in 1:polydeg
                brl = b[l, iz, ID, RhoPos] * invfac
                @simd for i in 1:polydeg
                    r3[i] -= Gblk[i, l, iz, ID] * brl
                end
            end
            for j in 1:polydeg
                r2j_sc      = r2[j] * A22_diag[j, iz, ID]
                dpdj_inv2dz = dpdRhoTh_cache[j, iz, ID] * inv2dz
                @simd for i in 1:polydeg
                    r3[i] -= DWS[i, j] * dpdj_inv2dz * r2j_sc
                end
            end
            r3[1] += (is_iz1 ? 2 : 1) * DWS[1, 1] * dpdRhoTh_cache[1, iz, ID] * inv2dz *
                      r2[1] * A22_diag[1, iz, ID]
            # A12 cross: r3[i] += invfac·A12_11·G_II[i,1]·A22⁻¹[1]·r2[1]  (column; pointwise row 1)
            a12c = invfac * A12_11[iz, ID] * r2[1] * A22_diag[1, iz, ID]
            @simd for i in 1:polydeg
                r3[i] += a12c * Gblk[i, 1, iz, ID]
            end
            ldivFull!(iz, ID, SA, r3, Val(nvnodes(cache_jacobian) - 1))
            for j in 1:polydeg
                r3j        = r3[j]
                thj_inv2dz = Th_cache[j, iz, ID] * inv2dz
                @simd for i in 1:polydeg
                    r2[i] -= DWS[i, j] * thj_inv2dz * r3j
                end
            end
            if !is_iz1
                r2[1] += DWS[1, 1] * Th_cache[1, iz, ID] * inv2dz * r3[1]
            end
            @simd for i in 1:polydeg
                r2[i] *= A22_diag[i, iz, ID]
            end
            rs1 = rs[sh+1, ID];  rs2 = rs[sh+2, ID];  rs3 = rs[sh+3, ID]
            @simd for j in 1:polydeg
                c1j = inv2dz * DWS[M, j]
                c3j = c1j * Th_cache[j, iz, ID]
                c2j = c1j * dpdRhoTh_cache[j, iz, ID]
                rs1 -= c1j * r3[j]
                rs3 -= c3j * r3[j]
                rs2 -= c2j * r2[j]
            end
            # Cgrav: rs[ρw_M] -= Σ_j G[M,j]·ρ_int[j]  (zero for pointwise)
            cg = zero(FT)
            for j in 1:polydeg
                rhoj = b[j, iz, ID, RhoPos]
                if j == 1
                    rhoj -= A12_11[iz, ID] * r2[1]
                end
                acc = zero(FT)
                @simd for l in 1:polydeg
                    acc += DWS[j, l] * r3[l]
                end
                rhoj -= inv2dz * acc
                if !is_iz1 && j == 1
                    rhoj += inv2dz * DWS[1, 1] * r3[1]   # A13[1,1]=0
                end
                cg += Gblk[M, j, iz, ID] * (invfac * rhoj)
            end
            rs2 -= cg
            rs[sh+1, ID] = rs1;  rs[sh+2, ID] = rs2;  rs[sh+3, ID] = rs3
            if iz > 1
                sh_prev = (iz - 2) * 3
                x3_1 = r3[1];  x2_1 = r2[1]
                invdzm1_invwB = invwB / dz[iz-1, ID]
                dpd1     = dpdRhoTh_cache[1, iz, ID]
                th1      = Th_cache[1, iz, ID]
                thM_prev = Th_cache[M, iz-1, ID]
                c2p1 = -dpd1 * invcS * invdzm1_invwB
                c2p2 =  dpd1 * invdzm1_invwB
                c2p3 = -thM_prev * dpd1 * invcS * invdzm1_invwB
                c3p1 =  invdzm1_invwB
                c3p2 = -cS * invdzm1_invwB
                c3p3 =  th1 * invdzm1_invwB
                rs[sh_prev+1, ID] -= c3p1*x3_1 + c2p1*x2_1
                rs[sh_prev+2, ID] -= c3p2*x3_1 + c2p2*x2_1
                rs[sh_prev+3, ID] -= c3p3*x3_1 + c2p3*x2_1
            end
        end
    end
end

ldivVerticalSKernel!(cache::JacobianCacheMarcoSplitNS) = solve!(cache, cache.rs)

function ldivVerticalBKernel!(cache_jacobian::JacobianCacheMarcoSplitNS, b, equations)
    (; A12_11, A22_diag, DWS, Th_cache, dpdRhoTh_cache, Gblk,
       SA, rs, M, NumG, jacobian_aux, FacGrav, ThPos, wPos, dz, cS, invwB) = cache_jacobian
    FT      = eltype(SA)
    invfac  = inv(cache_jacobian.fac)
    invcS   = 1 / cS
    nz      = nvelem(cache_jacobian)
    polydeg = M - 1
    RhoPos  = 1
    @threaded for ID in 1:NumG
        aux = jacobian_aux[Threads.threadid()]
        r2  = @view aux.r2[:, 1]
        r3  = @view aux.r3[:, 1]
        @inbounds for iz in 1:nz
            sh     = (iz - 1) * 3
            fgi    = FacGrav * invfac
            inv2dz = 2 / dz[iz, ID]
            is_iz1 = (iz == 1)
            b[M, iz, ID, RhoPos] = rs[sh+1, ID]
            b[M, iz, ID, wPos]   = rs[sh+2, ID]
            b[M, iz, ID, ThPos]  = rs[sh+3, ID]
            ThM_self  = Th_cache[M, iz, ID]
            dpdM_self = dpdRhoTh_cache[M, iz, ID]
            rsh2 = rs[sh+2, ID];  rsh3 = rs[sh+3, ID]
            @simd for i in 1:polydeg
                bdwsM                  = inv2dz * DWS[i, M]
                r2[i]                  = b[i, iz, ID, ThPos] - bdwsM*ThM_self*rsh2
                r3[i]                  = b[i, iz, ID, wPos]  - bdwsM*dpdM_self*rsh3
                b[i, iz, ID, RhoPos] -= bdwsM*rsh2
            end
            if iz > 1
                sh_prev = (iz - 2) * 3
                rsp2 = rs[sh_prev+2, ID];  rsp3 = rs[sh_prev+3, ID]
                invdz_cur = inv2dz / 2
                dpdM_prev = dpdRhoTh_cache[M, iz-1, ID]
                thM_prev  = Th_cache[M, iz-1, ID]
                th1       = Th_cache[1, iz, ID]
                b1m1 = -invwB * invdz_cur
                b1m2 = -dpdM_prev * invcS * invwB * invdz_cur
                b2m1 = -thM_prev * invdz_cur * invwB
                b2m2 = -th1 * dpdM_prev * invcS * invdz_cur * invwB
                b3m1 = -cS * invdz_cur * invwB
                b3m2 = -dpdM_prev * invdz_cur * invwB
                r2[1] -= b2m1*rsp2 + b2m2*rsp3
                r3[1] -= b3m1*rsp2 + b3m2*rsp3
                b[1, iz, ID, RhoPos] -= b1m1*rsp2 + b1m2*rsp3
            end
            # gravity rhs:  r3 -= invfac * G_II * b_rho_mod   (interior-interior block)
            #  plus Bgrav:  r3 -= G[:,M] * rho_M              (boundary rho_M = rs[sh+1])
            # both reduce to the pointwise  r3[i] -= fgi*b_rho[i]  when Gblk = g*I.
            rhoM = rs[sh+1, ID]
            for l in 1:polydeg
                brl = b[l, iz, ID, RhoPos] * invfac
                @simd for i in 1:polydeg
                    r3[i] -= Gblk[i, l, iz, ID] * brl
                end
            end
            @simd for i in 1:polydeg
                r3[i] -= Gblk[i, M, iz, ID] * rhoM
            end
            for j in 1:polydeg
                r2j_sc      = r2[j] * A22_diag[j, iz, ID]
                dpdj_inv2dz = dpdRhoTh_cache[j, iz, ID] * inv2dz
                @simd for i in 1:polydeg
                    r3[i] -= DWS[i, j] * dpdj_inv2dz * r2j_sc
                end
            end
            r3[1] += (is_iz1 ? 2 : 1) * DWS[1, 1] * dpdRhoTh_cache[1, iz, ID] * inv2dz *
                      r2[1] * A22_diag[1, iz, ID]
            # A12 cross-term:  r3 += invfac*A12_11 * G[:,1] * (A22 r2)_1   (column over all rows)
            a12c = invfac * A12_11[iz, ID] * r2[1] * A22_diag[1, iz, ID]
            @simd for i in 1:polydeg
                r3[i] += a12c * Gblk[i, 1, iz, ID]
            end
            ldivFull!(iz, ID, SA, r3, Val(nvnodes(cache_jacobian) - 1))
            for j in 1:polydeg
                r3j        = r3[j]
                thj_inv2dz = Th_cache[j, iz, ID] * inv2dz
                @simd for i in 1:polydeg
                    r2[i] -= DWS[i, j] * thj_inv2dz * r3j
                end
            end
            if !is_iz1
                r2[1] += DWS[1, 1] * Th_cache[1, iz, ID] * inv2dz * r3[1]
            end
            @simd for i in 1:polydeg
                r2[i] *= A22_diag[i, iz, ID]
            end
            b[1, iz, ID, RhoPos] -= A12_11[iz, ID] * r2[1]
            for j in 1:polydeg
                r3j_inv2dz = r3[j] * inv2dz
                @simd for i in 1:polydeg
                    b[i, iz, ID, RhoPos] -= DWS[i, j] * r3j_inv2dz
                end
            end
            if !is_iz1
                b[1, iz, ID, RhoPos] += DWS[1, 1] * inv2dz * r3[1]
            end
            @simd for i in 1:polydeg
                b[i, iz, ID, RhoPos] *= invfac
                b[i, iz, ID, ThPos]   = r2[i]
                b[i, iz, ID, wPos]    = r3[i]
            end
        end
    end
end

# ── fused permutation kernels (with phi gather) ──
@inline function rotate_wrap_to_vertical!(u_vert, u_wrap, cache::JacobianCacheMarcoSplitNS, semi,
                                          equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    (; solver) = semi
    (; ID_map, iz_map, rotation_matrix, phi_cache, wPos, ThPos) = cache
    # The vertical solve only needs the RADIAL momentum (east/north pass through unchanged).
    # So we project onto the radial unit vector r̂ = column 3 of the rotation matrix
    # (3 of 9 floats, one dot product) and store only (ϱ, ϱw_radial, ϱθ) — no east/north,
    # no full 3×3 rotation. FillJac reads only ϱ (slot 1) and ϱθ (ThPos), so this also
    # serves the state pass in update_jacobian!.
    @threaded for element in eachelement(solver, semi.cache)
        for k in eachnode(solver), j in eachnode(solver), i in eachnode(solver)
            ID  = ID_map[i, j, k, element]
            iz  = iz_map[i, j, k, element]
            r13 = rotation_matrix[1, 3, i, j, k, element]
            r23 = rotation_matrix[2, 3, i, j, k, element]
            r33 = rotation_matrix[3, 3, i, j, k, element]
            m1  = u_wrap[2, i, j, k, element]
            m2  = u_wrap[3, i, j, k, element]
            m3  = u_wrap[4, i, j, k, element]
            u_vert[k, iz, ID, 1]     = u_wrap[1, i, j, k, element]   # ϱ
            u_vert[k, iz, ID, wPos]  = r13*m1 + r23*m2 + r33*m3      # ϱw radial = r̂·ϱv
            u_vert[k, iz, ID, ThPos] = u_wrap[5, i, j, k, element]   # ϱθ
            phi_cache[k, iz, ID]     = u_wrap[6, i, j, k, element]
        end
    end
    return nothing
end

@inline function rotate_wrap_to_vertical!(u_vert, u_wrap, cache::JacobianCacheMarcoSplitNS, semi,
                                          equations::CompressibleEulerPotentialTemperatureEquationsWithGravity2D)
    (; solver) = semi
    (; ID_map, iz_map, phi_cache) = cache
    @threaded for element in eachelement(solver, semi.cache)
        for j in eachnode(solver), i in eachnode(solver)
            ID = ID_map[i, j, element]
            iz = iz_map[i, j, element]
            for var in 1:4
                u_vert[j, iz, ID, var] = u_wrap[var, i, j, element]
            end
            phi_cache[j, iz, ID] = u_wrap[5, i, j, element]
        end
    end
    return nothing
end

@inline function unwrap_scale_rotate!(u_vert, b_wrap, cache::JacobianCacheMarcoSplitNS, semi,
                                      equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D,
                                      rosenbrock_fac)
    (; solver) = semi
    (; ID_map, iz_map, rotation_matrix, wPos, ThPos) = cache
    # Only the radial momentum was solved; east/north are unchanged and merely scaled by fac.
    # Reconstruct the Cartesian momentum without re-rotating east/north, using r̂ (column 3):
    #   ϱv_new = fac·ϱv_orig + (ϱw_new − fac·(r̂·ϱv_orig))·r̂ .
    # ϱv_orig is the original RHS momentum, still in b_wrap (untouched by the solve).
    @threaded for element in eachelement(solver, semi.cache)
        for k in eachnode(solver), j in eachnode(solver), i in eachnode(solver)
            ID  = ID_map[i, j, k, element]
            iz  = iz_map[i, j, k, element]
            r13 = rotation_matrix[1, 3, i, j, k, element]
            r23 = rotation_matrix[2, 3, i, j, k, element]
            r33 = rotation_matrix[3, 3, i, j, k, element]
            m1  = b_wrap[2, i, j, k, element]
            m2  = b_wrap[3, i, j, k, element]
            m3  = b_wrap[4, i, j, k, element]
            c   = u_vert[k, iz, ID, wPos] - rosenbrock_fac * (r13*m1 + r23*m2 + r33*m3)
            b_wrap[1, i, j, k, element] = u_vert[k, iz, ID, 1]
            b_wrap[2, i, j, k, element] = rosenbrock_fac*m1 + c*r13
            b_wrap[3, i, j, k, element] = rosenbrock_fac*m2 + c*r23
            b_wrap[4, i, j, k, element] = rosenbrock_fac*m3 + c*r33
            b_wrap[5, i, j, k, element] = u_vert[k, iz, ID, ThPos]
            # b_wrap[6] (φ) is left unchanged.
        end
    end
    return nothing
end

@inline function unwrap_scale_rotate!(u_vert, b_wrap, cache::JacobianCacheMarcoSplitNS, semi,
                                      equations::CompressibleEulerPotentialTemperatureEquationsWithGravity2D,
                                      rosenbrock_fac)
    (; solver) = semi
    (; ID_map, iz_map) = cache
    @threaded for element in eachelement(solver, semi.cache)
        for j in eachnode(solver), i in eachnode(solver)
            ID = ID_map[i, j, element]
            iz = iz_map[i, j, element]
            b_wrap[1, i, j, element] = u_vert[j, iz, ID, 1]
            b_wrap[2, i, j, element] = u_vert[j, iz, ID, 2] * rosenbrock_fac
            b_wrap[3, i, j, element] = u_vert[j, iz, ID, 3]
            b_wrap[4, i, j, element] = u_vert[j, iz, ID, 4]
        end
    end
    return nothing
end

function update_jacobian!(u, semi, equations, cache_jacobian::JacobianCacheMarcoSplitNS, gamma)
    u_wrap = Trixi.wrap_array(u, semi)
    @trixi_timeit timer() "rotate+wrap" rotate_wrap_to_vertical!(cache_jacobian.u_vert, u_wrap, cache_jacobian, semi, equations)
    @trixi_timeit timer() "precompute gravity" precompute_gravity!(cache_jacobian)
    @trixi_timeit timer() "FillJacDGVertKernel" FillJacDGVertKernel!(cache_jacobian, cache_jacobian.u_vert, semi, equations)
    @trixi_timeit timer() "SchurBoundaryKernel" SchurBoundaryKernel!(cache_jacobian, equations)
    @trixi_timeit timer() "luBandkernel" luBandkernel!(cache_jacobian)
    return nothing
end

function solve_jacobian!(b, dt, gamma, semi, equations, cache_jacobian::JacobianCacheMarcoSplitNS, b_vertical_wrap)
    b_wrap = Trixi.wrap_array(b, semi)
    @trixi_timeit timer() "rotate+wrap" rotate_wrap_to_vertical!(b_vertical_wrap, b_wrap, cache_jacobian, semi, equations)
    @trixi_timeit timer() "ldiv F" ldivVerticalFKernel!(cache_jacobian, b_vertical_wrap, equations)
    @trixi_timeit timer() "ldiv S" ldivVerticalSKernel!(cache_jacobian)
    @trixi_timeit timer() "lidv B" ldivVerticalBKernel!(cache_jacobian, b_vertical_wrap, equations)
    rosenbrock_fac = gamma * dt
    @trixi_timeit timer() "unwrap+scale+rotate" unwrap_scale_rotate!(b_vertical_wrap, b_wrap, cache_jacobian, semi, equations, rosenbrock_fac)
    return nothing
end
