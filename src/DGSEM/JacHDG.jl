using StaticArrays
using KernelAbstractions
using LinearAlgebra

################################################################################
# Vertical HDG Jacobian / solver
#
# Memory layouts
#   Th, dpdRhoTh, FacLoc, InvSurf, Geo : (M, nz, NumI)
#   SA                                 : (nz, NumI) of SMatrix{M,M}  (LU factors)
#   Everything that the tridiagonal (Schur) kernels touch has the column index
#   ID as the FASTEST index so that one-thread-per-column access is coalesced:
#   SchurBand : (NumI, 3, nz-1)   factored tridiagonal  (see LinAlg.jl)
#   Blk       : (NumI, 4, nz)     per-element Schur contributions (Jac! time)
#   The RHS of the face system lives in shared memory of ldivHDGVerticalKernel!
#
# There are no atomics: every element writes its own contributions and the
# tridiagonal kernels gather them (also makes the result deterministic).
#
# The same linear algebra serves the conservative variables (rho, rhotheta,
# rhow) and the vector invariant variables (rho, rhotheta, w); only the node
# coefficients differ (VI = true for the vector invariant form):
#   Rm  : coefficient of the mass flux        1             | rho
#   Th  : coefficient of the rhotheta flux    theta         | rhotheta
#   dpdRhoTh : pressure derivative            dp/drhotheta  | dPi/drhotheta
#   Ps  : scaling of the pressure term        1             | cp theta
#   gravity                                   split A31     | linearized
#     cp theta grad(Pi), split form in delta theta (see build_GP)
#   Cs  : sound speed in the interface dissipation, constant cS | local sqrt(gamma p / rho),
#         face value 1/2 (c_L + c_R) as in the interface flux
################################################################################

mutable struct JacHDG{FT<:AbstractFloat,
                         AT3<:AbstractArray,
                         AT_SA<:AbstractArray,
                         SM<:StaticMatrix}
  CompTri::Bool                        
  grav_do::Bool
  NumI::Int
  M::Int
  nz::Int
  fac::FT
  FacGrav::FT
  cS::FT
  wB::FT                 # first 1D quadrature weight (w[1])

  DWS::SM                # DG.V.DWZ  as static matrix (passed by value to kernels)
  DWSS::SM               # DG.V.DWZM as static matrix

  Geo::AT3               # geopotential copy; A31 is rebuilt from it on the fly
  Th::AT3
  dpdRhoTh::AT3
  Rm::AT3
  Ps::AT3
  Pi::AT3                # Exner pressure (vector invariant form)
  Cs::AT3                # node sound speed (vector invariant form)
  VI::Bool
  FacLoc::AT3            # fac * Surf / J
  InvSurf::AT3           # 1 / Surf

  SA::AT_SA              # LU factors (inverse diagonal convention)
  SchurBand::AT3
  Blk::AT3
end  

function JacHDG(backend, FT, M, nz, DG) 
  CompTri = false
  grav_do = true
  fac = zero(FT)
  FacGrav = zero(FT)
  cS = zero(FT)
  NumI = DG.NumI

  # Type alias for M x M static matrix
  MatM = SMatrix{M, M, FT, M*M}

  # Differentiation matrices and first weight as kernel-argument constants
  DWS  = MatM(Tuple(FT.(vec(Array(DG.V.DWZ)[1:M, 1:M]))))
  DWSS = MatM(Tuple(FT.(vec(Array(DG.V.DWZM)[1:M, 1:M]))))
  wB   = FT(Array(DG.V.wZ)[1])

  SA = KernelAbstractions.zeros(backend, MatM, nz, NumI)

  Geo      = KernelAbstractions.zeros(backend, FT, M, nz, NumI)
  Th       = KernelAbstractions.zeros(backend, FT, M, nz, NumI)
  dpdRhoTh = KernelAbstractions.zeros(backend, FT, M, nz, NumI)
  Rm       = KernelAbstractions.ones(backend, FT, M, nz, NumI)
  Ps       = KernelAbstractions.ones(backend, FT, M, nz, NumI)
  Pi       = KernelAbstractions.zeros(backend, FT, M, nz, NumI)
  Cs       = KernelAbstractions.zeros(backend, FT, M, nz, NumI)
  FacLoc   = KernelAbstractions.zeros(backend, FT, M, nz, NumI)
  InvSurf  = KernelAbstractions.zeros(backend, FT, M, nz, NumI)

  SchurBand = KernelAbstractions.zeros(backend, FT, NumI, 3, nz-1)
  Blk       = KernelAbstractions.zeros(backend, FT, NumI, 4, nz)

  return JacHDG{FT,
                   typeof(Th),
                   typeof(SA),
                   typeof(DWS)}(
    CompTri,
    grav_do,
    NumI,
    M,
    nz,
    fac,
    FacGrav,
    cS,
    wB,
    DWS,
    DWSS,
    Geo,
    Th,
    dpdRhoTh,
    Rm,
    Ps,
    Pi,
    Cs,
    false,
    FacLoc,
    InvSurf,
    SA,
    SchurBand,
    Blk,
  )
end

# Gravity matrix A31 for one element, rebuilt from the geopotential
# (M^2 flops instead of loading M^2 values from global memory), with static
# vectors and broadcasts (loops over an MMatrix stay in local memory).
@inline function build_A31(Geo, DS, iz, ID, ::Val{M}) where {M}
  FT = eltype(Geo)
  G = SVector{M, FT}(ntuple(i -> Geo[i, iz, ID], Val(M)))
  DPhi = DS .* (transpose(G) .- G)
  acc = DPhi * SVector{M, FT}(ntuple(i -> one(FT), Val(M)))
  return FT(0.5) * (DPhi + diagm(acc))
end

# Sound speed of the lower (node 1) and upper (node M) face of element iz
@inline function FaceSoundSpeed(Cs, cS, iz, ID, nz, ::Val{M}, ::Val{false}) where {M}
  return cS, cS
end

@inline function FaceSoundSpeed(Cs, cS, iz, ID, nz, ::Val{M}, ::Val{true}) where {M}
  FT = eltype(Cs)
  c1 = iz > 1  ? FT(0.5) * (Cs[1, iz, ID] + Cs[M, iz-1, ID]) : Cs[1, iz, ID]
  cM = iz < nz ? FT(0.5) * (Cs[M, iz, ID] + Cs[1, iz+1, ID]) : Cs[M, iz, ID]
  return c1, cM
end

# Gravity matrix GM (acts on rho) and pressure matrix PrM (acts on rhotheta)
# of the vertical momentum row
@inline function build_GP(Geo, Pi, Rm, Ps, Th, dpdRhoTh, DWS, DWSS, iz, ID, ::Val{M}, ::Val{false}) where {M}
  FT = eltype(Geo)
  GM = build_A31(Geo, DWS, iz, ID, Val(M))
  dpV = SVector{M, FT}(ntuple(i -> dpdRhoTh[i, iz, ID], Val(M)))
  PrM = DWSS .* transpose(dpV)
  return GM, PrM
end

# Vector invariant: the pressure term sum_k D_ik cp theta^num_ik (Pi_k - Pi_i) is
# linearized with delta theta^num = (delta theta_i + delta theta_k) / 2, which gives
# the split form Bs_ik = cp (1/2 D_ik (Pi_k - Pi_i) + delta_ik 1/2 sum_m D_im (Pi_m - Pi_i))
# acting on delta theta_k = (delta rhotheta_k - theta_k delta rho_k) / rho_k
# Built with static vectors and broadcasts: the nested loops over an MMatrix were not
# fully unrolled by LLVM, GM, PrM, P and DWSS then stayed in local memory.
@inline function build_GP(Geo, Pi, Rm, Ps, Th, dpdRhoTh, DWS, DWSS, iz, ID, ::Val{M}, ::Val{true}) where {M}
  FT = eltype(Geo)
  P   = SVector{M, FT}(ntuple(i -> Pi[i, iz, ID], Val(M)))
  Rho = SVector{M, FT}(ntuple(i -> Rm[i, iz, ID], Val(M)))
  ThV = SVector{M, FT}(ntuple(i -> Th[i, iz, ID], Val(M)))
  PsV = SVector{M, FT}(ntuple(i -> Ps[i, iz, ID], Val(M)))
  dpV = SVector{M, FT}(ntuple(i -> dpdRhoTh[i, iz, ID], Val(M)))
  rRho = one(FT) ./ Rho
  # cp = Ps_i / theta_i
  cp = PsV .* Rho ./ ThV
  # DP_ik = D_ik (Pi_k - Pi_i), Bs = cp (1/2 DP + 1/2 diag(sum_m DP_im))
  DP = DWSS .* (transpose(P) .- P)
  acc = DP * SVector{M, FT}(ntuple(i -> one(FT), Val(M)))
  Bs = cp .* (FT(0.5) * (DP + diagm(acc)))
  GM = -Bs .* transpose(ThV .* rRho .* rRho)
  PrM = PsV .* DWSS .* transpose(dpV) + Bs .* transpose(rRho)
  return GM, PrM
end

################################################################################
# Jac! : build and factor the local matrices, per-element Schur contributions
################################################################################
@kernel inbounds = true function FillJacHDGKernel!(
    @Const(Th), @Const(dpdRhoTh), @Const(Rm), @Const(Ps), @Const(Pi), @Const(Cs), FacLoc, InvSurf, SA_gpu, Blk,
    @Const(Geo), @Const(J), @Const(Surf),
    DWS, DWSS, wB, fac, cS, ::Val{M}, ::Val{VI}
) where {M, VI}

  iz, ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  ND = @uniform @ndrange()[2]

  if ID <= ND
    FT = eltype(Th)
    RhoPos = 1
    ThPos = 5
    invwB = inv(wB)
    has_below = iz > 1
    has_above = iz < nz

    # Node-local quantities
    ThL   = MVector{M, FT}(undef)
    dpL   = MVector{M, FT}(undef)
    facL  = MVector{M, FT}(undef)

    @unroll for i = 1:M
      Sf  = Surf[i, iz, ID]
      ThL[i]  = Th[i, iz, ID]
      dpL[i]  = dpdRhoTh[i, iz, ID]
      facL[i] = fac * Sf / J[i, iz, ID]
      FacLoc[i, iz, ID]   = facL[i]
      InvSurf[i, iz, ID]  = inv(Sf)
    end

    RmL = MVector{M, FT}(undef)
    @unroll for i = 1:M
      RmL[i] = Rm[i, iz, ID]
    end

    # Gravity and pressure matrices (not stored)
    GM, PrM = build_GP(Geo, Pi, Rm, Ps, Th, dpdRhoTh, DWS, DWSS, iz, ID, Val(M), Val(VI))

    # Local matrix SA (in registers)
    SA_loc = MMatrix{M, M, FT, M*M}(undef)
    @unroll for i = 1:M
      @unroll for j = 1:M
        val = zero(FT)
        @unroll for k = 1:M
          val = muladd(-(PrM[i,k] * ThL[j] + GM[i,k] * RmL[j]), DWS[k,j], val)
        end
        SA_loc[i,j] = val * facL[i]
      end
      SA_loc[i,i] += inv(facL[i])
    end
    c1, cM = FaceSoundSpeed(Cs, cS, iz, ID, nz, Val(M), Val(VI))
    SA_loc[1,1] += c1 * invwB
    SA_loc[M,M] += cM * invwB

    # LU in registers (inverse-diagonal convention), then store for Solve!
    LUInv!(SA_loc, Val(M))
    SA_gpu[iz, ID] = SA_loc

    # Neighbouring coefficients (only used when the neighbour exists)
    Thm = has_below ? Th[M, iz-1, ID] : ThL[1]
    Thp = has_above ? Th[1, iz+1, ID] : ThL[M]
    Rmm = has_below ? Rm[M, iz-1, ID] : RmL[1]
    Rmp = has_above ? Rm[1, iz+1, ID] : RmL[M]

    # Trace couplings: lower face (node 1) and upper face (node M)
    B1_1 = -FT(0.5) * (RmL[1] + Rmm) * invwB
    B2_1 = -FT(0.5) * (ThL[1] + Thm) * invwB
    B3_1 = -invwB * c1
    B1_2 =  FT(0.5) * (RmL[M] + Rmp) * invwB
    B2_2 =  FT(0.5) * (ThL[M] + Thp) * invwB
    B3_2 = -invwB * cM
    C2_1 =  Ps[1, iz, ID] * dpL[1]
    C2_2 = -Ps[M, iz, ID] * dpL[M]
    C3_1 = -c1
    C3_2 = -cM

    # RHS for the lower-face column (ra) and the upper-face column (rb);
    # both are solved with the same factors in one sweep.
    ra = MVector{M, FT}(undef)
    rb = MVector{M, FT}(undef)
    @unroll for i = 1:M
      ra[i] = -(GM[i,1] * B1_1 + PrM[i,1] * B2_1) * facL[i]
      rb[i] = -(GM[i,M] * B1_2 + PrM[i,M] * B2_2) * facL[i]
    end
    ra[1] += B3_1
    rb[M] += B3_2
    ldivInv2!(SA_loc, ra, rb, Val(M))

    r2_1a = B2_1
    r2_Ma = zero(FT)
    r2_1b = zero(FT)
    r2_Mb = B2_2
    @unroll for k = 1:M
      s1 = DWS[1,k] * ThL[k]
      sM = DWS[M,k] * ThL[k]
      r2_1a -= s1 * ra[k]
      r2_Ma -= sM * ra[k]
      r2_1b -= s1 * rb[k]
      r2_Mb -= sM * rb[k]
    end
    r2_1a *= facL[1]
    r2_1b *= facL[1]
    r2_Ma *= facL[M]
    r2_Mb *= facL[M]

    # Per-element contributions to the Schur tridiagonal (values to be ADDED,
    # i.e. minus the former "t"); gathered later by luTriBlkKernel!.
    # VI: the face diagonal c_L + c_R is added by the two elements (cS2 = 0)
    if has_below
      Blk[ID, 1, iz] = -(r2_1a * C2_1 + ra[1] * C3_1)    # diag  of face iz-1
      if VI
        Blk[ID, 1, iz] += c1
      end
    end
    if has_below & has_above
      Blk[ID, 2, iz] = -(r2_Ma * C2_2 + ra[M] * C3_2)    # sub   (iz, iz-1)
      Blk[ID, 3, iz] = -(r2_1b * C2_1 + rb[1] * C3_1)    # super (iz-1, iz)
    end
    if has_above
      Blk[ID, 4, iz] = -(r2_Mb * C2_2 + rb[M] * C3_2)    # diag  of face iz
      if VI
        Blk[ID, 4, iz] += cM
      end
    end
  end
end

################################################################################
# Host side
################################################################################

function FillJacHDG!(Jac::JacHDG,DG,Metric,fac,Phys)
  
  backend = get_backend(Jac.Th)
  FT = eltype(Jac.Th)
  
  M = Jac.M
  nz = Jac.nz
  ND = Jac.NumI
  
  Jac.fac = fac
  Jac.cS = Phys.cS

  KFillJacHDGKernel! = FillJacHDGKernel!(backend)
  KFillJacHDGKernel!(Jac.Th, Jac.dpdRhoTh, Jac.Rm, Jac.Ps, Jac.Pi, Jac.Cs, Jac.FacLoc, Jac.InvSurf, Jac.SA, Jac.Blk,
    Jac.Geo, Metric.J, Metric.VolSurfVV,
    Jac.DWS, Jac.DWSS, Jac.wB, Jac.fac, Jac.cS, Val(M), Val(Jac.VI);
    ndrange=(nz, ND))
end   

# Assemble the Schur tridiagonal from the element blocks and factor it.
function SchurBoundary!(Jac::JacHDG)

  backend = get_backend(Jac.Blk)
  FT = eltype(Jac.Blk)

  N  = Jac.nz - 1
  ND = Jac.NumI

  KluTriBlkKernel! = luTriBlkKernel!(backend, (128,))
  cS2 = Jac.VI ? zero(FT) : FT(2) * Jac.cS
  KluTriBlkKernel!(Jac.SchurBand, Jac.Blk, cS2, N; ndrange=(ND,))
end  

# bout = J^{-1} b   (bout and b may be the same array)
function SolveHDG!(Jac::JacHDG, bout, b, Metric)

  backend = get_backend(b)

  M  = Jac.M
  nz = Jac.nz
  ND = Jac.NumI
  N  = nz - 1

  fac = Jac.fac
  cS  = Jac.cS

  # forward, face solve and backward in one kernel, a workgroup holds NC complete columns
  NC = max(div(128, nz), 1)
  KldivVerticalKernel! = ldivHDGVerticalKernel!(backend, (nz, NC))
  KldivVerticalKernel!(Jac.Geo, Jac.Th, Jac.dpdRhoTh, Jac.Rm, Jac.Ps, Jac.Pi, Jac.Cs, Jac.FacLoc, Jac.InvSurf, Jac.SA,
    Jac.SchurBand, b, bout, Metric.J, Metric.NVV, Jac.DWS, Jac.DWSS, Jac.wB, fac, cS, Val(M), Val(Jac.VI);
    ndrange=(nz, ND))
end

# Solve! of one column in one kernel: forward elimination (element-local), solve of
# the scalar tridiagonal face system, back-substitution. The face contributions and
# the face solution are kept in shared memory, the element values of the forward
# elimination (r1, r2, b3, Th, dpdRhoTh and the matrices GM, PrM) in private memory
# for the back-substitution. b (input) and bout (output) may be the same array: every
# thread reads its own nodes before writing them, so they are not marked @Const.
@kernel inbounds = true function ldivHDGVerticalKernel!(
    @Const(Geo), @Const(Th), @Const(dpdRhoTh), @Const(Rm), @Const(Ps), @Const(Pi), @Const(Cs),
    @Const(FacLoc), @Const(InvSurf), @Const(SA_gpu), @Const(SchurBand),
    b, bout, @Const(J), @Const(NV), DWS, DWSS, wB, fac, cS, ::Val{M}, ::Val{VI}
) where {M, VI}

  iz, iC = @index(Local, NTuple)
  _, ID = @index(Global, NTuple)

  nz = @uniform @groupsize()[1]
  NC = @uniform @groupsize()[2]
  ND = @uniform @ndrange()[2]

  RsCLoc = @localmem eltype(Th) (2, nz, NC)
  rsLoc  = @localmem eltype(Th) (nz, NC)

  r1  = @private eltype(Th) (M,)
  r2  = @private eltype(Th) (M,)
  b3  = @private eltype(Th) (M,)
  ThL = @private eltype(Th) (M,)
  dpL = @private eltype(Th) (M,)
  GMP = @private eltype(Th) (M, M)
  PrP = @private eltype(Th) (M, M)

  # ---- forward elimination of the element ----
  if ID <= ND
    FT = eltype(Th)
    RhoPos = 1; uPos = 2; vPos = 3; wPos = 4; ThPos = 5

    r3 = MVector{M, FT}(undef)
    @unroll for i = 1:M
      Ji = J[i, iz, ID]
      ThL[i] = Th[i, iz, ID]
      dpL[i] = dpdRhoTh[i, iz, ID]
      r1[i] = b[i, iz, ID, RhoPos] * Ji
      r2[i] = b[i, iz, ID, ThPos] * Ji
      b3[i] = NV[1,i,iz,ID] * b[i, iz, ID, uPos] +
              NV[2,i,iz,ID] * b[i, iz, ID, vPos] +
              NV[3,i,iz,ID] * b[i, iz, ID, wPos]
      r3[i] = b3[i] * Ji
    end

    GM, PrM = build_GP(Geo, Pi, Rm, Ps, Th, dpdRhoTh, DWS, DWSS, iz, ID, Val(M), Val(VI))
    @unroll for j = 1:M
      @unroll for i = 1:M
        GMP[i,j] = GM[i,j]
        PrP[i,j] = PrM[i,j]
      end
    end
    @unroll for i = 1:M
      acc = zero(FT)
      @unroll for j = 1:M
        acc = muladd(GMP[i,j], r1[j], acc)
        acc = muladd(PrP[i,j], r2[j], acc)
      end
      r3[i] -= acc * FacLoc[i, iz, ID]
    end

    SA_static = SA_gpu[iz, ID]
    ldivInv!(SA_static, r3, Val(M))

    r21 = r2[1]
    r2M = r2[M]
    @unroll for j = 1:M
      r21 -= DWS[1, j] * ThL[j] * r3[j]
      r2M -= DWS[M, j] * ThL[j] * r3[j]
    end
    r21 *= FacLoc[1, iz, ID]
    r2M *= FacLoc[M, iz, ID]

    # RsCLoc[1,iz] -> face iz-1 (node 1), RsCLoc[2,iz] -> face iz (node M)
    c1, cM = FaceSoundSpeed(Cs, cS, iz, ID, nz, Val(M), Val(VI))
    if iz > 1
      RsCLoc[1, iz, iC] = -(r21 * (Ps[1, iz, ID] * dpL[1]) + r3[1] * (-c1))
    end
    if iz < nz
      RsCLoc[2, iz, iC] = -(r2M * (-Ps[M, iz, ID] * dpL[M]) + r3[M] * (-cM))
    end
  end

  @synchronize

  # ---- face solve of the column with the factors of LUTriBlk!, face values in shared memory ----
  if iz == 1 && ID <= ND && nz > 1
    N = nz - 1
    y = RsCLoc[2, 1, iC] + RsCLoc[1, 2, iC]
    rsLoc[1, iC] = y
    for i in 2:N
      bi = RsCLoc[2, i, iC] + RsCLoc[1, i+1, iC]
      y  = muladd(-SchurBand[ID,3,i-1], y, bi)
      rsLoc[i, iC] = y
    end
    x = y * SchurBand[ID,2,N]
    rsLoc[N, iC] = x
    for i in N-1:-1:1
      x = muladd(-SchurBand[ID,1,i+1], x, rsLoc[i, iC] * SchurBand[ID,2,i])
      rsLoc[i, iC] = x
    end
  end

  @synchronize

  # ---- backward substitution of the element ----
  if ID <= ND
    FT = eltype(Th)
    RhoPos = 1; uPos = 2; vPos = 3; wPos = 4; ThPos = 5
    invwB = inv(wB)

    r3 = MVector{M, FT}(undef)
    @unroll for i = 1:M
      r3[i] = b3[i] * J[i, iz, ID]
    end

    # Face (trace) corrections
    c1, cM = FaceSoundSpeed(Cs, cS, iz, ID, nz, Val(M), Val(VI))
    if iz < nz
      Thp  = Th[1, iz+1, ID]
      B1_2 = FT(0.5) * (Rm[M, iz, ID] + Rm[1, iz+1, ID]) * invwB
      B2_2 = FT(0.5) * (ThL[M] + Thp) * invwB
      B3_2 = -invwB * cM
      rs_val = rsLoc[iz, iC]
      r1[M] -= B1_2 * rs_val
      r2[M] -= B2_2 * rs_val
      r3[M] -= B3_2 * rs_val
    end
    if iz > 1
      Thm  = Th[M, iz-1, ID]
      B1_1 = -FT(0.5) * (Rm[1, iz, ID] + Rm[M, iz-1, ID]) * invwB
      B2_1 = -FT(0.5) * (ThL[1] + Thm) * invwB
      B3_1 = -invwB * c1
      rs_val = rsLoc[iz-1, iC]
      r1[1] -= B1_1 * rs_val
      r2[1] -= B2_1 * rs_val
      r3[1] -= B3_1 * rs_val
    end

    @unroll for i = 1:M
      acc = zero(FT)
      @unroll for j = 1:M
        acc = muladd(GMP[i,j], r1[j], acc)
        acc = muladd(PrP[i,j], r2[j], acc)
      end
      r3[i] -= acc * FacLoc[i, iz, ID]
    end

    SA_static = SA_gpu[iz, ID]
    ldivInv!(SA_static, r3, Val(M))

    @unroll for i = 1:M
      s1 = r1[i]
      s2 = r2[i]
      @unroll for j = 1:M
        s1 -= DWS[i, j] * (Rm[j, iz, ID] * r3[j])
        s2 -= DWS[i, j] * ThL[j] * r3[j]
      end
      fl = FacLoc[i, iz, ID]
      r1[i] = s1 * fl
      r2[i] = s2 * fl
    end

    @unroll for i = 1:M
      isf = InvSurf[i, iz, ID]
      bb  = r3[i] * isf - fac * b3[i]
      bout[i, iz, ID, RhoPos] = r1[i] * isf
      bout[i, iz, ID, ThPos]  = r2[i] * isf
      bout[i, iz, ID, uPos]   = fac * b[i, iz, ID, uPos] + NV[1,i,iz,ID] * bb
      bout[i, iz, ID, vPos]   = fac * b[i, iz, ID, vPos] + NV[2,i,iz,ID] * bb
      bout[i, iz, ID, wPos]   = fac * b[i, iz, ID, wPos] + NV[3,i,iz,ID] * bb
    end
  end
end

function Jac!(U,fac,DG,Metric,Phys,Cache,JCache::JacHDG,Global)
  InitGeo!(JCache,Cache,DG,Global)
  if Global.Model.ModelType == "VectorInvariant"
    JCache.VI = true
    FillCoeffVI!(JCache.Th,JCache.dpdRhoTh,JCache.Rm,JCache.Ps,JCache.Pi,JCache.Cs,U,Phys,Global.Model.SoundSpeedDG)
  else
    FillThetaAlpha!(JCache.Th,JCache.dpdRhoTh,U,Phys)
  end
  FillJacHDG!(JCache,DG,Metric,fac,Phys)
  SchurBoundary!(JCache)
end

# k = J^{-1} v : the copy k = v is fused into the kernels (v is read, k is written)
function Solve!(k,v,Jac::JacHDG,fac,DG::FiniteElements.DGElement,Metric,Global)
  SolveHDG!(Jac, k, v, Metric)
end
