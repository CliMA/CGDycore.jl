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
#   RsC       : (NumI, 2, nz)     per-element RHS contributions   (Solve! time)
#   rs        : (NumI, nz-1)      face RHS / solution
#
# There are no atomics: every element writes its own contributions and the
# tridiagonal kernels gather them (also makes the result deterministic).
################################################################################

mutable struct JacHDGVert{FT<:AbstractFloat,
                         AT2<:AbstractArray,
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

  DWS::SM                # DG.DWZ  as static matrix (passed by value to kernels)
  DWSS::SM               # DG.DWZM as static matrix

  Geo::AT3               # geopotential copy; A31 is rebuilt from it on the fly
  Th::AT3
  dpdRhoTh::AT3
  FacLoc::AT3            # fac * Surf / J
  InvSurf::AT3           # 1 / Surf

  SA::AT_SA              # LU factors (inverse diagonal convention)
  SchurBand::AT3
  Blk::AT3
  RsC::AT3
  rs::AT2
end  

function JacHDGVert(backend, FT, M, nz, DG) 
  CompTri = false
  grav_do = true
  fac = zero(FT)
  FacGrav = zero(FT)
  cS = zero(FT)
  NumI = DG.NumI

  # Type alias for M x M static matrix
  MatM = SMatrix{M, M, FT, M*M}

  # Differentiation matrices and first weight as kernel-argument constants
  DWS  = MatM(Tuple(FT.(vec(Array(DG.DWZ)[1:M, 1:M]))))
  DWSS = MatM(Tuple(FT.(vec(Array(DG.DWZM)[1:M, 1:M]))))
  wB   = FT(Array(DG.wZ)[1])

  SA = KernelAbstractions.zeros(backend, MatM, nz, NumI)

  Geo      = KernelAbstractions.zeros(backend, FT, M, nz, NumI)
  Th       = KernelAbstractions.zeros(backend, FT, M, nz, NumI)
  dpdRhoTh = KernelAbstractions.zeros(backend, FT, M, nz, NumI)
  FacLoc   = KernelAbstractions.zeros(backend, FT, M, nz, NumI)
  InvSurf  = KernelAbstractions.zeros(backend, FT, M, nz, NumI)

  SchurBand = KernelAbstractions.zeros(backend, FT, NumI, 3, nz-1)
  Blk       = KernelAbstractions.zeros(backend, FT, NumI, 4, nz)
  RsC       = KernelAbstractions.zeros(backend, FT, NumI, 2, nz)
  rs        = KernelAbstractions.zeros(backend, FT, NumI, nz-1)

  return JacHDGVert{FT,
                   typeof(rs),
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
    FacLoc,
    InvSurf,
    SA,
    SchurBand,
    Blk,
    RsC,
    rs,
  )
end

# Gravity matrix A31 for one element, rebuilt from the geopotential
# (M^2 flops instead of loading M^2 values from global memory).
@inline function build_A31(Geo, DS, iz, ID, ::Val{M}) where {M}
  FT = eltype(Geo)
  G = ntuple(i -> Geo[i, iz, ID], Val(M))
  A = MMatrix{M, M, FT, M*M}(undef)
  @unroll for i = 1:M
    acc = zero(FT)
    @unroll for k = 1:M
      dphi = G[k] - G[i]
      A[i, k] = FT(0.5) * DS[i, k] * dphi
      acc = muladd(DS[i, k], dphi, acc)
    end
    A[i, i] += FT(0.5) * acc
  end
#=
  @unroll for i = 1:M
    @. A[i,:] = FT(0)
    for j = 1 : M
      A[i,i] += DS[i,j] * G[j]  
    end   
  end  
=#  
  return SMatrix(A)
end

################################################################################
# Jac! : build and factor the local matrices, per-element Schur contributions
################################################################################
@kernel inbounds = true function FillJacHDGVertKernel!(
    Th, dpdRhoTh, FacLoc, InvSurf, SA_gpu, Blk,
    @Const(Geo), @Const(U), @Const(J), @Const(Surf),
    DWS, DWSS, wB, fac, cS, kexp, kfac, Rdp0, ::Val{M}
) where {M}

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
    InvSfL  = MVector{M, FT}(undef)

    @unroll for i = 1:M
      rth = U[i, iz, ID, ThPos]
      InvSfL[i]  = inv(Surf[i, iz, ID])
      ThL[i]  = rth / U[i, iz, ID, RhoPos]
      dpL[i]  = kfac * (Rdp0 * rth)^kexp
#     facL[i] = fac * Sf / J[i, iz, ID]
      facL[i] = fac / J[i, iz, ID]
      Th[i, iz, ID]       = ThL[i]
      dpdRhoTh[i, iz, ID] = dpL[i]
      FacLoc[i, iz, ID]   = facL[i]
#     InvSurf[i, iz, ID]  = inv(Sf)
      InvSurf[i, iz, ID]  = InvSfL[i]
    end

    # Gravity matrix from the geopotential (not stored)
    A31 = build_A31(Geo, DWS, iz, ID, Val(M))

    # Local matrix SA (in registers)
    SA_loc = MMatrix{M, M, FT, M*M}(undef)
    @unroll for i = 1:M
      @unroll for j = 1:M
        val = zero(FT)
        @unroll for k = 1:M
          val = muladd(-(DWSS[i,k] * dpL[k] * ThL[j] + A31[i,k]), DWS[k,j], val)
        end
        SA_loc[i,j] = val * facL[i]
      end
      SA_loc[i,i] += inv(facL[i]) * InvSfL[i]^2
    end
    SA_loc[1,1] += cS * invwB * InvSfL[1]
    SA_loc[M,M] += cS * invwB * InvSfL[M]

    # LU in registers (inverse-diagonal convention), then store for Solve!
    LUFull!(SA_loc, Val(M))
    SA_gpu[iz, ID] = SA_loc

    # Neighbouring theta (only used when the neighbour exists)
    Thm = has_below ? U[M, iz-1, ID, ThPos] / U[M, iz-1, ID, RhoPos] : ThL[1]
    Thp = has_above ? U[1, iz+1, ID, ThPos] / U[1, iz+1, ID, RhoPos] : ThL[M]

    # Trace couplings: lower face (node 1) and upper face (node M)
    B1_1 = -invwB
    B2_1 = -FT(0.5) * (ThL[1] + Thm) * invwB
    B3_1 = -invwB * cS
    B1_2 =  invwB
    B2_2 =  FT(0.5) * (ThL[M] + Thp) * invwB
    B3_2 = -invwB * cS
    C2_1 =  dpL[1] * Surf[1, iz, ID]
    C2_2 = -dpL[M] * Surf[M, iz, ID]
    C3   = -cS

    # RHS for the lower-face column (ra) and the upper-face column (rb);
    # both are solved with the same factors in one sweep.
    ra = MVector{M, FT}(undef)
    rb = MVector{M, FT}(undef)
    @unroll for i = 1:M
      ra[i] = -(A31[i,1] * B1_1 + DWSS[i,1] * dpL[1] * B2_1) * facL[i]
      rb[i] = -(A31[i,M] * B1_2 + DWSS[i,M] * dpL[M] * B2_2) * facL[i]
    end
    ra[1] += B3_1 * InvSfL[1]
    rb[M] += B3_2 * InvSfL[M]
    ldivFull2!(SA_loc, ra, rb, Val(M))

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
    if has_below
      Blk[ID, 1, iz] = -(r2_1a * C2_1 + ra[1] * C3)    # diag  of face iz-1
    end
    if has_below & has_above
      Blk[ID, 2, iz] = -(r2_Ma * C2_2 + ra[M] * C3)    # sub   (iz, iz-1)
      Blk[ID, 3, iz] = -(r2_1b * C2_1 + rb[1] * C3)    # super (iz-1, iz)
    end
    if has_above
      Blk[ID, 4, iz] = -(r2_Mb * C2_2 + rb[M] * C3)    # diag  of face iz
    end
  end
end

################################################################################
# Solve! : forward elimination (element-local), tridiagonal solve, back-substitution
################################################################################
@kernel inbounds = true function ldivHDGVerticalFKernel!(
    @Const(Geo), @Const(Th), @Const(dpdRhoTh), @Const(FacLoc), @Const(SA_gpu),
    @Const(b), RsC, @Const(J), @Const(InvSurf), @Const(NV), DWS, DWSS, cS, ::Val{M}
) where {M}

  iz, ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  ND = @uniform @ndrange()[2]

  if ID <= ND
    FT = eltype(Th)
    RhoPos = 1; uPos = 2; vPos = 3; wPos = 4; ThPos = 5

    r1  = MVector{M, FT}(undef)
    r2  = MVector{M, FT}(undef)
    r3  = MVector{M, FT}(undef)
    ThL = MVector{M, FT}(undef)
    dpL = MVector{M, FT}(undef)

    @unroll for i = 1:M
      Ji = J[i, iz, ID] 
      ThL[i] = Th[i, iz, ID]
      dpL[i] = dpdRhoTh[i, iz, ID]
      r1[i] = b[i, iz, ID, RhoPos] * Ji 
      r2[i] = b[i, iz, ID, ThPos] * Ji
      r3[i] = (NV[1,i,iz,ID] * b[i, iz, ID, uPos] +
               NV[2,i,iz,ID] * b[i, iz, ID, vPos] +
               NV[3,i,iz,ID] * b[i, iz, ID, wPos]) * Ji * InvSurf[i, iz, ID]
    end

    A31 = build_A31(Geo, DWS, iz, ID, Val(M))
    @unroll for i = 1:M
      acc = zero(FT)
      @unroll for j = 1:M
        acc = muladd(A31[i,j], r1[j], acc)
        acc = muladd(DWSS[i,j] * dpL[j], r2[j], acc)
      end
      r3[i] -= acc * FacLoc[i, iz, ID]
    end

    # Solve in registers with the stored factors
    SA_static = SA_gpu[iz, ID]
    ldivFull!(SA_static, r3, Val(M))

    r21 = r2[1]
    r2M = r2[M]
    @unroll for j = 1:M
      r21 -= DWS[1, j] * ThL[j] * r3[j]
      r2M -= DWS[M, j] * ThL[j] * r3[j]
    end
    r21 *= FacLoc[1, iz, ID]
    r2M *= FacLoc[M, iz, ID]

    # Contributions to the face RHS (gathered by ldivTriBlkKernel!).
    # RsC[ID,1,iz] -> face iz-1 (node 1), RsC[ID,2,iz] -> face iz (node M).
    if iz > 1
      RsC[ID, 1, iz] = -(r21 * dpL[1] / InvSurf[1, iz, ID] + r3[1] * (-cS))
    end
    if iz < nz
      RsC[ID, 2, iz] = -(r2M * (-dpL[M] / InvSurf[M, iz, ID]) + r3[M] * (-cS))
    end
  end
end

# b (input) and bout (output) may be the same array: every thread reads its own
# nodes before writing them, so they are deliberately NOT marked @Const.
@kernel inbounds = true function ldivHDGVerticalBKernel!(
    @Const(Geo), @Const(Th), @Const(dpdRhoTh), @Const(FacLoc), @Const(InvSurf),
    @Const(SA_gpu), b, bout, @Const(rs), @Const(J), @Const(NV),
    DWS, DWSS, wB, fac, cS, ::Val{M}
) where {M}

  iz, ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  ND = @uniform @ndrange()[2]

  if ID <= ND
    FT = eltype(Th)
    RhoPos = 1; uPos = 2; vPos = 3; wPos = 4; ThPos = 5
    invwB = inv(wB)

    r1  = MVector{M, FT}(undef)
    r2  = MVector{M, FT}(undef)
    r3  = MVector{M, FT}(undef)
    b3  = MVector{M, FT}(undef)
    ThL = MVector{M, FT}(undef)
    dpL = MVector{M, FT}(undef)

    @unroll for i = 1:M
      Ji = J[i, iz, ID] 
      ThL[i] = Th[i, iz, ID]
      dpL[i] = dpdRhoTh[i, iz, ID]
      r1[i] = b[i, iz, ID, RhoPos] * Ji
      r2[i] = b[i, iz, ID, ThPos] * Ji
      b3[i] = NV[1,i,iz,ID] * b[i, iz, ID, uPos] +
              NV[2,i,iz,ID] * b[i, iz, ID, vPos] +
              NV[3,i,iz,ID] * b[i, iz, ID, wPos]
      r3[i] = b3[i] * Ji * InvSurf[i, iz, ID]
    end

    # Face (trace) corrections
    if iz < nz
      Thp  = Th[1, iz+1, ID]
      B1_2 = invwB
      B2_2 = FT(0.5) * (ThL[M] + Thp) * invwB
      B3_2 = -invwB * cS
      rs_val = rs[ID, iz]
      r1[M] -= B1_2 * rs_val
      r2[M] -= B2_2 * rs_val
      r3[M] -= B3_2 * rs_val * InvSurf[M, iz, ID]
    end
    if iz > 1
      Thm  = Th[M, iz-1, ID]
      B1_1 = -invwB
      B2_1 = -FT(0.5) * (ThL[1] + Thm) * invwB
      B3_1 = -invwB * cS
      rs_val = rs[ID, iz-1]
      r1[1] -= B1_1 * rs_val
      r2[1] -= B2_1 * rs_val
      r3[1] -= B3_1 * rs_val * InvSurf[1, iz, ID]
    end

    A31 = build_A31(Geo, DWS, iz, ID, Val(M))
    @unroll for i = 1:M
      acc = zero(FT)
      @unroll for j = 1:M
        acc = muladd(A31[i,j], r1[j], acc)
        acc = muladd(DWSS[i,j] * dpL[j], r2[j], acc)
      end
      r3[i] -= acc * FacLoc[i, iz, ID]          # r3 now holds the solve RHS
    end

    SA_static = SA_gpu[iz, ID]
    ldivFull!(SA_static, r3, Val(M))

    @unroll for i = 1:M
      s1 = r1[i]
      s2 = r2[i]
      @unroll for j = 1:M
        s1 -= DWS[i, j] * r3[j]
        s2 -= DWS[i, j] * ThL[j] * r3[j]
      end
      fl = FacLoc[i, iz, ID]
      r1[i] = s1 * fl
      r2[i] = s2 * fl
    end

    @unroll for i = 1:M
      isf = InvSurf[i, iz, ID]
      bb  = r3[i] * isf - fac * b3[i]
      bout[i, iz, ID, RhoPos] = r1[i] 
      bout[i, iz, ID, ThPos]  = r2[i] 
      bout[i, iz, ID, uPos]   = fac * b[i, iz, ID, uPos] + NV[1,i,iz,ID] * bb
      bout[i, iz, ID, vPos]   = fac * b[i, iz, ID, vPos] + NV[2,i,iz,ID] * bb
      bout[i, iz, ID, wPos]   = fac * b[i, iz, ID, wPos] + NV[3,i,iz,ID] * bb
    end
  end
end

################################################################################
# Host side
################################################################################

# The geopotential is static: copy it once, A31 is rebuilt inside the kernels.
# (GeoPot, DWZ, NumberThreadGPU are kept for call compatibility.)
function precompute_gravity!(GeoPot, DWZ, Jac::JacHDGVert, NumberThreadGPU)
  Jac.Geo .= GeoPot
  return nothing
end

function FillJacHDGVert!(Jac::JacHDGVert,U,DG,Metric,fac,Phys)
  
  backend = get_backend(U)
  FT = eltype(Jac.Th)
  
  M = Jac.M
  nz = Jac.nz
  ND = Jac.NumI
  
  Jac.fac = fac
  Jac.cS = Phys.cS

  # Thermodynamic constants as plain scalars of the working precision
  kappa = FT(Phys.kappa)
  kexp  = kappa / (one(FT) - kappa)
  kfac  = FT(Phys.Rd) / (one(FT) - kappa)
  Rdp0  = FT(Phys.Rd) / FT(Phys.p0)

  KFillJacHDGVertKernel! = FillJacHDGVertKernel!(backend)
  KFillJacHDGVertKernel!(Jac.Th, Jac.dpdRhoTh, Jac.FacLoc, Jac.InvSurf, Jac.SA, Jac.Blk,
    Jac.Geo, U, Metric.J, Metric.VolSurfVV,
    Jac.DWS, Jac.DWSS, Jac.wB, Jac.fac, Jac.cS, kexp, kfac, Rdp0, Val(M);
    ndrange=(nz, ND))
end   

# Assemble the Schur tridiagonal from the element blocks and factor it.
function SchurBoundary!(Jac::JacHDGVert)

  backend = get_backend(Jac.rs)
  FT = eltype(Jac.rs)

  N  = Jac.nz - 1
  ND = Jac.NumI

  KluTriBlkKernel! = luTriBlkKernel!(backend, (128,))
  KluTriBlkKernel!(Jac.SchurBand, Jac.Blk, FT(2) * Jac.cS, N; ndrange=(ND,))
end  

# bout = J^{-1} b   (bout and b may be the same array)
function SolveHDGVert!(Jac::JacHDGVert, bout, b, Metric)

  backend = get_backend(b)

  M  = Jac.M
  nz = Jac.nz
  ND = Jac.NumI
  N  = nz - 1

  fac = Jac.fac
  cS  = Jac.cS

  KldivVerticalFKernel! = ldivHDGVerticalFKernel!(backend)
  KldivVerticalFKernel!(Jac.Geo, Jac.Th, Jac.dpdRhoTh, Jac.FacLoc, Jac.SA,
    b, Jac.RsC, Metric.J, Jac.InvSurf, Metric.NVV, Jac.DWS, Jac.DWSS, cS, Val(M);
    ndrange=(nz, ND))

  KldivTriBlkKernel! = ldivTriBlkKernel!(backend, (128,))
  KldivTriBlkKernel!(Jac.SchurBand, Jac.RsC, Jac.rs, N; ndrange=(ND,))

  KldivVerticalBKernel! = ldivHDGVerticalBKernel!(backend)
  KldivVerticalBKernel!(Jac.Geo, Jac.Th, Jac.dpdRhoTh, Jac.FacLoc, Jac.InvSurf, Jac.SA,
    b, bout, Jac.rs, Metric.J, Metric.NVV, Jac.DWS, Jac.DWSS, Jac.wB, fac, cS, Val(M);
    ndrange=(nz, ND))
end

# In-place variant (same signature as before)

function Solve!(Jac::JacHDGVert,b,DG,Metric,NumberThreadGPU)
  SolveHDGVert!(Jac, b, b, Metric)
end

function Jac!(U,fac,DG,Metric,Phys,Cache,JCache::JacHDGVert,Global,VelForm)
  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU
  if JCache.grav_do
    @views Geo = Cache.Aux[:,:,1:DG.NumI,2]
    precompute_gravity!(Geo,DG.DWZ,JCache,NumberThreadGPU)
    JCache.grav_do = false
  end  
  FillJacHDGVert!(JCache,U,DG,Metric,fac,Phys)
  SchurBoundary!(JCache)
end

# k = J^{-1} v : the copy k = v is fused into the kernels (v is read, k is written)
function Solve!(k,v,Jac::JacHDGVert,fac,DG::FiniteElements.DGElement,Metric,Global,VelForm)
  SolveHDGVert!(Jac, k, v, Metric)
end
