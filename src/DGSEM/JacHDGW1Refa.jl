using StaticArrays
using KernelAbstractions
using LinearAlgebra

mutable struct JacHDGVert{FT<:AbstractFloat,
                         AT2<:AbstractArray,
                         AT3<:AbstractArray,
                         AT_SA<:AbstractArray}
  CompTri::Bool                        
  grav_do::Bool
  NumI::Int
  M::Int
  nz::Int
  fac::FT
  FacGrav::FT
  cS::FT

  # A31 and SA are stored as 2D GPU Arrays of SMatrix: Size (nz, NumI)
  A31::AT_SA
  Th::AT3
  dpdRhoTh::AT3

  D::AT2

  SA::AT_SA

  SchurBand::AT3
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

  # Allocate A31 and SA directly as 2D CuArrays of Static Arrays
  A31 = KernelAbstractions.zeros(backend, MatM, nz, NumI)
  SA  = KernelAbstractions.zeros(backend, MatM, nz, NumI)

  Th = KernelAbstractions.zeros(backend, FT, M, nz, NumI)
  dpdRhoTh = KernelAbstractions.zeros(backend, FT, M, nz, NumI)
  D = KernelAbstractions.zeros(backend, FT, nz-1, NumI)
  SchurBand = KernelAbstractions.zeros(backend, FT, 3, nz-1, NumI)
  rs = KernelAbstractions.zeros(backend, FT, nz-1, NumI)

  return JacHDGVert{FT,
                   typeof(rs),
                   typeof(Th),
                   typeof(SA)}(
    CompTri,
    grav_do,
    NumI,
    M,
    nz,
    fac,
    FacGrav,
    cS,
    A31,
    Th,
    dpdRhoTh,
    D,
    SA,
    SchurBand,
    rs,
  )
end

@inline function update_schur!(SchurBand, C2_val, C3_val, r2_val, r3_val, i, j, ID)
  t = r2_val * C2_val + r3_val * C3_val
  iB = 2 + i - j
  @atomic SchurBand[iB, j, ID] -= t
end

@inline function update_rs!(rs, C2_val, C3_val, r2_val, r3_val, i, ID)
  t = r2_val * C2_val + r3_val * C3_val
  @atomic rs[i, ID] -= t
end

@kernel inbounds = true function FillJacHDGVertKernel!(
    Th, dpdRhoTh, @Const(A31), SA_gpu, SchurBand, @Const(U),
    @Const(dz), @Const(DWS), @Const(DWSS), @Const(w), fac, cS, Phys, ::Val{M}
) where {M}

  iz, ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  ND = @uniform @ndrange()[2]

  if ID <= ND
    FT = eltype(Th)

    # Static vectors for local node evaluations
    ThL = MVector{M, FT}(undef)
    dpdRhoThL = MVector{M, FT}(undef)

    RhoPos = 1
    ThPos = 5
    wB = w[1]
    invwB = FT(1) / wB

    kappa = Phys.kappa
    kexp = kappa / (FT(1) - kappa)
    kfac = FT(1) / (FT(1) - kappa) * Phys.Rd
    Rdp0 = Phys.Rd / Phys.p0
    facLoc = fac / (dz[iz] * FT(0.5))
    invfac = FT(1) / facLoc

    @unroll for i = 1:M
      ThL[i] = U[i, iz, ID, ThPos] / U[i, iz, ID, RhoPos]
      dpdRhoThL[i] = kfac * (Rdp0 * U[i, iz, ID, ThPos])^kexp
      Th[i, iz, ID] = ThL[i]
      dpdRhoTh[i, iz, ID] = dpdRhoThL[i]
    end

    # 1. Load A31 static matrix into register
    A31_loc = A31[iz, ID]

    # 2. Build local matrix SA in mutable thread memory (MMatrix)
    SA_loc = MMatrix{M, M, FT, M*M}(undef)
    @unroll for i = 1:M
      @unroll for j = 1:M
        val = zero(FT)
        @unroll for k = 1:M
          val -= (DWSS[i,k] * dpdRhoThL[k] * ThL[j] + A31_loc[i,k]) * DWS[k,j]
        end
        SA_loc[i,j] = val * facLoc
      end
      SA_loc[i,i] += invfac
    end
    SA_loc[1,1] += cS * invwB
    SA_loc[M,M] += cS * invwB

    # 3. Save static matrix to global memory for Ldiv kernels

    # Compute LU factorization in thread registers
#   SA_lu = lu(SA_static)
    LUFull!(SA_loc, Val(M))
    SA_static = SMatrix(SA_loc)
    SA_gpu[iz, ID] = SA_static

    # 4. Handle Structural Boundary/Trace RHS vectors
    i_m1 = iz - 1
    i_p0 = iz

    # CASE 1: iz == 1
    if iz == 1
      Thp = U[1, iz + 1, ID, ThPos] / U[1, iz + 1, ID, RhoPos]
      B1_2 = invwB
      B2_2 = FT(0.5) * (ThL[M] + Thp) * invwB
      B3_2 = -invwB * cS
      C2_2 = -dpdRhoThL[M]
      C3_2 = -cS

      r1M = B1_2
      r2M = B2_2
      r3_loc = MVector{M, FT}(undef)
      @unroll for i = 1:M
        a32iM = DWSS[i,M] * dpdRhoThL[M]
        r3_loc[i] = -(A31_loc[i,M] * r1M + a32iM * r2M) * facLoc
      end
      r3_loc[M] += B3_2

      # Solved completely in register space!
#     r3_sol = SA_lu \ SVector(r3_loc)
      ldivFull!(SA_static, r3_loc, Val(M))

      @unroll for k = 1:M
        a23Mk = DWS[M,k] * ThL[k]
        r2M -= a23Mk * r3_loc[k]
      end
      r2M *= facLoc
      update_schur!(SchurBand, C2_2, C3_2, r2M, r3_loc[M], i_p0, i_p0, ID)
    end

    # CASE 2: iz > 1 && iz < nz
    if iz > 1 && iz < nz
      Thm = U[M, iz - 1, ID, ThPos] / U[M, iz - 1, ID, RhoPos]
      Thp = U[1, iz + 1, ID, ThPos] / U[1, iz + 1, ID, RhoPos]
      B1_1 = -invwB;        B1_2 = invwB
      B2_1 = -FT(0.5) * (ThL[1] + Thm) * invwB
      B2_2 =  FT(0.5) * (ThL[M] + Thp) * invwB
      B3_1 = -invwB * cS;   B3_2 = -invwB * cS
      C2_1 = dpdRhoThL[1];           C2_2 = -dpdRhoThL[M]
      C3_1 = -cS;                    C3_2 = -cS

      # Column j = iz - 1
      r11 = B1_1; r21 = B2_1
      r3_loc = MVector{M, FT}(undef)
      @unroll for i = 1:M
        a32i1 = DWSS[i,1] * dpdRhoThL[1]
        r3_loc[i] = -(A31_loc[i,1] * r11 + a32i1 * r21) * facLoc
      end
      r3_loc[1] += B3_1
#     r3_sol = SA_lu \ SVector(r3_loc)
      ldivFull!(SA_static, r3_loc, Val(M))

      r2M = zero(FT)
      @unroll for k = 1:M
        r21 -= DWS[1,k] * ThL[k] * r3_loc[k]
        r2M -= DWS[M,k] * ThL[k] * r3_loc[k]
      end
      r21 *= facLoc; r2M *= facLoc
      update_schur!(SchurBand, C2_1, C3_1, r21, r3_loc[1], i_m1, i_m1, ID)
      update_schur!(SchurBand, C2_2, C3_2, r2M, r3_loc[M], i_p0, i_m1, ID)

      # Column j = iz
      r1M = B1_2; r2M = B2_2
      @unroll for i = 1:M
        a32iM = DWSS[i,M] * dpdRhoThL[M]
        r3_loc[i] = -(A31_loc[i,M] * r1M + a32iM * r2M) * facLoc
      end
      r3_loc[M] += B3_2
#     r3_sol = SA_lu \ SVector(r3_loc)
      ldivFull!(SA_static, r3_loc, Val(M))

      r21 = zero(FT)
      @unroll for k = 1:M
        r21 -= DWS[1,k] * ThL[k] * r3_loc[k]
        r2M -= DWS[M,k] * ThL[k] * r3_loc[k]
      end
      r21 *= facLoc; r2M *= facLoc
      update_schur!(SchurBand, C2_1, C3_1, r21, r3_loc[1], i_m1, i_p0, ID)
      update_schur!(SchurBand, C2_2, C3_2, r2M, r3_loc[M], i_p0, i_p0, ID)
    end

    # CASE 3: iz == nz
    if iz == nz
      Thm = U[M, iz - 1, ID, ThPos] / U[M, iz - 1, ID, RhoPos]
      B1_1 = -invwB
      B2_1 = -FT(0.5) * (ThL[1] + Thm) * invwB
      B3_1 = -invwB * cS
      C2_1 = dpdRhoThL[1]
      C3_1 = -cS

      r11 = B1_1; r21 = B2_1
      r3_loc = MVector{M, FT}(undef)
      @unroll for i = 1:M
        a32i1 = DWSS[i,1] * dpdRhoThL[1]
        r3_loc[i] = -(A31_loc[i,1] * r11 + a32i1 * r21) * facLoc
      end
      r3_loc[1] += B3_1
#     r3_sol = SA_lu \ SVector(r3_loc)
      ldivFull!(SA_static, r3_loc, Val(M))

      @unroll for k = 1:M
        r21 -= DWS[1,k] * ThL[k] * r3_loc[k]
      end
      r21 *= facLoc
      update_schur!(SchurBand, C2_1, C3_1, r21, r3_loc[1], i_m1, i_m1, ID)
    end
  end
end

@kernel inbounds = true function ldivHDGVerticalFKernel!(
    @Const(A31_gpu), @Const(Th), @Const(dpdRhoTh), @Const(SA_gpu),
    @Const(b), rs, @Const(dz), @Const(DWS), @Const(DWSS), @Const(w),
    fac, cS, ::Val{M}
) where {M}

  iz, ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  ND = @uniform @ndrange()[2]

  if ID <= ND
    FT = eltype(Th)
    RhoPos = 1; wPos = 4; ThPos = 5

    dz2 = dz[iz, ID] * FT(0.5)
    facLoc = fac / dz2

    # Build local RHS vectors in register memory
    r1 = MVector{M, FT}(undef)
    r2 = MVector{M, FT}(undef)
    r3 = MVector{M, FT}(undef)

    @unroll for i = 1:M
      r1[i] = b[i, iz, ID, RhoPos] * dz2
      r2[i] = b[i, iz, ID, ThPos] * dz2
    end

    A31_loc = A31_gpu[iz, ID]
    @unroll for i = 1:M
      r3[i] = b[i, iz, ID, wPos] * dz2
      @unroll for j = 1:M
        a31ij = A31_loc[i, j]
        a32ij = DWSS[i, j] * dpdRhoTh[j, iz, ID]
        r3[i] -= (a31ij * r1[j] + a32ij * r2[j]) * facLoc
      end
    end

    # Load SA static matrix and solve in registers
    SA_static = SA_gpu[iz, ID]
#   r3_sol = SA_static \ SVector(r3)
    ldivFull!(SA_static, r3, Val(M))

    r21 = r2[1]
    r2M = r2[M]
    @unroll for j = 1:M
      a231j = DWS[1, j] * Th[j, iz, ID]
      a23Mj = DWS[M, j] * Th[j, iz, ID]
      r21 -= a231j * r3[j]
      r2M -= a23Mj * r3[j]
    end
    r21 *= facLoc
    r2M *= facLoc

    if iz < nz
      C2_2 = -dpdRhoTh[M, iz, ID]
      C3_2 = -cS
      update_rs!(rs, C2_2, C3_2, r2M, r3[M], iz, ID)
    end
    if iz > 1
      C2_1 = dpdRhoTh[1, iz, ID]
      C3_1 = -cS
      update_rs!(rs, C2_1, C3_1, r21, r3[1], iz-1, ID)
    end
  end
end

@kernel inbounds = true function ldivHDGVerticalBKernel!(
    @Const(A31_gpu), @Const(Th), @Const(dpdRhoTh), @Const(SA_gpu),
    b, @Const(rs), @Const(dz), @Const(DWS), @Const(DWSS), @Const(w),
    fac, cS, ::Val{M}
) where {M}

  iz, ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  ND = @uniform @ndrange()[2]

  if ID <= ND
    FT = eltype(Th)
    RhoPos = 1; wPos = 4; ThPos = 5
    wB = w[1]
    invwB = FT(1) / wB

    dz2 = dz[iz, ID] * FT(0.5)
    facLoc = fac / dz2
    Thm = iz > 1  ? Th[M, iz-1, ID] : zero(FT)
    Thp = iz < nz ? Th[1, iz+1, ID] : zero(FT)

    r1 = MVector{M, FT}(undef)
    r2 = MVector{M, FT}(undef)
    r3 = MVector{M, FT}(undef)

    @unroll for i = 1:M
      r1[i] = b[i, iz, ID, RhoPos] * dz2
      r2[i] = b[i, iz, ID, ThPos] * dz2
      r3[i] = b[i, iz, ID, wPos] * dz2
    end

    if iz < nz
      B1_2 = invwB
      B2_2 = FT(0.5) * (Th[M, iz, ID] + Thp) * invwB
      B3_2 = -invwB * cS
      rs_val = rs[iz, ID]
      r1[M] -= B1_2 * rs_val
      r2[M] -= B2_2 * rs_val
      r3[M] -= B3_2 * rs_val
    end
    if iz > 1
      B1_1 = -invwB
      B2_1 = -FT(0.5) * (Th[1, iz, ID] + Thm) * invwB
      B3_1 = -invwB * cS
      rs_val = rs[iz-1, ID]
      r1[1] -= B1_1 * rs_val
      r2[1] -= B2_1 * rs_val
      r3[1] -= B3_1 * rs_val
    end

    A31_loc = A31_gpu[iz, ID]
    r3_rhs = MVector{M, FT}(undef)
    @unroll for i = 1:M
      r3i = zero(FT)
      @unroll for j = 1:M
        a31ij = A31_loc[i, j]
        a32ij = DWSS[i, j] * dpdRhoTh[j, iz, ID]
        r3i += (a31ij * r1[j] + a32ij * r2[j])
      end
      r3_rhs[i] = r3[i] - r3i * facLoc
    end

    # Solve in registers
    SA_static = SA_gpu[iz, ID]
#   r3_sol = SA_static \ SVector(r3_rhs)
    ldivFull!(SA_static, r3_rhs, Val(M))

    @unroll for i = 1:M
      @unroll for j = 1:M
        a23ij = DWS[i, j] * Th[j, iz, ID]
        r1[i] = (r1[i] - DWS[i, j] * r3_rhs[j])
        r2[i] = (r2[i] - a23ij * r3_rhs[j])
      end
      r1[i] *= facLoc
      r2[i] *= facLoc
    end

    @unroll for i = 1:M
      b[i, iz, ID, RhoPos] = r1[i]
      b[i, iz, ID, ThPos]  = r2[i]
      b[i, iz, ID, wPos]   = r3_rhs[i]
    end
  end
end

@kernel inbounds = true function precompute_gravityKernel!(
    A31_gpu, @Const(Geo), @Const(DS), @Const(dz), ::Val{M}
) where {M}

  iz, ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  ND = @uniform @ndrange()[2]

  if ID <= ND
    FT = eltype(dz)
    A31_loc = MMatrix{M, M, FT, M*M}(undef)

    @unroll for i = 1:M
      acc = zero(FT)
      Geoi = Geo[i, iz, ID]
      @unroll for k = 1:M
        dphi = Geo[k, iz, ID] - Geoi
        A31_loc[i, k] = FT(0.5) * DS[i, k] * dphi
        acc += DS[i, k] * dphi
      end
      A31_loc[i, i] += FT(0.5) * acc
    end

    # Write back 2D static matrix element
    A31_gpu[iz, ID] = SMatrix(A31_loc)
  end
end

function precompute_gravity!(GeoPot, dz, DWZ, Jac::JacHDGVert, NumberThreadGPU)
  (; A31) = Jac
  backend = get_backend(dz)
  M = Jac.M
  nz = Jac.nz
  NumI = Jac.NumI

  group = (nz, 1)
  ndrange = (nz, NumI)
  Kprecompute_gravityKernel! = precompute_gravityKernel!(backend, group)
  Kprecompute_gravityKernel!(A31, GeoPot, DWZ, dz, Val(M); ndrange=ndrange)
end


function FillJacHDGVert!(Jac::JacHDGVert,U,DG,dz,fac,Phys)
  
  backend = get_backend(U)
  FTB = eltype(U)
  
  M = Jac.M
  nz = Jac.nz
  DoF  = DG.NumI 
  
  Jac.fac = fac
  Jac.cS = Phys.cS
  
  DWZ = DG.DWZ
  DWZM = DG.DWZM
  
  DoFG = 10
  group = (nz, DoFG)
  ndrange = (nz, DoF)
  @. Jac.SchurBand = 0
  @. Jac.SchurBand[2,:,:] = 2.0 * P.cS
  KFillJacHDGVertKernel! = FillJacHDGVertKernel!(backend,group)
  KFillJacHDGVertKernel!(Jac.Th,Jac.dpdRhoTh,Jac.A31,Jac.SA,Jac.SchurBand,U,dz,
  DWZ,DWZM,DG.wZ,fac,Phys.cS,Phys,Val(M);ndrange=ndrange)

end   

function SchurBoundary!(Jac::JacHDGVert)

  backend = get_backend(Jac.rs)
  FTB = eltype(Jac.rs)

  M = Jac.M
  nz = Jac.nz
  ND = Jac.NumI

  NDG = 32
  group = (NDG)
  ndrange = (ND)
  KluBandKernel! = luBandKernel!(backend,group)
  KluBandKernel!(Jac.SchurBand,Val(1),Val(1),Val(nz-1),ndrange=ndrange)
end  


function Solve!(Jac::JacHDGVert,b,DG,Metric)

  backend = get_backend(b)
  FTB = eltype(b)

  M = Jac.M
  nz = Jac.nz
  ND = Jac.NumI

  invfac = FTB(1) / Jac.fac
  fac = Jac.fac
  cS = Jac.cS
  dz = Metric.dz
  DWZ = DG.DWZ
  DWZM = DG.DWZM
  wZ = DG.wZ

  NDG = 32
  group = (nz, NDG)
  ndrange = (nz, ND)
  @. Jac.rs = 0.0
  KldivVerticalFKernel! = ldivHDGVerticalFKernel!(backend,group)
  KldivVerticalFKernel!(Jac.A31,Jac.Th,Jac.dpdRhoTh,
    Jac.SA,b,Jac.rs,dz,DWZ,DWZM,wZ,fac,cS,Val(M);ndrange=ndrange)

  group = (NDG)
  ndrange = (ND)
  KldivVerticalSKernel! = ldivVerticalSKernel!(backend,group)
  KldivVerticalSKernel!(Jac.SchurBand,Jac.rs,Val(1),Val(1),Val(nz-1);ndrange=ndrange)

  group = (nz, NDG)
  ndrange = (nz, ND)
  KldivVerticalBKernel! = ldivHDGVerticalBKernel!(backend,group)
  KldivVerticalBKernel!(Jac.A31,Jac.Th,Jac.dpdRhoTh,
    Jac.SA,b,Jac.rs,dz,DWZ,DWZM,wZ,fac,cS,Val(M);ndrange=ndrange)

end

function Jac!(U,fac,DG,Metric,Phys,Cache,JCache::JacHDGVert,Global,VelForm)
  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU
  if JCache.grav_do
    @views Geo = Cache.Aux[:,:,:,2]
#   @views Geo = Cache[:,:,:,2]
    precompute_gravity!(Geo,Metric.dz,DG.DWZ,JCache,NumberThreadGPU)
    JCache.grav_do = false
  end  
  dz = Metric.dz
  FillJacHDGVert!(JCache,U,DG,dz,fac,Phys)
  SchurBoundary!(JCache)
end

function Solve!(k,v,Jac::JacHDGVert,fac,DG::FiniteElements.DGElement,Metric,Global,VelForm)

  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU
  @. k = v
  @views TendVCart2VSp!(k,DG,Metric,NumberThreadGPU,VelForm)
  Solve!(Jac,k,DG,Metric)
  @views @. k[:,:,:,2:3] *= fac
  @views TendVSp2VCart!(k,DG,Metric,NumberThreadGPU,VelForm)
end
