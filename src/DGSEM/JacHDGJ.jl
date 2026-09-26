mutable struct JacHDGJVert{FT<:AbstractFloat,
                         AT2<:AbstractArray,
                         AT3<:AbstractArray,
                         AT4<:AbstractArray}
  CompTri::Bool                        
  grav_do::Bool
  NumI::Int
  M::Int
  nz::Int
  fac::FT
  FacGrav::FT
  cS::FT
  A31::AT4

# On-the-fly quantities: Th and dpdRhoTh are the only per-node state kept.
# A13, A23, A32, A33, B1, B2, B3, C2, C3 are all cheap closed-form
# expressions of Th/dpdRhoTh (+ dz, DWZ, wZ, cS) and are recomputed
# wherever they're needed instead of being materialized as arrays.
  Th::AT3
  dpdRhoTh::AT3

  D::AT2

  SA::AT4

  SchurBand::AT3
  rs::AT2
end  

function JacHDGJVert(backend,FT,M,nz,DG) 
  CompTri = false
  grav_do = true
  fac = 0
  FacGrav = 0
  cS = 0
  NumI = DG.NumI
  A31 = KernelAbstractions.zeros(backend,FT,M,M,nz,NumI)
  Th = KernelAbstractions.zeros(backend,FT,M,nz,NumI)
  dpdRhoTh = KernelAbstractions.zeros(backend,FT,M,nz,NumI)
  D = KernelAbstractions.zeros(backend,FT,nz-1,NumI)
  SA = KernelAbstractions.zeros(backend,FT,M,M,nz,NumI)
  SchurBand = KernelAbstractions.zeros(backend,FT,3,nz-1,NumI)
  rs = KernelAbstractions.zeros(backend,FT,nz-1,NumI)

  return JacHDGJVert{FT,
                   typeof(rs),
                   typeof(Th),
                   typeof(A31)}(

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

# ---------------------------------------------------------------------------
# Shared on-the-fly reconstruction helpers.
# These replace what used to be array reads into A13, A23, A32, A33, B1, B2,
# B3, C2, C3. They are pure, @inline scalar functions so they can be called
# from any kernel (Fill or Solve) without needing extra global memory.
# ---------------------------------------------------------------------------


@inline function calcA32(DWS,dpdRhoThj,i,j,::Val{M}) where M
  a = DWS[i,j] * dpdRhoThj
  return (i == 1 && j == 1) || (i == M && j == M) ? -a : a
end


# Helper inline function to handle the Schur complement atomic update cleanly
@inline function update_schur!(SchurBand, C2_val, C3_val, r2_val, r3_val, i, j, ID)
  t = r2_val * C2_val + r3_val * C3_val
  iB = 2 + i - j
  @atomic SchurBand[iB, j, ID] -= t
end

@kernel inbounds = true function FillJacHDGJVertKernel!(Th,dpdRhoTh,@Const(A31),
  SA,SchurBand,@Const(U),@Const(J),@Const(invJ),@Const(Surf),
  @Const(DWS),@Const(w),fac,cS,Phys, ::Val{M}) where {M}

  iz, ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  DoF = @uniform @ndrange()[2]
  ThL = @private eltype(SA) (M,)
  dpdRhoThL = @private eltype(SA) (M,)
  r3 = @private eltype(SA) (M,)

  @uniform RhoPos = 1
  @uniform ThPos = 5
  @uniform invcS = eltype(SA)(1) / cS
  @uniform wB = w[1]
  @uniform invcSwB = eltype(SA)(1) / (cS * wB)
  @uniform invwB = eltype(SA)(1) / wB
  @uniform FTe = eltype(SA)

  if ID <= DoF
    kappa = Phys.kappa
    kexp = kappa / (eltype(SA)(1) - kappa)
    kfac = eltype(SA)(1) / (eltype(SA)(1) - kappa) * Phys.Rd
    Rdp0 = Phys.Rd / Phys.p0
    invfac = FTe(1) / facLoc

    @unroll for i = 1 : M
      ThL[i] = U[i,iz,ID,ThPos] / U[i,iz,ID,RhoPos]
      dpdRhoThL[i] = kfac * (Rdp0 * U[i,iz,ID,ThPos])^kexp
      Th[i,iz,ID] = ThL[i]
      dpdRhoTh[i,iz,ID] = dpdRhoThL[i]
    end

    @unroll for i = 1 : M
      @unroll for j = 1 : M
        val = zero(eltype(SA))
        @unroll for k = 1 : M
          a32ik = calcA32(DWS,dpdRhoThL[k],i,k,Val(M))
          a23kj = DWS[k,j] * ThL[j] 
          val -= a32ik * a23kj + A31[i,k,iz,ID] * DWS[k,j]
        end
        SA[i,j,iz,ID] = val * fac * Surf[k,iz,ID]^2 * invJ[k,iz,ID]
      end
      SA[i,i,iz,ID] += invfac * J[i,iz,ID]
    end
    SA[1,1,iz,ID] += cS * invwB
    SA[M,M,iz,ID] += cS * invwB

    LUFull!(iz,ID,SA, Val(M))
  end

  # Cache common column indices
  i_m1 = iz - 1
  i_p0 = iz
  # CASE 1: iz == 1
  if iz == 1
    Thp = U[1,iz + 1,ID,ThPos] / U[1,iz + 1,ID,RhoPos]
    B1_1 = zero(FTe);             B1_2 = invwB
    B2_1 = zero(FTe);             B2_2 = FTe(0.5) * (ThL[M] + Thp) * invwB
    B3_1 = zero(FTe);             B3_2 = -invwB * cS
    C2_1 = zero(FTe);             C2_2 = -dpdRhoThL[M] * Surf[M,iz,ID]
    C3_1 = zero(FTe);             C3_2 = -cS
      
    # Column j = iz 
    j = i_p0
    r1M = B1_2
    r2M = B2_2
    @unroll for i = 1 : M
      a32iM = calcA32(DWS,dpdRhoThL[M],i,M,Val(M))
      r3[i] = -(A31[i,M,iz,ID] * r1M + a32iM * r2M) * facLoc
    end
    r3[M] += B3_2
    ldivFull!(iz,ID,SA, r3, Val(M))

    @unroll for k = 1 : M
      a23Mk = DWS[M,k] * ThL[k] 
      r2M -= a23Mk * r3[k]
    end
    r2M *= facLoc
    update_schur!(SchurBand, C2_2, C3_2, r2M, r3[M], i_p0, j, ID)
  end

  # CASE 2: iz > 1 && iz < nz
  if iz > 1 && iz < nz
    Thm = U[M,iz - 1,ID,ThPos] / U[M,iz - 1,ID,RhoPos]
    Thp = U[1,iz + 1,ID,ThPos] / U[1,iz + 1,ID,RhoPos]
    B1_1 = -invwB;       B1_2 = invwB
    B2_1 = -FTe(0.5) * (ThL[1] + Thm) * invwB
    B2_2 =  FTe(0.5) * (ThL[M] + Thp) * invwB
    B3_1 = -invwB * cS;  B3_2 = -invwB * cS
    C2_1 = dpdRhoThL[1] * Surf[1,iz,ID]
    C2_2 = -dpdRhoThL[M] * Surf[M,iz,ID]
    C3_1 = -cS;                   C3_2 = -cS

    # Column j = iz - 1
    j = i_m1
    r11 = B1_1
    r21 = B2_1
    @unroll for i = 1 : M
      a32i1 = calcA32(DWS,dpdRhoThL[1],i,1,Val(M))
      r3[i] = -(A31[i,1,iz,ID] * r11 + a32i1 * r21) * facLoc
    end
    r3[1] += B3_1
    ldivFull!(iz,ID,SA, r3, Val(M))
    r2M = eltype(SA)(0)
    @unroll for k = 1 : M
      a231k = DWS[1,k] * ThL[k] 
      a23Mk = DWS[M,k] * ThL[k] 
      r21 -= a231k * r3[k]
      r2M -= a23Mk * r3[k]
    end
    r21 *= facLoc
    r2M *= facLoc
    update_schur!(SchurBand, C2_1, C3_1, r21, r3[1], i_m1, j, ID)
    update_schur!(SchurBand, C2_2, C3_2, r2M, r3[M], i_p0, j, ID)

    # Column j = iz 
    j = i_p0
    r1M = B1_2
    r2M = B2_2
    @unroll for i = 1 : M
      a32iM = calcA32(DWS,dpdRhoThL[M],i,M,Val(M))
      r3[i] = -(A31[i,M,iz,ID] * r1M + a32iM * r2M) * facLoc
    end
    r3[M] += B3_2
    ldivFull!(iz,ID,SA, r3, Val(M))

    r21 = eltype(SA)(0)
    @unroll for k = 1 : M
      a231k = DWS[1,k] * ThL[k] 
      a23Mk = DWS[M,k] * ThL[k] 
      r21 -= a231k * r3[k]
      r2M -= a23Mk * r3[k]
    end
    r21 *= facLoc
    r2M *= facLoc
    update_schur!(SchurBand, C2_1, C3_1, r21, r3[1], i_m1, j, ID)
    update_schur!(SchurBand, C2_2, C3_2, r2M, r3[M], i_p0, j, ID)
  end
  if iz == nz
    Thm = U[M,iz - 1,ID,ThPos] / U[M,iz - 1,ID,RhoPos]
    B1_1 = -invwB;       B1_2 = zero(FTe)
    B2_1 = -FTe(0.5) * (ThL[1] + Thm) * invwB; B2_2 = zero(FTe)
    B3_1 = -invwB * cS;  B3_2 = zero(FTe)
    C2_1 = dpdRhoThL[1] * Surf[1,iz,ID]
    C2_2 = zero(FTe)
    C3_1 = -cS;                   C3_2 = zero(FTe)
    # Column j = iz - 1
    j = i_m1
    r11 = B1_1
    r21 = B2_1
    @unroll for i = 1 : M
      a32i1 = calcA32(DWS,dpdRhoThL[1],i,1,Val(M))
      r3[i] = -(A31[i,1,iz,ID] * r11 + a32i1 * r21) * facLoc
    end
    r3[1] += B3_1
    ldivFull!(iz,ID,SA, r3, Val(M))
    @unroll for k = 1 : M
      r21 -=  DWS[1,k] * ThL[k] * r3[k]
    end
    r21 *= facLoc
    update_schur!(SchurBand, C2_1, C3_1, r21, r3[1], i_m1, j, ID)
  end  
end


# Helper inline function to handle atomic updates cleanly
@inline function update_rs!(rs, C2_val, C3_val, r2_val, r3_val, i, ID)
  t = r2_val * C2_val + r3_val * C3_val
  @atomic rs[i, ID] -= t
end


@kernel inbounds = true function ldivHDGVerticalFKernel!(@Const(A31),@Const(Th),@Const(dpdRhoTh),
  @Const(SA),@Const(b),rs,@Const(dz),@Const(DW),@Const(w),fac,cS, ::Val{M}) where {M}

  iz, ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  ND = @uniform @ndrange()[2]

  @uniform FT = eltype(SA)

  r1 = @private FT (M,)
  r2 = @private FT (M,)
  r3 = @private FT (M,)
  DWS = @localmem FT (M,M)

  @uniform RhoPos = 1
  @uniform wPos = 4
  @uniform ThPos = 5
  @uniform wB = w[1]
  @uniform invwB = FT(1) / wB

  if iz == 1
    @. DWS = DW
  end
  @synchronize

  if ID <= ND
    dz2 = dz[iz,ID] * FT(0.5)
    facLoc = fac / dz2

    @unroll for i = 1 : M
      r1[i] = b[i,iz,ID,RhoPos] * dz2
      r2[i] = b[i,iz,ID,ThPos] * dz2
    end  

    for i = 1 : M
      r3[i] = b[i,iz,ID,wPos] * dz2
      for j = 1 : M
        a31ij = A31[i,j,iz,ID]
        a32ij = calcA32(DWS,dpdRhoTh[j,iz,ID],i,j,Val(M))
        r3[i] -= (a31ij * r1[j] + a32ij * r2[j]) * facLoc
      end
    end

    ldivFull!(iz, ID, SA, r3, Val(M))

    r21 = r2[1]
    r2M = r2[M]
    @unroll for j = 1 : M
      a231j = DWS[1,j] * Th[j,iz,ID] 
      a23Mj = DWS[M,j] * Th[j,iz,ID] 
      r21 -= a231j * r3[j]
      r2M -= a23Mj * r3[j]
    end
    r21 *= facLoc
    r2M *= facLoc

    # Pre-cache target row indices
    i_m1 = iz - 1
    i_p0 = iz

    if iz < nz
      C2_2 = -dpdRhoTh[M,iz,ID]
      C3_2 = -cS
      update_rs!(rs, C2_2, C3_2, r2M, r3[M], iz, ID)
    end  
    if iz > 1
      C2_1 = dpdRhoTh[1,iz,ID];
      C3_1 = -cS;
      update_rs!(rs, C2_1, C3_1, r21, r3[1], iz-1, ID)
    end
  end
end

@kernel inbounds = true function ldivHDGVerticalBKernel!(@Const(A31),@Const(Th),@Const(dpdRhoTh),
  @Const(SA),b,@Const(rs),@Const(dz),@Const(DW),@Const(w),fac,cS, ::Val{M}) where {M}

  iz,ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  ND = @uniform @ndrange()[2]

  @uniform FT = eltype(SA)

  r1 = @private FT (M,)
  r2 = @private FT (M,)
  r3 = @private FT (M,)
  DWS = @localmem FT (M,M)

  @uniform RhoPos = 1
  @uniform wPos = 4
  @uniform ThPos = 5
  @uniform wB = w[1]
  @uniform invwB = FT(1) / wB

  if iz == 1
    @. DWS = DW
  end
  @synchronize

  if ID <= ND
    dz2 = dz[iz,ID] * FT(0.5)
    facLoc = fac / dz2
    Thm = iz > 1  ? Th[M,iz-1,ID] : zero(FT)
    Thp = iz < nz ? Th[1,iz+1,ID] : zero(FT)

    @unroll for i = 1 : M
      r1[i] = b[i,iz,ID,RhoPos] * dz2
      r2[i] = b[i,iz,ID,ThPos] * dz2
      r3[i] = b[i,iz,ID,wPos] * dz2
    end  
    if iz <  nz
      B1_2 = invwB
      B2_2 = FT(0.5) * (Th[M,iz,ID] + Thp) * invwB
      B3_2 = -invwB * cS  
      j = iz
      r1[M] -= B1_2 * rs[j,ID]
      r2[M] -= B2_2 * rs[j,ID]
      r3[M] -= B3_2 * rs[j,ID]
    end
    if iz > 1
      B1_1 = -invwB
      B2_1 = -FT(0.5) * (Th[1,iz,ID] + Thm) * invwB
      B3_1 = -invwB * cS
      j = iz - 1
      r1[1] -= B1_1 * rs[j,ID]
      r2[1] -= B2_1 * rs[j,ID]
      r3[1] -= B3_1 * rs[j,ID]
    end    

    for i = 1 : M
      r3i = zero(eltype(SA))
      @unroll for j = 1 : M
        a31ij = A31[i,j,iz,ID]
        a32ij = calcA32(DWS,dpdRhoTh[j,iz,ID],i,j,Val(M))
        r3i += (a31ij * r1[j] + a32ij * r2[j])
      end
      r3[i] -= r3i * facLoc
    end

    ldivFull!(iz,ID,SA,r3,Val(M))
    @unroll for i = 1 : M
      @unroll for j = 1 : M
        a23ij = DWS[i,j] * Th[j,iz,ID]
        r1[i] = (r1[i] - DWS[i,j] * r3[j])
        r2[i] = (r2[i] - a23ij * r3[j])
      end  
      r1[i] *= facLoc
      r2[i] *= facLoc
    end
    @unroll for i = 1 : M
      b[i,iz,ID,RhoPos] = r1[i]
      b[i,iz,ID,ThPos] = r2[i]
      b[i,iz,ID,wPos] = r3[i]
    end
  end  
end  

function FillJacHDGJVert!(Jac::JacHDGJVert,U,DG,dz,fac,Phys)
  
  backend = get_backend(U)
  FTB = eltype(U)
  
  M = Jac.M
  nz = Jac.nz
  DoF  = DG.NumI 
  
  Jac.fac = fac
  Jac.cS = Phys.cS
  
  DWZ = DG.DWZ
  
  DoFG = 10
  group = (nz, DoFG)
  ndrange = (nz, DoF)
  @. Jac.SchurBand = 0
  @. Jac.SchurBand[2,:,:] = 2.0 * P.cS
  KFillJacHDGJVertKernel! = FillJacHDGJVertKernel!(backend,group)
  KFillJacHDGJVertKernel!(Jac.Th,Jac.dpdRhoTh,Jac.A31,Jac.SA,Jac.SchurBand,U,dz,
  DWZ,DG.wZ,fac,Phys.cS,Phys,Val(M);ndrange=ndrange)

end   

function SchurBoundary!(Jac::JacHDGJVert)

  backend = get_backend(Jac.SA)
  FTB = eltype(Jac.SA)

  M = Jac.M
  nz = Jac.nz
  ND = Jac.NumI

  NDG = 32
  group = (NDG)
  ndrange = (ND)
  KluBandKernel! = luBandKernel!(backend,group)
  KluBandKernel!(Jac.SchurBand,Val(1),Val(1),Val(nz-1),ndrange=ndrange)
end  


function Solve!(Jac::JacHDGJVert,b,DG,Metric)

  backend = get_backend(Jac.SA)
  FTB = eltype(Jac.SA)

  M = Jac.M
  nz = Jac.nz
  ND = Jac.NumI

  invfac = FTB(1) / Jac.fac
  fac = Jac.fac
  cS = Jac.cS
  dz = Metric.dz
  DWZ = DG.DWZ
  wZ = DG.wZ

  NDG = 32
  group = (nz, NDG)
  ndrange = (nz, ND)
  @. Jac.rs = 0.0
  KldivVerticalFKernel! = ldivHDGVerticalFKernel!(backend,group)
  KldivVerticalFKernel!(Jac.A31,Jac.Th,Jac.dpdRhoTh,
    Jac.SA,b,Jac.rs,dz,DWZ,wZ,fac,cS,Val(M);ndrange=ndrange)

  group = (NDG)
  ndrange = (ND)
  KldivVerticalSKernel! = ldivVerticalSKernel!(backend,group)
  KldivVerticalSKernel!(Jac.SchurBand,Jac.rs,Val(1),Val(1),Val(nz-1);ndrange=ndrange)

  group = (nz, NDG)
  ndrange = (nz, ND)
  KldivVerticalBKernel! = ldivHDGVerticalBKernel!(backend,group)
  KldivVerticalBKernel!(Jac.A31,Jac.Th,Jac.dpdRhoTh,
    Jac.SA,b,Jac.rs,dz,DWZ,wZ,fac,cS,Val(M);ndrange=ndrange)

end

function Jac!(U,fac,DG,Metric,Phys,Cache,JCache::JacHDGJVert,Global,VelForm)
  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU
  if JCache.grav_do
    @views Geo = Cache.Aux[:,:,:,2]
#   @views Geo = Cache[:,:,:,2]
    precompute_gravity!(Geo,Metric.dz,DG.DWZ,JCache,NumberThreadGPU)
    JCache.grav_do = false
  end  
  dz = Metric.dz
  FillJacHDGJVert!(JCache,U,DG,dz,fac,Phys)
  SchurBoundary!(JCache)
end

function Solve!(k,v,Jac::JacHDGJVert,fac,DG::FiniteElements.DGElement,Metric,Global,VelForm)

  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU
  @. k = v
  @views TendVCart2VSp!(k,DG,Metric,NumberThreadGPU,VelForm)
  Solve!(Jac,k,DG,Metric)
  @views @. k[:,:,:,2:3] *= fac
  @views TendVSp2VCart!(k,DG,Metric,NumberThreadGPU,VelForm)
end


function precompute_gravity!(GeoPot,dz, DWZ,Jac::JacHDGJVert, NumberThreadGPU)
  (; A31) = Jac
  backend = get_backend(dz)
  M = Jac.M
  nz = Jac.nz
  NumI = Jac.NumI
  NumIG = min(div(NumberThreadGPU,M*nz),NumI)
  group = (M,nz,NumIG)
  ndrange = (M,nz,NumI)
  Kprecompute_gravityKernel! = precompute_gravityKernel!(backend,group)
  Kprecompute_gravityKernel!(A31,GeoPot,DWZ,dz,Val(M),ndrange=ndrange)
end

@kernel inbounds = true function precompute_gravityKernel!(A31,@Const(Geo),@Const(DS),@Const(dz), ::Val{M}) where {M}

  i,iz,ID = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]
  if ID <= ND
    acc = eltype(dz)(0)
    Geoi = Geo[i, iz, ID]
    @unroll for k in 1:M
      dphi = Geo[k, iz, ID] - Geoi
      A31[i,k,iz,ID] = 0.5 * DS[i, k] * dphi
      acc += DS[i, k] * dphi
    end
    A31[i, i, iz, ID] += 0.5 * acc
  end
end
