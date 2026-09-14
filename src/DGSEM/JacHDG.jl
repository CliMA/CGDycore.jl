mutable struct JacHDGVert{FT<:AbstractFloat,
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
  A13::AT4
  A23::AT4
  A31::AT4
  A32::AT4
# B part  
  B1::AT3
  B2::AT3
  B3::AT3
# C part
  C2::AT4
  C3::AT4

  SA::AT4

  SchurBand::AT3
  rs::AT2
end  

function JacHDGVert(backend,FT,M,nz,DG) 
  CompTri = false
  grav_do = true
  fac = 0
  FacGrav = 0
  NumI = DG.NumI
  A13 = KernelAbstractions.zeros(backend,FT,M,M,nz,NumI)
  A23 = KernelAbstractions.zeros(backend,FT,M,M,nz,NumI)
  A31 = KernelAbstractions.zeros(backend,FT,M,M,nz,NumI)
  A32 = KernelAbstractions.zeros(backend,FT,M,M,nz,NumI)
  B1 = KernelAbstractions.zeros(backend,FT,2,nz,NumI)
  B2 = KernelAbstractions.zeros(backend,FT,2,nz,NumI)
  B3 = KernelAbstractions.zeros(backend,FT,2,nz,NumI)
  C2 = KernelAbstractions.zeros(backend,FT,2,2,nz,NumI)
  C3 = KernelAbstractions.zeros(backend,FT,2,2,nz,NumI)
  SA = KernelAbstractions.zeros(backend,FT,M,M,nz,NumI)
  SchurBand = KernelAbstractions.zeros(backend,FT,7,2*nz,NumI)
  rs = KernelAbstractions.zeros(backend,FT,2*nz,NumI)

  return JacHDGVert{FT,
                   typeof(rs),
                   typeof(B1),
                   typeof(A13)}(

    CompTri,
    grav_do,
    NumI,
    M,
    nz,
    fac,
    FacGrav,
    A13,
    A23,
    A31,
    A32,
    B1,
    B2,
    B3,
    C2,
    C3,
    SA,
    SchurBand,
    rs,
  )
end  

# Helper inline function to handle the Schur complement atomic update cleanly
@inline function update_schur!(SchurBand, C2_val, C3_val, r2_val, r3_val, i, j, ID)
  t = r2_val * C2_val + r3_val * C3_val
  iB = 4 + i - j
  @atomic SchurBand[iB, j, ID] -= t
end

@kernel inbounds = true function FillJacHDGVertKernel!(A13,A23,@Const(A31),A32,
  B1,B2,B3,C2,C3,SA,SchurBand,@Const(U),@Const(dz),
  @Const(DW),@Const(w),fac,cS,Phys, ::Val{M}) where {M}

  iz, ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  DoF = @uniform @ndrange()[2]
  Th = @private eltype(SA) (M,)
  dpdRhoTh = @private eltype(SA) (M,)
  r3 = @private eltype(SA) (M,)
  DWS = @localmem eltype(SA) (M,M)
  SAL = @localmem eltype(SA) (M,M)

  @uniform RhoPos = 1
  @uniform ThPos = 5
  @uniform invcS = eltype(SA)(1) / cS
  @uniform wB = w[1]
  @uniform invcSwB = eltype(SA)(1) / (cS * wB)
  @uniform invwB = eltype(SA)(1) / wB

  if iz == 1
    @. DWS = DW
  end
  @synchronize

  if ID <= DoF
    kappa = Phys.kappa
    kexp = kappa / (eltype(SA)(1) - kappa)
    kfac = eltype(SA)(1) / (eltype(SA)(1) - kappa) * Phys.Rd
    inv2dz = eltype(SA)(2) / dz[iz,ID]
    invdz = eltype(SA)(1) / dz[iz,ID]
    Rdp0 = Phys.Rd / Phys.p0
    invfac = eltype(SA)(1) / fac

    @unroll for i = 1 : M
      Th[i] = U[i,iz,ID,ThPos] / U[i,iz,ID,RhoPos]
      dpdRhoTh[i] = kfac * (Rdp0 * U[i,iz,ID,ThPos])^kexp
    end

    @unroll for i = 1 : M
      @unroll for j = 1 : M
        A13[i,j,iz,ID] = inv2dz * DWS[i,j]
        A23[i,j,iz,ID] = inv2dz * DWS[i,j] * Th[j]
        A32[i,j,iz,ID] = inv2dz * DWS[i,j] * dpdRhoTh[j]
      end
    end

    if iz == 1
      Thp = U[1,iz + 1,ID,ThPos] / U[1,iz + 1,ID,RhoPos]
      B1[1,iz,ID] = zero(eltype(SA))
      B1[2,iz,ID] = invdz
      B2[1,iz,ID] = zero(eltype(SA))
      B2[2,iz,ID] = eltype(SA)(0.5) * (Th[M] + Thp) * invdz
      B3[1,iz,ID] = -inv2dz
      B3[2,iz,ID] = invdz
#w      
      C2[1,1,iz,ID] = zero(eltype(SA))
      C2[1,2,iz,ID] = -dpdRhoTh[M] * invcSwB
      C3[1,1,iz,ID] = zero(eltype(SA))
      C3[1,2,iz,ID] = -invwB
#p      
      C2[2,1,iz,ID] = -dpdRhoTh[1] * invwB
      C2[2,2,iz,ID] = -dpdRhoTh[M] * invwB
      C3[2,1,iz,ID] = invwB * cS
      C3[2,2,iz,ID] = -invwB * cS
    elseif iz == nz
      Thm = U[M,iz - 1,ID,ThPos] / U[M,iz - 1,ID,RhoPos]
      B1[1,iz,ID] = -invdz
      B1[2,iz,ID] = zero(eltype(SA))
      B2[1,iz,ID] = -eltype(SA)(0.5) * (Th[1] + Thm) * invdz
      B2[2,iz,ID] = zero(eltype(SA))
      B3[1,iz,ID] = -invdz
      B3[2,iz,ID] = inv2dz
#w      
      C2[1,1,iz,ID] = dpdRhoTh[1] * invcSwB
      C2[1,2,iz,ID] = zero(eltype(SA))
      C3[1,1,iz,ID] = -invwB
      C3[1,2,iz,ID] = zero(eltype(SA))
#p      
      C2[2,1,iz,ID] = -dpdRhoTh[1] * invwB
      C2[2,2,iz,ID] = -dpdRhoTh[M] * invwB
      C3[2,1,iz,ID] = invwB * cS
      C3[2,2,iz,ID] = -invwB * cS
    else
      # p = 0.5*(pL+pR)-0.5*cS*(wR-wL)  
      # w = 0.5*(wL+wR)-0.5/cS*(pR-pL)  
      Thm = U[M,iz - 1,ID,ThPos] / U[M,iz - 1,ID,RhoPos]
      Thp = U[1,iz + 1,ID,ThPos] / U[1,iz + 1,ID,RhoPos]
      B1[1,iz,ID] = -invdz
      B1[2,iz,ID] = invdz
      B2[1,iz,ID] = -eltype(SA)(0.5) * (Th[1] + Thm) * invdz
      B2[2,iz,ID] = eltype(SA)(0.5) * (Th[M] + Thp) * invdz
      B3[1,iz,ID] = -invdz
      B3[2,iz,ID] = invdz
#w
      C2[1,1,iz,ID] = dpdRhoTh[1] * invcSwB
      C2[1,2,iz,ID] = -dpdRhoTh[M] * invcSwB
      C3[1,1,iz,ID] = -invwB
      C3[1,2,iz,ID] = -invwB
#p      
      C2[2,1,iz,ID] = -dpdRhoTh[1] * invwB
      C2[2,2,iz,ID] = -dpdRhoTh[M] * invwB
      C3[2,1,iz,ID] = invwB * cS
      C3[2,2,iz,ID] = -invwB * cS
    end

    @unroll for i = 1 : M
      @unroll for j = 1 : M
        SAL[i,j] = zero(eltype(SA))
        @unroll for k = 1 : M
          SAL[i,j] -= A32[i,k,iz,ID] * A23[k,j,iz,ID] + A31[i,k,iz,ID] * A13[k,j,iz,ID]
        end
        SAL[i,j] *= fac
      end
      SAL[i,i] += invfac
    end

    LUFull!(SAL, Val(M))
    @unroll for i = 1 : M
      @unroll for j = 1 : M
        SA[i,j,iz,ID] = SAL[i,j]
      end
    end
  end

  # Cache common column indices
  i_m2 = 2 * iz - 2
  i_m1 = 2 * iz - 1
  i_p0 = 2 * iz
  i_p1 = 2 * iz + 1

  # CASE 1: iz == 1
  if iz == 1
    # Column j = 2 * iz - 1
    j = i_m1
    r3[1] = B3[1,iz,ID]
    @unroll for i = 2 : M; r3[i] = zero(eltype(SA)); end
    ldivFull!(iz, ID, SA, r3, Val(M))

    r21 = zero(eltype(SA))
    r2M = zero(eltype(SA))
    @unroll for k = 1 : M
      r21 -= A23[1,k,iz,ID] * r3[k]
      r2M -= A23[M,k,iz,ID] * r3[k]
    end
    r21 *= fac; r2M *= fac

    update_schur!(SchurBand, C2[2,1,iz,ID], C3[2,1,iz,ID], r21, r3[1], i_m1, j, ID)
    update_schur!(SchurBand, C2[1,2,iz,ID], C3[1,2,iz,ID], r2M, r3[M], i_p0, j, ID)
    update_schur!(SchurBand, C2[2,2,iz,ID], C3[2,2,iz,ID], r2M, r3[M], i_p1, j, ID)

    # Column j = 2 * iz
    j = i_p0
    r1M = B1[2,iz,ID]
    r2M_init = B2[2,iz,ID]
    @unroll for i = 1 : M
      r3[i] = -(A31[i,M,iz,ID] * r1M + A32[i,M,iz,ID] * r2M_init) * fac
    end
    ldivFull!(iz, ID, SA, r3, Val(M))

    r21 = zero(eltype(SA))
    r2M = r2M_init
    @unroll for k = 1 : M
      r21 -= A23[1,k,iz,ID] * r3[k]
      r2M -= A23[M,k,iz,ID] * r3[k]
    end
    r21 *= fac; r2M *= fac

    update_schur!(SchurBand, C2[2,1,iz,ID], C3[2,1,iz,ID], r21, r3[1], i_m1, j, ID)
    update_schur!(SchurBand, C2[1,2,iz,ID], C3[1,2,iz,ID], r2M, r3[M], i_p0, j, ID)
    update_schur!(SchurBand, C2[2,2,iz,ID], C3[2,2,iz,ID], r2M, r3[M], i_p1, j, ID)

    # Column j = 2 * iz + 1
    j = i_p1
    r3[M] = B3[2,iz,ID]
    @unroll for i = 1 : M - 1; r3[i] = zero(eltype(SA)); end
    ldivFull!(iz, ID, SA, r3, Val(M))

    r21 = zero(eltype(SA))
    r2M = zero(eltype(SA))
    @unroll for k = 1 : M
      r21 -= A23[1,k,iz,ID] * r3[k]
      r2M -= A23[M,k,iz,ID] * r3[k]
    end
    r21 *= fac; r2M *= fac

    update_schur!(SchurBand, C2[2,1,iz,ID], C3[2,1,iz,ID], r21, r3[1], i_m1, j, ID)
    update_schur!(SchurBand, C2[1,2,iz,ID], C3[1,2,iz,ID], r2M, r3[M], i_p0, j, ID)
    update_schur!(SchurBand, C2[2,2,iz,ID], C3[2,2,iz,ID], r2M, r3[M], i_p1, j, ID)
  end

  # CASE 2: iz > 1 && iz < nz
  if iz > 1 && iz < nz
    # Column j = 2 * iz - 2
    j = i_m2
    r11 = B1[1,iz,ID]
    r21_init = B2[1,iz,ID]
    @unroll for i = 1 : M
      r3[i] = -(A31[i,1,iz,ID] * r11 + A32[i,1,iz,ID] * r21_init) * fac
    end
    ldivFull!(iz, ID, SA, r3, Val(M))

    r21 = r21_init
    r2M = zero(eltype(SA))
    @unroll for k = 1 : M
      r21 -= A23[1,k,iz,ID] * r3[k]
      r2M -= A23[M,k,iz,ID] * r3[k]
    end
    r21 *= fac; r2M *= fac

    update_schur!(SchurBand, C2[1,1,iz,ID], C3[1,1,iz,ID], r21, r3[1], i_m2, j, ID)
    update_schur!(SchurBand, C2[2,1,iz,ID], C3[2,1,iz,ID], r21, r3[1], i_m1, j, ID)
    update_schur!(SchurBand, C2[1,2,iz,ID], C3[1,2,iz,ID], r2M, r3[M], i_p0, j, ID)
    update_schur!(SchurBand, C2[2,2,iz,ID], C3[2,2,iz,ID], r2M, r3[M], i_p1, j, ID)

    # Column j = 2 * iz - 1
    j = i_m1
    r3[1] = B3[1,iz,ID]
    @unroll for i = 2 : M; r3[i] = zero(eltype(SA)); end
    ldivFull!(iz, ID, SA, r3, Val(M))

    r21 = zero(eltype(SA))
    r2M = zero(eltype(SA))
    @unroll for k = 1 : M
      r21 -= A23[1,k,iz,ID] * r3[k]
      r2M -= A23[M,k,iz,ID] * r3[k]
    end
    r21 *= fac; r2M *= fac

    update_schur!(SchurBand, C2[1,1,iz,ID], C3[1,1,iz,ID], r21, r3[1], i_m2, j, ID)
    update_schur!(SchurBand, C2[2,1,iz,ID], C3[2,1,iz,ID], r21, r3[1], i_m1, j, ID)
    update_schur!(SchurBand, C2[1,2,iz,ID], C3[1,2,iz,ID], r2M, r3[M], i_p0, j, ID)
    update_schur!(SchurBand, C2[2,2,iz,ID], C3[2,2,iz,ID], r2M, r3[M], i_p1, j, ID)

    # Column j = 2 * iz
    j = i_p0
    r1M = B1[2,iz,ID]
    r2M_init = B2[2,iz,ID]
    @unroll for i = 1 : M
      r3[i] = -(A31[i,M,iz,ID] * r1M + A32[i,M,iz,ID] * r2M_init) * fac
    end
    ldivFull!(iz, ID, SA, r3, Val(M))

    r21 = zero(eltype(SA))
    r2M = r2M_init
    @unroll for k = 1 : M
      r21 -= A23[1,k,iz,ID] * r3[k]
      r2M -= A23[M,k,iz,ID] * r3[k]
    end
    r21 *= fac; r2M *= fac

    update_schur!(SchurBand, C2[1,1,iz,ID], C3[1,1,iz,ID], r21, r3[1], i_m2, j, ID)
    update_schur!(SchurBand, C2[2,1,iz,ID], C3[2,1,iz,ID], r21, r3[1], i_m1, j, ID)
    update_schur!(SchurBand, C2[1,2,iz,ID], C3[1,2,iz,ID], r2M, r3[M], i_p0, j, ID)
    update_schur!(SchurBand, C2[2,2,iz,ID], C3[2,2,iz,ID], r2M, r3[M], i_p1, j, ID)

    # Column j = 2 * iz + 1
    j = i_p1
    r3[M] = B3[2,iz,ID]
    @unroll for i = 1 : M - 1; r3[i] = zero(eltype(SA)); end
    ldivFull!(iz, ID, SA, r3, Val(M))

    r21 = zero(eltype(SA))
    r2M = zero(eltype(SA))
    @unroll for k = 1 : M
      r21 -= A23[1,k,iz,ID] * r3[k]
      r2M -= A23[M,k,iz,ID] * r3[k]
    end
    r21 *= fac; r2M *= fac

    update_schur!(SchurBand, C2[1,1,iz,ID], C3[1,1,iz,ID], r21, r3[1], i_m2, j, ID)
    update_schur!(SchurBand, C2[2,1,iz,ID], C3[2,1,iz,ID], r21, r3[1], i_m1, j, ID)
    update_schur!(SchurBand, C2[1,2,iz,ID], C3[1,2,iz,ID], r2M, r3[M], i_p0, j, ID)
    update_schur!(SchurBand, C2[2,2,iz,ID], C3[2,2,iz,ID], r2M, r3[M], i_p1, j, ID)
  end

  # CASE 3: iz == nz
  if iz == nz
    # Column j = 2 * iz - 2
    j = i_m2
    r11 = B1[1,iz,ID]
    r21_init = B2[1,iz,ID]
    @unroll for i = 1 : M
      r3[i] = -(A31[i,1,iz,ID] * r11 + A32[i,1,iz,ID] * r21_init) * fac
    end
    ldivFull!(iz, ID, SA, r3, Val(M))

    r21 = r21_init
    r2M = zero(eltype(SA))
    @unroll for k = 1 : M
      r21 -= A23[1,k,iz,ID] * r3[k]
      r2M -= A23[M,k,iz,ID] * r3[k]
    end
    r21 *= fac; r2M *= fac

    update_schur!(SchurBand, C2[1,1,iz,ID], C3[1,1,iz,ID], r21, r3[1], i_m2, j, ID)
    update_schur!(SchurBand, C2[2,1,iz,ID], C3[2,1,iz,ID], r21, r3[1], i_m1, j, ID)
    update_schur!(SchurBand, C2[2,2,iz,ID], C3[2,2,iz,ID], r2M, r3[M], i_p0, j, ID)

    # Column j = 2 * iz - 1
    j = i_m1
    r3[1] = B3[1,iz,ID]
    @unroll for i = 2 : M; r3[i] = zero(eltype(SA)); end
    ldivFull!(iz, ID, SA, r3, Val(M))

    r21 = zero(eltype(SA))
    r2M = zero(eltype(SA))
    @unroll for k = 1 : M
      r21 -= A23[1,k,iz,ID] * r3[k]
      r2M -= A23[M,k,iz,ID] * r3[k]
    end
    r21 *= fac; r2M *= fac

    update_schur!(SchurBand, C2[1,1,iz,ID], C3[1,1,iz,ID], r21, r3[1], i_m2, j, ID)
    update_schur!(SchurBand, C2[2,1,iz,ID], C3[2,1,iz,ID], r21, r3[1], i_m1, j, ID)
    update_schur!(SchurBand, C2[2,2,iz,ID], C3[2,2,iz,ID], r2M, r3[M], i_p0, j, ID)

    # Column j = 2 * iz
    j = i_p0
    r3[M] = B3[2,iz,ID]
    @unroll for i = 1 : M - 1; r3[i] = zero(eltype(SA)); end
    ldivFull!(iz, ID, SA, r3, Val(M))

    r21 = zero(eltype(SA))
    r2M = zero(eltype(SA))
    @unroll for k = 1 : M
      r21 -= A23[1,k,iz,ID] * r3[k]
      r2M -= A23[M,k,iz,ID] * r3[k]
    end
    r21 *= fac; r2M *= fac

    update_schur!(SchurBand, C2[1,1,iz,ID], C3[1,1,iz,ID], r21, r3[1], i_m2, j, ID)
    update_schur!(SchurBand, C2[2,1,iz,ID], C3[2,1,iz,ID], r21, r3[1], i_m1, j, ID)
    update_schur!(SchurBand, C2[2,2,iz,ID], C3[2,2,iz,ID], r2M, r3[M], i_p0, j, ID)
  end
end


# Helper inline function to handle atomic updates cleanly
@inline function update_rs!(rs, C2_val, C3_val, r2_val, r3_val, i, ID)
  t = r2_val * C2_val + r3_val * C3_val
  @atomic rs[i, ID] -= t
end


@kernel inbounds = true function ldivHDGVerticalFKernel!(@Const(A13),@Const(A23),@Const(A31),@Const(A32),
  @Const(C2), @Const(C3),@Const(SA),@Const(b),rs,fac, ::Val{M}) where {M}

  iz, ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  ND = @uniform @ndrange()[2]

  @uniform FT = eltype(SA)

  r1 = @private FT (M,)
  r2 = @private FT (M,)
  r3 = @private FT (M,)

  @uniform RhoPos = 1
  @uniform wPos = 4
  @uniform ThPos = 5

  if ID <= ND

    @unroll for i = 1 : M
      r1[i] = b[i,iz,ID,RhoPos]
      r2[i] = b[i,iz,ID,ThPos]
#     r3[i] = b[i,iz,ID,wPos]
    end  

    @unroll for i = 1 : M
      r3i = zero(eltype(SA))
      @unroll for j = 1 : M
        r3i += (A31[i,j,iz,ID] * r1[j] + A32[i,j,iz,ID] * r2[j])
      end
#     r3[i] -= r3i * fac
      r3[i] = b[i,iz,ID,wPos] - r3i * fac
    end
#=
    for i = 1 : M
      r3i = zero(eltype(SA))
      for j = 1 : M
        r3i += (A31[i,j,iz,ID] * r1[j] + A32[i,j,iz,ID] * r2[j])
      end
      r3[i] = b[i,iz,ID,wPos] - r3i * fac
    end
    for i = 1 : M
      r3[i] = b[i,iz,ID,wPos]
      for j = 1 : M
        r3[i] -= (A31[i,j,iz,ID] * r1[j] + A32[i,j,iz,ID] * r2[j]) * fac
      end
    end
=#    

    ldivFull!(iz, ID, SA, r3, Val(M))

    r21 = r2[1]
    r2M = r2[M]
    @unroll for j = 1 : M
      r21 -= A23[1,j,iz,ID] * r3[j]
      r2M -= A23[M,j,iz,ID] * r3[j]
    end
    r21 *= fac
    r2M *= fac

    # Pre-cache target row indices
    i_m2 = 2 * iz - 2
    i_m1 = 2 * iz - 1
    i_p0 = 2 * iz
    i_p1 = 2 * iz + 1

    if iz == 1
      update_rs!(rs, C2[2,1,iz,ID], C3[2,1,iz,ID], r21, r3[1], i_m1, ID)
      update_rs!(rs, C2[1,2,iz,ID], C3[1,2,iz,ID], r2M, r3[M], i_p0, ID)
      update_rs!(rs, C2[2,2,iz,ID], C3[2,2,iz,ID], r2M, r3[M], i_p1, ID)
    elseif iz > 1 && iz < nz
      update_rs!(rs, C2[1,1,iz,ID], C3[1,1,iz,ID], r21, r3[1], i_m2, ID)
      update_rs!(rs, C2[2,1,iz,ID], C3[2,1,iz,ID], r21, r3[1], i_m1, ID)
      update_rs!(rs, C2[1,2,iz,ID], C3[1,2,iz,ID], r2M, r3[M], i_p0, ID)
      update_rs!(rs, C2[2,2,iz,ID], C3[2,2,iz,ID], r2M, r3[M], i_p1, ID)
    elseif iz == nz
      update_rs!(rs, C2[1,1,iz,ID], C3[1,1,iz,ID], r21, r3[1], i_m2, ID)
      update_rs!(rs, C2[2,1,iz,ID], C3[2,1,iz,ID], r21, r3[1], i_m1, ID)
      update_rs!(rs, C2[2,2,iz,ID], C3[2,2,iz,ID], r2M, r3[M], i_p0, ID)
    end
  end
end

@kernel inbounds = true function ldivHDGVerticalBKernel!(@Const(A13),@Const(A23),@Const(A31),@Const(A32),
  @Const(B1),@Const(B2),@Const(B3),@Const(SA),b,@Const(rs),fac,  ::Val{M}) where {M}

  iz,ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  ND = @uniform @ndrange()[2]

  @uniform FT = eltype(SA)

  r1 = @private FT (M,)
  r2 = @private FT (M,)
  r3 = @private FT (M,)

  @uniform RhoPos = 1
  @uniform wPos = 4
  @uniform ThPos = 5

  if ID <= ND

    @unroll for i = 1 : M
      r1[i] = b[i,iz,ID,RhoPos]
      r2[i] = b[i,iz,ID,ThPos]
      r3[i] = b[i,iz,ID,wPos]
    end  
    if iz == 1
      j = 2 * iz - 1
      r3[1] -= B3[1,iz,ID] * rs[j,ID]

      j = 2 * iz 
      r1[M] -= B1[2,iz,ID] * rs[j,ID]
      r2[M] -= B2[2,iz,ID] * rs[j,ID]

      j = 2 * iz + 1
      r3[M] -= B3[2,iz,ID] * rs[j,ID]

    end
    if iz > 1 && iz < nz
      j = 2 * iz - 2
      r1[1] -= B1[1,iz,ID] * rs[j,ID]
      r2[1] -= B2[1,iz,ID] * rs[j,ID]


      j = 2 * iz - 1
      r3[1] -= B3[1,iz,ID] * rs[j,ID]

      j = 2 * iz
      r1[M] -= B1[2,iz,ID] * rs[j,ID]
      r2[M] -= B2[2,iz,ID] * rs[j,ID]

      j = 2 * iz + 1
      r3[M] -= B3[2,iz,ID] * rs[j,ID]
    end  
    if iz == nz
      j = 2 * iz - 2
      r1[1] -= B1[1,iz,ID] * rs[j,ID]
      r2[1] -= B2[1,iz,ID] * rs[j,ID]

      j = 2 * iz - 1
      r3[1] -= B3[1,iz,ID] * rs[j,ID]

      j = 2 * iz
      r3[M] -= B3[2,iz,ID] * rs[j,ID]
    end    

    for i = 1 : M
      r3i = zero(eltype(SA))
      @unroll for j = 1 : M
        r3i += (A31[i,j,iz,ID] * r1[j] + A32[i,j,iz,ID] * r2[j])
      end
      r3[i] -= r3i * fac
    end

    ldivFull!(iz,ID,SA,r3,Val(M))
    @unroll for i = 1 : M
      @unroll for j = 1 : M
        r1[i] = (r1[i] - A13[i,j,iz,ID] * r3[j])
        r2[i] = (r2[i] - A23[i,j,iz,ID] * r3[j])
      end  
      r1[i] *= fac
      r2[i] *= fac
    end
    @unroll for i = 1 : M
      b[i,iz,ID,RhoPos] = r1[i]
      b[i,iz,ID,ThPos] = r2[i]
      b[i,iz,ID,wPos] = r3[i]
    end
  end  
end  

function FillJacHDGVert!(Jac::JacHDGVert,U,DG,dz,fac,Phys)
  
  backend = get_backend(U)
  FTB = eltype(U)
  
  M = Jac.M
  nz = Jac.nz
  DoF  = DG.NumI 
  
  Jac.fac = fac
  
  DWZ = DG.DWZ
  
  DoFG = 10
  group = (nz, DoFG)
  ndrange = (nz, DoF)
  @. Jac.SchurBand = 0
  @. Jac.SchurBand[4,:,:] = 1.0 # Diagonal for the extended part
  KFillJacHDGVertKernel! = FillJacHDGVertKernel!(backend,group)
  KFillJacHDGVertKernel!(Jac.A13,Jac.A23,Jac.A31,Jac.A32,Jac.B1,Jac.B2,Jac.B3,
  Jac.C2,Jac.C3,Jac.SA,Jac.SchurBand,U,dz,
  DWZ,DG.wZ,fac,Phys.cS,Phys,Val(M);ndrange=ndrange)

end   

function SchurBoundary!(Jac::JacHDGVert)

  backend = get_backend(Jac.SA)
  FTB = eltype(Jac.SA)

  M = Jac.M
  nz = Jac.nz
  ND = Jac.NumI

  NDG = 32
  group = (NDG)
  ndrange = (ND)
  KluBandKernel! = luBandKernel!(backend,group)
  KluBandKernel!(Jac.SchurBand,Val(2*nz),ndrange=ndrange)
end  


function Solve!(Jac::JacHDGVert,b)

  backend = get_backend(Jac.SA)
  FTB = eltype(Jac.SA)

  M = Jac.M
  nz = Jac.nz
  ND = Jac.NumI

  invfac = FTB(1) / Jac.fac
  fac = Jac.fac

  NDG = 32
  group = (nz, NDG)
  ndrange = (nz, ND)
  @. Jac.rs = 0.0
  KldivVerticalFKernel! = ldivHDGVerticalFKernel!(backend,group)
  KldivVerticalFKernel!(Jac.A13,Jac.A23,Jac.A31,Jac.A32,Jac.C2,Jac.C3,
    Jac.SA,b,Jac.rs,fac,Val(M);ndrange=ndrange)

  group = (NDG)
  ndrange = (ND)
  KldivVerticalSKernel! = ldivVerticalSKernel!(backend,group)
  KldivVerticalSKernel!(Jac.SchurBand,Jac.rs,Val(2*nz);ndrange=ndrange)

  group = (nz, NDG)
  ndrange = (nz, ND)
  KldivVerticalBKernel! = ldivHDGVerticalBKernel!(backend,group)
  KldivVerticalBKernel!(Jac.A13,Jac.A23,Jac.A31,Jac.A32,Jac.B1,Jac.B2,Jac.B3,
    Jac.SA,b,Jac.rs,fac,Val(M);ndrange=ndrange)

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
  Solve!(Jac,k)
  @views @. k[:,:,:,2:3] *= fac
  @views TendVSp2VCart!(k,DG,Metric,NumberThreadGPU,VelForm)
end


function precompute_gravity!(GeoPot,dz, DWZ,Jac::JacHDGVert, NumberThreadGPU)
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
    inv2dz = eltype(A31)(2) / dz[iz, ID]
    @unroll for k in 1:M
      dphi = Geo[k, iz, ID] - Geoi
      A31[i,k,iz,ID] = inv2dz * 0.5 * DS[i, k] * dphi
      acc += DS[i, k] * dphi
    end
    A31[i, i, iz, ID] += inv2dz * 0.5 * acc
  end
end

