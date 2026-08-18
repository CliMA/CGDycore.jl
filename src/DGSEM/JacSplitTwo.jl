@inline function LUFull!(A,::Val{n}) where {n}

  @inbounds for j = 1 : n - 1
    @inbounds for i = j + 1 : n
      A[i,j] /= A[j,j]
      @inbounds for k = j + 1 : n
        A[i,k] -= A[i,j] * A[j,k]
      end
    end
  end
end
mutable struct JacSplitDGVert{FT<:AbstractFloat,
                         AT2<:AbstractArray,
                         AT3<:AbstractArray,
                         AT4<:AbstractArray}
  CompTri::Bool                        
  grav_do::Bool
  M::Int
  nz::Int
  fac::FT
  FacGrav::FT
  A13::AT4
  A23::AT4
  A31::AT4
  A32::AT4
# B part  
  B1m_34::AT3
  B1_1::AT2
  B1_23::AT4
  B1_4::AT2
  B1p_12::AT3
  B2_23::AT4
  B3_14::AT4
# C part
  C23_1::AT4
  C23_2::AT4
  C14_3::AT4

  SA::AT4

  DL::AT4
  DD::AT4
  DU::AT4
  rs::AT3
end  

function JacSplitDGVert(backend,FT,M,nz,DG) 
  CompTri = false
  grav_do = true
  M2 = M - 2
  fac = 0
  FacGrav = 0
  NumI = DG.NumI
  A13 = KernelAbstractions.zeros(backend,FT,M,M2,nz,NumI)
  A23 = KernelAbstractions.zeros(backend,FT,M2,M2,nz,NumI)
  A31 = KernelAbstractions.zeros(backend,FT,M2,M,nz,NumI)
  A32 = KernelAbstractions.zeros(backend,FT,M2,M2,nz,NumI)
  B1m_34 = KernelAbstractions.zeros(backend,FT,2,nz,NumI)
  B1_1 = KernelAbstractions.zeros(backend,FT,nz,NumI)
  B1_23 = KernelAbstractions.zeros(backend,FT,M,2,nz,NumI)
  B1_4 = KernelAbstractions.zeros(backend,FT,nz,NumI)
  B2_23 = KernelAbstractions.zeros(backend,FT,M2,2,nz,NumI)
  B3_14 = KernelAbstractions.zeros(backend,FT,M2,2,nz,NumI)
  B1p_12 = KernelAbstractions.zeros(backend,FT,2,nz,NumI)
  C23_1 = KernelAbstractions.zeros(backend,FT,2,M,nz,NumI)
  C23_2 = KernelAbstractions.zeros(backend,FT,2,M2,nz,NumI)
  C14_3 = KernelAbstractions.zeros(backend,FT,2,M2,nz,NumI)
  SA = KernelAbstractions.zeros(backend,FT,M2,M2,nz,NumI)
  DL = KernelAbstractions.zeros(backend,FT,NumI,4,2,nz-1)
  DD = KernelAbstractions.zeros(backend,FT,NumI,4,4,nz)
  DU = KernelAbstractions.zeros(backend,FT,NumI,4,2,nz-1)
  rs = KernelAbstractions.zeros(backend,FT,NumI,4,nz)

  return JacSplitDGVert{FT,
                   typeof(B1_1),
                   typeof(B1m_34),
                   typeof(A13)}(

    CompTri,
    grav_do,
    M,
    nz,
    fac,
    FacGrav,
    A13,
    A23,
    A31,
    A32,
    B1m_34,
    B1_1,
    B1_23,
    B1_4,
    B1p_12,
    B2_23,
    B3_14,
    C23_1,
    C23_2,
    C14_3,
    SA,
    DL,
    DD,
    DU,
    rs,
  )
end  

@kernel inbounds = true function FillJacSplitDGVertKernel!(A13,A23,@Const(A31),A32,B1m_34,B1_1,
  B1_23,B1_4, B2_23,B3_14,B1p_12,@Const(C23_1),C23_2,C14_3,SA,DL,DD,DU,@Const(U),@Const(dz),
  @Const(DW),@Const(w),fac,cS,Phys, ::Val{M}) where {M}

  iz,ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  DoF = @uniform @ndrange()[2]
  Th = @private eltype(U) (M,)
  dpdRhoTh = @private eltype(U) (M,)
  DWS = @localmem eltype(U) (M,M)
  SAL = @localmem eltype(U) (M-2,M-2)
  @uniform RhoPos = 1
  @uniform ThPos = 5
  @uniform invcS = eltype(U)(1) / cS
  @uniform wB = w[1]

  if iz == 1
    @. DWS = DW   
  end
  @synchronize 


  if ID <= DoF
    kappa = Phys.kappa
    kexp = kappa / (eltype(U)(1) - kappa)
    kfac = eltype(U)(1) / (eltype(U)(1) - kappa) * Phys.Rd
    sh = (iz - 1) * 4 + 1
    inv2dz = eltype(U)(2) / dz[iz,ID]
    invdz = eltype(U)(1) / dz[iz,ID]
    invwBdz = invdz / wB
    Rdp0 = Phys.Rd / Phys.p0
    invfac = eltype(U)(1) / fac
    @unroll for i = 1 : M
      Th[i] = U[i,iz,ID,ThPos] / U[i,iz,ID,RhoPos]
      dpdRhoTh[i] = kfac * (Rdp0 * U[i,iz,ID,ThPos])^kexp
    end      
    @unroll for i = 1 : M 
      @unroll for j = 1 : M - 2
        A13[i,j,iz,ID] = inv2dz * DWS[i,j+1]
      end
    end  
    @unroll for i = 1 : M - 2
      @unroll for j = 1 : M - 2
        A23[i,j,iz,ID] = inv2dz * DWS[i+1,j+1] * Th[j+1]  
        A32[i,j,iz,ID] = inv2dz * DWS[i+1,j+1] * dpdRhoTh[j+1]  
      end
    end  
#   @unroll for i = 1 : M - 2 
#     @unroll for j = 1 : M 
#       A31[i,j,iz,ID] = Gblk[i+1,j,iz,ID]
#     end
#   end  

    @unroll for i = 1 : M - 2
      @unroll for j = 1 : M - 2
        SAL[i,j] = eltype(U)(0)
        @unroll for k = 1 : M - 2
          SAL[i,j] += A32[i,k,iz,ID] * A23[k,j,iz,ID]
        end
        @unroll for k = 1 : M
          SAL[i,j] += A31[i,k,iz,ID] * A13[k,j,iz,ID]
        end
        if i == j
          SAL[i,j] = fac - SAL[i,j] * invfac
        else
          SAL[i,j] *= -invfac 
        end
      end
    end
    LUFull!(SAL,Val(M - 2))
    @unroll for i = 1 : M - 2
      @unroll for j = 1 : M - 2
        SA[i,j,iz,ID] = SAL[i,j]
      end
    end  
    if iz > 1
      B1m_34[1,iz,ID] = -invwBdz

      dpdRhoThM   = kfac * (Phys.Rd * U[M,iz-1,ID,ThPos] / Phys.p0)^kexp
      B1m_34[2,iz,ID] = -dpdRhoThM * invcS * invwBdz
    end
    @unroll for i = 1 : M
      B1_23[i,1,iz,ID] = eltype(U)(2) * DW[i,1] * invdz
      B1_23[i,2,iz,ID] = eltype(U)(2) * DW[i,M] * invdz
    end  

    if iz == 1
      B1_1[iz,ID] = eltype(U)(0)
    else
      B1_23[1,1,iz,ID] = eltype(U)(0)
      B1_1[iz,ID] = dpdRhoTh[1] * invcS * invwBdz
    end
    if iz == nz
      B1_4[iz,ID] = eltype(U)(0)
    else
      B1_23[M,2,iz,ID] = eltype(U)(0)
      B1_4[iz,ID] = dpdRhoTh[M] * invcS * invwBdz 
    end
    if iz < nz
      dpdRhoThP   = kfac * (Rdp0 * U[1,iz+1,ID,ThPos])^kexp
      B1p_12[1,iz,ID] = -dpdRhoThP * invcS * invwBdz
      B1p_12[2,iz,ID] = invwBdz
    end

#   @unroll for i = 1 : M
#     C23_1[1,i,iz,ID] = Gblk[1,i,iz,ID]
#     C23_1[2,i,iz,ID] = Gblk[M,i,iz,ID]
#   end  
    @unroll for i = 1 : M - 2
      B2_23[i,1,iz,ID] = inv2dz * DWS[i+1,1] * Th[1]
      B2_23[i,2,iz,ID] = inv2dz * DWS[i+1,M] * Th[M]
      B3_14[i,1,iz,ID] = inv2dz * DWS[i+1,1] * dpdRhoTh[1]
      B3_14[i,2,iz,ID] = inv2dz * DWS[i+1,M] * dpdRhoTh[M]
      C23_2[1,i,iz,ID] = inv2dz * DWS[1,i+1] * dpdRhoTh[i+1]
      C23_2[2,i,iz,ID] = inv2dz * DWS[M,i+1] * dpdRhoTh[i+1]
      C14_3[1,i,iz,ID] = inv2dz * DWS[1,i+1] * Th[i+1]
      C14_3[2,i,iz,ID] = inv2dz * DWS[M,i+1] * Th[i+1]
    end


    DD[ID,1,2,iz] = eltype(U)(0)
    DD[ID,1,4,iz] = eltype(U)(0)
    DD[ID,2,1,iz] = eltype(U)(0)
    DD[ID,2,3,iz] = eltype(U)(0)
    DD[ID,3,2,iz] = eltype(U)(0)
    DD[ID,3,4,iz] = eltype(U)(0)
    DD[ID,4,1,iz] = eltype(U)(0)
    DD[ID,4,3,iz] = eltype(U)(0)

    if iz > 1
      DU[ID,1,1,iz-1] = eltype(U)(0)
      DU[ID,1,2,iz-1] = eltype(U)(0)
      DU[ID,2,1,iz-1] = eltype(U)(0)
      DU[ID,2,2,iz-1] = eltype(U)(0)
      DL[ID,3,1,iz-1] = eltype(U)(0)
      DL[ID,3,2,iz-1] = eltype(U)(0)
      DL[ID,4,1,iz-1] = eltype(U)(0)
      DL[ID,4,2,iz-1] = eltype(U)(0)
      ThM = U[M,iz-1,ID,ThPos] / U[M,iz-1,ID,RhoPos]
      dpdRhoThM = eltype(U)(1) / (eltype(U)(1) - Phys.kappa) * Phys.Rd *
        (Phys.Rd * U[M,iz-1,ID,ThPos] / Phys.p0)^(Phys.kappa / (eltype(U)(1) - Phys.kappa))
      invwBdzm = eltype(U)(1) / (wB * dz[iz-1,ID])

      DD[ID,1,1,iz] = fac + Th[1] * dpdRhoTh[1] * invcS * invwBdz
      DD[ID,2,2,iz] = fac + cS * invwBdz

      DD[ID,3,3,iz-1] = fac + cS * invwBdzm
      DD[ID,4,4,iz-1] = fac + ThM * dpdRhoThM * invcS * invwBdzm

      DL[ID,1,1,iz-1] = -ThM * invwBdz
      DL[ID,1,2,iz-1] = -Th[1] * dpdRhoThM * invcS * invwBdz
      DL[ID,2,1,iz-1] = -cS * invwBdz
      DL[ID,2,2,iz-1] = -dpdRhoThM * invwBdz

      DU[ID,3,1,iz-1] = dpdRhoTh[1] * invwBdzm
      DU[ID,4,1,iz-1] = -ThM * dpdRhoTh[1] * invcS * invwBdzm
      DU[ID,3,2,iz-1] = -cS * invwBdzm
      DU[ID,4,2,iz-1] = Th[1] * invwBdzm
    end

    if iz == 1
      DD[ID,1,1,1] = fac 
      DD[ID,1,2,1] = inv2dz * DWS[1,1] * Th[1]
      DD[ID,2,1,1] = -inv2dz * DWS[1,1] * dpdRhoTh[1]
      DD[ID,2,2,1] = fac + inv2dz * cS / wB
    end
    DD[ID,1,3,iz] = inv2dz * DWS[1,M] * Th[M]
    DD[ID,3,1,iz] = inv2dz * DWS[M,1] * dpdRhoTh[1]
    DD[ID,4,2,iz] = inv2dz * DWS[M,1] * Th[1]
    DD[ID,2,4,iz] = inv2dz * DWS[1,M] * dpdRhoTh[M]
    if iz == nz
      DD[ID,3,3,nz] = fac + inv2dz * cS / wB
      DD[ID,3,4,nz] = -inv2dz * DWS[M,M] * dpdRhoTh[M]
      DD[ID,4,3,nz] = inv2dz * DWS[M,M] * Th[M]
      DD[ID,4,4,nz] = fac
    end
  end
end  


@kernel inbounds = true function SchurBoundaryKernel1!(@Const(A13),@Const(A23),@Const(A31),@Const(A32),
  @Const(B1m_34),@Const(B1_1),@Const(B1_23),@Const(B1_4), @Const(B2_23),@Const(B3_14),
  @Const(B1p_12),@Const(C23_1),@Const(C23_2),@Const(C14_3),@Const(SA),DU,DD,DL,fac, ::Val{M}) where {M}

  iz,ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  ND = @uniform @ndrange()[2]

  @uniform FT = eltype(DD)

  r1 = @private FT (M,2,)
  r2 = @private FT (M-2,2,)
  r3 = @private FT (M-2,2,)
  r11 = @private FT (2,)
  r1M = @private FT (2,)
  s = @private FT (4,2,)
  invfac = FT(1) / fac

  if ID <= ND
      
#   Column 1 and 4
    @unroll for i = 1 : M
      r1[i,1] = FT(0)  
      r1[i,2] = FT(0)  
    end  
    r1[1,1] = B1_1[iz,ID]
    r1[M,2] = B1_4[iz,ID]
    @unroll for i = 1 : M - 2
      r3[i,1] = B3_14[i,1,iz,ID] - invfac * A31[i,1,iz,ID] * r1[1,1]
      r3[i,2] = B3_14[i,2,iz,ID] - invfac * A31[i,M,iz,ID] * r1[M,2]
    end  
  

    @unroll for k = 1 : M - 3
      @unroll for i = k + 1 : M - 2
       r3[i,1] -= SA[i,k,iz,ID] * r3[k,1]
       r3[i,2] -= SA[i,k,iz,ID] * r3[k,2]
      end
    end
#   Backward loop
    @unroll for k = M - 2 : -1 : 1
      r3[k,1] /= SA[k,k,iz,ID]
      r3[k,2] /= SA[k,k,iz,ID]
      @unroll for i = 1 : k - 1
        r3[i,1] -= SA[i,k,iz,ID] * r3[k,1]
        r3[i,2] -= SA[i,k,iz,ID] * r3[k,2]
      end
    end

    @unroll for i = 1 : M
      @unroll for j = 1 : M - 2
        r1[i,1] -= A13[i,j,iz,ID] * r3[j,1]
        r1[i,2] -= A13[i,j,iz,ID] * r3[j,2]
      end
      r1[i,1] *= invfac
      r1[i,2] *= invfac
    end  

    @unroll for i = 1 : M - 2
      r2[i,1] = FT(0)
      r2[i,2] = FT(0)
      @unroll for j = 1 : M - 2
        r2[i,1] -= A23[i,j,iz,ID] * r3[j,1]
        r2[i,2] -= A23[i,j,iz,ID] * r3[j,2]
      end
      r2[i,1] *= invfac
      r2[i,2] *= invfac
    end

    s[1,1] = FT(0)
    s[2,1] = FT(0)
    s[3,1] = FT(0)
    s[4,1] = FT(0)
    s[1,2] = FT(0)
    s[2,2] = FT(0)
    s[3,2] = FT(0)
    s[4,2] = FT(0)
    @unroll for j = 1 : M - 2
      s[1,1] += -C14_3[1,j,iz,ID] * r3[j,1]
      s[1,2] += -C14_3[1,j,iz,ID] * r3[j,2]
      s[4,1] += -C14_3[2,j,iz,ID] * r3[j,1]
      s[4,2] += -C14_3[2,j,iz,ID] * r3[j,2]
      s[2,1] += -C23_2[1,j,iz,ID] * r2[j,1]  
      s[2,2] += -C23_2[1,j,iz,ID] * r2[j,2]
      s[3,1] += -C23_2[2,j,iz,ID] * r2[j,1]
      s[3,2] += -C23_2[2,j,iz,ID] * r2[j,2]
    end
    @unroll for j = 1 : M
      s[2,1] += -C23_1[1,j,iz,ID] * r1[j,1]  
      s[2,2] += -C23_1[1,j,iz,ID] * r1[j,2]
      s[3,1] += -C23_1[2,j,iz,ID] * r1[j,1]
      s[3,2] += -C23_1[2,j,iz,ID] * r1[j,2]
    end
    DD[ID,1,1,iz] += s[1,1]
    DD[ID,2,1,iz] += s[2,1]
    DD[ID,3,1,iz] += s[3,1]
    DD[ID,4,1,iz] += s[4,1]
    DD[ID,1,4,iz] += s[1,2]
    DD[ID,2,4,iz] += s[2,2]
    DD[ID,3,4,iz] += s[3,2]
    DD[ID,4,4,iz] += s[4,2]
  
    @unroll for i = 1 : M
      r1[i,1] = B1_23[i,1,iz,ID]
      r1[i,2] = B1_23[i,2,iz,ID]
    end  
    @unroll for i = 1 : M - 2
      r2[i,1] = B2_23[i,1,iz,ID]
      r2[i,2] = B2_23[i,2,iz,ID]
    end  
    @unroll for i = 1 : M - 2
      r3[i,1] = FT(0)
      r3[i,2] = FT(0)
      @unroll for j = 1 : M 
        r3[i,1] -= A31[i,j,iz,ID] * r1[j,1]
        r3[i,2] -= A31[i,j,iz,ID] * r1[j,2]
      end
      @unroll for j = 1 : M - 2
        r3[i,1] -= A32[i,j,iz,ID] * r2[j,1]
        r3[i,2] -= A32[i,j,iz,ID] * r2[j,2]
      end
      r3[i,1] *= invfac
      r3[i,2] *= invfac
    end

    @unroll for k = 1 : M - 3
      @unroll for i = k + 1 : M - 2
       r3[i,1] -= SA[i,k,iz,ID] * r3[k,1]
       r3[i,2] -= SA[i,k,iz,ID] * r3[k,2]
      end
    end
#   Backward loop
    @unroll for k = M - 2 : -1 : 1
      r3[k,1] /= SA[k,k,iz,ID]
      r3[k,2] /= SA[k,k,iz,ID]
      @unroll for i = 1 : k - 1
        r3[i,1] -= SA[i,k,iz,ID] * r3[k,1]
        r3[i,2] -= SA[i,k,iz,ID] * r3[k,2]
      end
    end

    @unroll for i = 1 : M 
      @unroll for j = 1 : M - 2
        r1[i,1] += -A13[i,j,iz,ID] * r3[j,1]
        r1[i,2] += -A13[i,j,iz,ID] * r3[j,2]
      end
      r1[i,1] *= invfac
      r1[i,2] *= invfac
    end
    @unroll for i = 1 : M - 2
      @unroll for j = 1 : M - 2
        r2[i,1] += -A23[i,j,iz,ID] * r3[j,1]
        r2[i,2] += -A23[i,j,iz,ID] * r3[j,2]
      end
      r2[i,1] *= invfac
      r2[i,2] *= invfac
    end

    s[1,1] = FT(0)
    s[2,1] = FT(0)
    s[3,1] = FT(0)
    s[4,1] = FT(0)
    s[1,2] = FT(0)
    s[2,2] = FT(0)
    s[3,2] = FT(0)
    s[4,2] = FT(0)
    @unroll for j = 1 : M - 2
      s[1,1] += - C14_3[1,j,iz,ID] * r3[j,1]
      s[1,2] += - C14_3[1,j,iz,ID] * r3[j,2]
      s[4,1] += - C14_3[2,j,iz,ID] * r3[j,1]
      s[4,2] += - C14_3[2,j,iz,ID] * r3[j,2]
      s[2,1] += - C23_2[1,j,iz,ID] * r2[j,1]
      s[2,2] += - C23_2[1,j,iz,ID] * r2[j,2]
      s[3,1] += - C23_2[2,j,iz,ID] * r2[j,1]
      s[3,2] += - C23_2[2,j,iz,ID] * r2[j,2]
    end
    @unroll for j = 1 : M
      s[2,1] += -C23_1[1,j,iz,ID] * r1[j,1]  
      s[2,2] += -C23_1[1,j,iz,ID] * r1[j,2]
      s[3,1] += -C23_1[2,j,iz,ID] * r1[j,1]
      s[3,2] += -C23_1[2,j,iz,ID] * r1[j,2]
    end
    DD[ID,1,2,iz] += s[1,1]
    DD[ID,2,2,iz] += s[2,1]
    DD[ID,3,2,iz] += s[3,1]
    DD[ID,4,2,iz] += s[4,1]
    DD[ID,1,3,iz] += s[1,2]
    DD[ID,2,3,iz] += s[2,2]
    DD[ID,3,3,iz] += s[3,2]
    DD[ID,4,3,iz] += s[4,2]

    if iz > 1
#     @atomic :monotonic SchurBand[7,sh-1,ID] -= FacGrav * invfac * B1m_34[1,iz,ID]
#     @atomic :monotonic SchurBand[6,sh,ID] -= FacGrav * invfac * B1m_34[2,iz,ID]
      @unroll for i = 1 : M
        r1[i,1] = FT(0)
        r1[i,2] = FT(0)
      end
      r1[1,1] = B1m_34[1,iz,ID]
      r1[1,2] = B1m_34[2,iz,ID]
      @unroll for i = 1 : M - 2
        r3[i,1] = - invfac * A31[i,1,iz,ID] * r1[1,1]
        r3[i,2] = - invfac * A31[i,1,iz,ID] * r1[1,2]
      end

      @unroll for k = 1 : M - 3
        @unroll for i = k + 1 : M - 2
          r3[i,1] -= SA[i,k,iz,ID] * r3[k,1]
          r3[i,2] -= SA[i,k,iz,ID] * r3[k,2]
        end
      end
#     Backward loop
      @unroll for k = M - 2 : -1 : 1
        r3[k,1] /= SA[k,k,iz,ID]
        r3[k,2] /= SA[k,k,iz,ID]
        @unroll for i = 1 : k - 1
          r3[i,1] -= SA[i,k,iz,ID] * r3[k,1]
          r3[i,2] -= SA[i,k,iz,ID] * r3[k,2]
        end
      end
      @unroll for i = 1 : M
        @unroll for j = 1 : M - 2
          r1[i,1] += -A13[i,j,iz,ID] * r3[j,1]
          r1[i,2] += -A13[i,j,iz,ID] * r3[j,2]
        end
        r1[i,1] *= invfac
        r1[i,2] *= invfac
      end
      @unroll for i = 1 : M - 2
        r2[i,1] = FT(0) 
        r2[i,2] = FT(0) 
        @unroll for j = 1 : M - 2
          r2[i,1] += -A23[i,j,iz,ID] * r3[j,1]
          r2[i,2] += -A23[i,j,iz,ID] * r3[j,2]
        end
        r2[i,1] *= invfac
        r2[i,2] *= invfac
      end

      s[1,1] = FT(0)
      s[2,1] = FT(0)
      s[3,1] = FT(0)
      s[4,1] = FT(0)
      s[1,2] = FT(0)
      s[2,2] = FT(0)
      s[3,2] = FT(0)
      s[4,2] = FT(0)
      @unroll for j = 1 : M - 2
        s[1,1] += - C14_3[1,j,iz,ID] * r3[j,1]
        s[1,2] += - C14_3[1,j,iz,ID] * r3[j,2]
        s[4,1] += - C14_3[2,j,iz,ID] * r3[j,1]
        s[4,2] += - C14_3[2,j,iz,ID] * r3[j,2]
        s[2,1] += - C23_2[1,j,iz,ID] * r2[j,1]
        s[2,2] += - C23_2[1,j,iz,ID] * r2[j,2]
        s[3,1] += - C23_2[2,j,iz,ID] * r2[j,1]
        s[3,2] += - C23_2[2,j,iz,ID] * r2[j,2]
      end
      @unroll for j = 1 : M
        s[2,1] += -C23_1[1,j,iz,ID] * r1[j,1]
        s[2,2] += -C23_1[1,j,iz,ID] * r1[j,2]
        s[3,1] += -C23_1[2,j,iz,ID] * r1[j,1]
        s[3,2] += -C23_1[2,j,iz,ID] * r1[j,2]
      end
      DL[ID,1,1,iz-1] += s[1,1]
      DL[ID,2,1,iz-1] += s[2,1]
      DL[ID,3,1,iz-1] += s[3,1]
      DL[ID,4,1,iz-1] += s[4,1]
      DL[ID,1,2,iz-1] += s[1,2]
      DL[ID,2,2,iz-1] += s[2,2]
      DL[ID,3,2,iz-1] += s[3,2]
      DL[ID,4,2,iz-1] += s[4,2]
    end
    if iz < nz
      #@atomic :monotonic SchurBand[2,sh+5,ID] -= FacGrav * invfac * B1p_12[1,iz,ID]
      #@atomic :monotonic SchurBand[1,sh+6,ID] -= FacGrav * invfac * B1p_12[2,iz,ID]
      @unroll for i = 1 : M
        r1[i,1] = FT(0)
        r1[i,2] = FT(0)
      end
      r1[M,1] = B1p_12[1,iz,ID]
      r1[M,2] = B1p_12[2,iz,ID]
      @unroll for i = 1 : M - 2
        r3[i,1] = - invfac * A31[i,M,iz,ID] * r1[M,1]
        r3[i,2] = - invfac * A31[i,M,iz,ID] * r1[M,2]
      end

      @unroll for k = 1 : M - 3
        @unroll for i = k + 1 : M - 2
          r3[i,1] -= SA[i,k,iz,ID] * r3[k,1]
          r3[i,2] -= SA[i,k,iz,ID] * r3[k,2]
        end
      end
#     Backward loop
      @unroll for k = M - 2 : -1 : 1
        r3[k,1] /= SA[k,k,iz,ID]
        r3[k,2] /= SA[k,k,iz,ID]
        @unroll for i = 1 : k - 1
          r3[i,1] -= SA[i,k,iz,ID] * r3[k,1]
          r3[i,2] -= SA[i,k,iz,ID] * r3[k,2]
        end
      end
      @unroll for i = 1 : M
        @unroll for j = 1 : M - 2
          r1[i,1] += -A13[i,j,iz,ID] * r3[j,1]
          r1[i,2] += -A13[i,j,iz,ID] * r3[j,2]
        end
        r1[i,1] *= invfac
        r1[i,2] *= invfac
      end
      @unroll for i = 1 : M - 2
        r2[i,1] = FT(0) 
        r2[i,2] = FT(0) 
        @unroll for j = 1 : M - 2
          r2[i,1] += -A23[i,j,iz,ID] * r3[j,1]
          r2[i,2] += -A23[i,j,iz,ID] * r3[j,2]
        end
        r2[i,1] *= invfac
        r2[i,2] *= invfac
      end

      s[1,1] = FT(0)
      s[2,1] = FT(0)
      s[3,1] = FT(0)
      s[4,1] = FT(0)
      s[1,2] = FT(0)
      s[2,2] = FT(0)
      s[3,2] = FT(0)
      s[4,2] = FT(0)
      @unroll for j = 1 : M - 2
        s[1,1] += - C14_3[1,j,iz,ID] * r3[j,1]
        s[1,2] += - C14_3[1,j,iz,ID] * r3[j,2]
        s[4,1] += - C14_3[2,j,iz,ID] * r3[j,1]
        s[4,2] += - C14_3[2,j,iz,ID] * r3[j,2]
        s[2,1] += - C23_2[1,j,iz,ID] * r2[j,1]
        s[2,2] += - C23_2[1,j,iz,ID] * r2[j,2]
        s[3,1] += - C23_2[2,j,iz,ID] * r2[j,1]
        s[3,2] += - C23_2[2,j,iz,ID] * r2[j,2]
      end
      @unroll for j = 1 : M
        s[2,1] += -C23_1[1,j,iz,ID] * r1[j,1]
        s[2,2] += -C23_1[1,j,iz,ID] * r1[j,2]
        s[3,1] += -C23_1[2,j,iz,ID] * r1[j,1]
        s[3,2] += -C23_1[2,j,iz,ID] * r1[j,2]
      end
      DU[ID,1,1,iz] += s[1,1]
      DU[ID,2,1,iz] += s[2,1]
      DU[ID,3,1,iz] += s[3,1]
      DU[ID,4,1,iz] += s[4,1]
      DU[ID,1,2,iz] += s[1,2]
      DU[ID,2,2,iz] += s[2,2]
      DU[ID,3,2,iz] += s[3,2]
      DU[ID,4,2,iz] += s[4,2]
    end
  end
end      

@kernel inbounds = true function ldivSplitVerticalFKernel!(@Const(A13),@Const(A23),@Const(A31),@Const(A32),
  @Const(B1m_34),@Const(B1_1),@Const(B1_23),@Const(B1_4), @Const(B2_23),@Const(B3_14),
  @Const(B1p_12),@Const(C23_1),@Const(C23_2),@Const(C14_3),@Const(SA),
  b,rs,invfac, ::Val{M}) where {M}

  iz,ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  ND = @uniform @ndrange()[2]

  @uniform FT = eltype(SA)

  r1 = @private FT (M,)
  r2 = @private FT (M-2,)
  r3 = @private FT (M-2,)
  rsS = @private FT (4,)

  @uniform RhoPos = 1
  @uniform wPos = 4
  @uniform ThPos = 5

  if ID <= ND
    sh = (iz - 1) * 4
    rsS[1] = b[1,iz,ID,ThPos]
    rsS[2] = b[1,iz,ID,wPos]
    rsS[3] = b[M,iz,ID,wPos]
    rsS[4] = b[M,iz,ID,ThPos]
    @unroll for i = 1 : M
      r1[i] = b[i,iz,ID,RhoPos]
    end  
    @unroll for i = 1 : M - 2
      r2[i] = b[i+1,iz,ID,ThPos]
      r3[i] = b[i+1,iz,ID,wPos]
    end  
    @unroll for i = 1 : M - 2
      for j = 1 : M
        r3[i] += -invfac * A31[i,j,iz,ID] * r1[j]
      end  
      @unroll for j = 1 : M - 2
        r3[i] += -invfac * A32[i,j,iz,ID] * r2[j]
      end
    end  

    ldivFull!(iz,ID,SA,r3,Val(M - 2))

    @unroll for i = 1 : M 
      @unroll for j = 1 : M - 2
        r1[i] += -A13[i,j,iz,ID] * r3[j]
      end
      r1[i] *= invfac
    end

    @unroll for i = 1 : M - 2
      @unroll for j = 1 : M - 2
        r2[i] += -A23[i,j,iz,ID] * r3[j]
      end
      r2[i] *= invfac
    end
    @unroll for j = 1 : M - 2
      rsS[1] += -C14_3[1,j,iz,ID] * r3[j]
      rsS[4] += -C14_3[2,j,iz,ID] * r3[j]
      rsS[2] += -C23_2[1,j,iz,ID] * r2[j]
      rsS[3] += -C23_2[2,j,iz,ID] * r2[j]
    end
    @unroll for j = 1 : M 
      rsS[2] += -C23_1[1,j,iz,ID] * r1[j]
      rsS[3] += -C23_1[2,j,iz,ID] * r1[j]
    end
    rs[ID,1,iz] = rsS[1]
    rs[ID,2,iz] = rsS[2]
    rs[ID,3,iz] = rsS[3]
    rs[ID,4,iz] = rsS[4]
  end
end

@kernel inbounds = true function BlockThomasKernel!(
    @Const(DL),         # Untere Diagonale (gepackt 4x2):   (num_systems, 4, 2, nz-1)
    @Const(DD),         # Hauptdiagonale (4x4):            (num_systems, 4, 4, nz)
    @Const(DU),         # Obere Diagonale (gepackt 4x2):    (num_systems, 4, 2, nz-1)
    D,                  # Eingang & Ausgang (RHS / Lösung): (num_systems, 4, nz)
    ::Val{nz}, 
) where {nz} 
    iD = @index(Local, Linear)
    ID = @index(Global, Linear)
    NDloc = @uniform @groupsize()[1]
    ND = @uniform @ndrange()[1]

    
    # Lokaler On-Chip-Speicher via @localmem
    C_local = @localmem eltype(DD) (NDloc, 4, 4, nz)
    D_local = @localmem eltype(DD) (NDloc, 4, nz)

    
    if ID <= ND
        
        # --- 1. VORWÄRTSELIMINATION (i = 1) ---
        B_1 = SMatrix{4, 4, eltype(DD)}(
            DD[ID, 1, 1, 1], DD[ID, 2, 1, 1], DD[ID, 3, 1, 1], DD[ID, 4, 1, 1],
            DD[ID, 1, 2, 1], DD[ID, 2, 2, 1], DD[ID, 3, 2, 1], DD[ID, 4, 2, 1],
            DD[ID, 1, 3, 1], DD[ID, 2, 3, 1], DD[ID, 3, 3, 1], DD[ID, 4, 3, 1],
            DD[ID, 1, 4, 1], DD[ID, 2, 4, 1], DD[ID, 3, 4, 1], DD[ID, 4, 4, 1]
        )
        B_inv = inv(B_1)
        
        # DU hat nur N-1 Blöcke. Block 1 existiert immer, da N >= 2 vorausgesetzt wird.
        # Spalten 1-2 besetzt, Spalten 3-4 sind Null (Padding)
        C_1 = SMatrix{4, 4, eltype(DD)}(
            DU[ID, 1, 1, 1], DU[ID, 2, 1, 1], DU[ID, 3, 1, 1], DU[ID, 4, 1, 1],
            DU[ID, 1, 2, 1], DU[ID, 2, 2, 1], DU[ID, 3, 2, 1], DU[ID, 4, 2, 1],
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0
        )
        
        C_prime = B_inv * C_1
        D_1 = SVector{4, eltype(DD)}(D[ID, 1, 1], D[ID, 2, 1], D[ID, 3, 1], D[ID, 4, 1])
        D_prime = B_inv * D_1
        
        for j in 1:4, i in 1:4 C_local[iD, i, j, 1] = C_prime[i, j] end
        for i in 1:4 D_local[iD, i, 1] = D_prime[i] end
        
        # --- Schleife für i = 2 bis N ---
        for iz in 2:nz
            B_i = SMatrix{4, 4, eltype(DD)}(
                DD[ID, 1, 1, iz], DD[ID, 2, 1, iz], DD[ID, 3, 1, iz], DD[ID, 4, 1, iz],
                DD[ID, 1, 2, iz], DD[ID, 2, 2, iz], DD[ID, 3, 2, iz], DD[ID, 4, 2, iz],
                DD[ID, 1, 3, iz], DD[ID, 2, 3, iz], DD[ID, 3, 3, iz], DD[ID, 4, 3, iz],
                DD[ID, 1, 4, iz], DD[ID, 2, 4, iz], DD[ID, 3, 4, iz], DD[ID, 4, 4, iz]
            )
            D_i = SVector{4, eltype(DD)}(D[ID, 1, iz], D[ID, 2, iz], D[ID, 3, iz], D[ID, 4, iz])
            
            # DL_idx greift auf iz - 1 zu (läuft von 1 bis N-1)
            izm1 = iz - 1
            
            # A_i aufbauen aus DL (Spalten 1-2 sind 0, Spalten 3-4 gepackt)
            A_i = SMatrix{4, 4, eltype(DD)}(
                0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0,
                DL[ID, 1, 1, izm1], DL[ID, 2, 1, izm1], DL[ID, 3, 1, izm1], DL[ID, 4, 1, izm1],
                DL[ID, 1, 2, izm1], DL[ID, 2, 2, izm1], DL[ID, 3, 2, izm1], DL[ID, 4, 2, izm1]
            )
            
            B_mod_inv = inv(B_i - A_i * C_prime)
            D_prime = B_mod_inv * (D_i - A_i * D_prime)
            for i in 1:4 D_local[iD, i, iz] = D_prime[i] end
            
            # C_i / DU_i existiert nur bis einschließlich iz == N-1
            if iz < nz
                C_i = SMatrix{4, 4, eltype(DD)}(
                    DU[ID, 1, 1, iz], DU[ID, 2, 1, iz], DU[ID, 3, 1, iz], DU[ID, 4, 1, iz],
                    DU[ID, 1, 2, iz], DU[ID, 2, 2, iz], DU[ID, 3, 2, iz], DU[ID, 4, 2, iz],
                    0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0
                )
                C_prime = B_mod_inv * C_i
                for j in 1:4, i in 1:4 C_local[iD, i, j, iz] = C_prime[i, j] end
            end
        end
        
        # --- 2. RÜCKWÄRTSSUBSTITUTION (In-Place in D) ---
        X_curr = SVector{4, eltype(DD)}(D_local[iD, 1, nz], D_local[iD, 2, nz], D_local[iD, 3, nz], D_local[iD, 4, nz])
        for i in 1:4 D[ID, i, nz] = X_curr[i] end  
        
        for iz in (nz-1):-1:1
            C_prime = SMatrix{4, 4, eltype(DD)}(
                C_local[iD, 1, 1, iz], C_local[iD, 2, 1, iz], C_local[iD, 3, 1, iz], C_local[iD, 4, 1, iz],
                C_local[iD, 1, 2, iz], C_local[iD, 2, 2, iz], C_local[iD, 3, 2, iz], C_local[iD, 4, 2, iz],
                C_local[iD, 1, 3, iz], C_local[iD, 2, 3, iz], C_local[iD, 3, 3, iz], C_local[iD, 4, 3, iz],
                C_local[iD, 1, 4, iz], C_local[iD, 2, 4, iz], C_local[iD, 3, 4, iz], C_local[iD, 4, 4, iz]
            )
            D_prime = SVector{4, eltype(DD)}(D_local[iD, 1, iz], D_local[iD, 2, iz], D_local[iD, 3, iz], D_local[iD, 4, iz])
            X_curr = D_prime - C_prime * X_curr
            for i in 1:4 D[ID, i, iz] = X_curr[i] end  
        end
    end
end

@kernel inbounds = true function ldivSplitVerticalBKernel!(@Const(A13),@Const(A23),@Const(A31),@Const(A32),
  @Const(B1m_34),@Const(B1_1),@Const(B1_23),@Const(B1_4), @Const(B2_23),@Const(B3_14),
  @Const(B1p_12),@Const(SA),b,rs,invfac,  ::Val{M}) where {M}

  iz,ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  ND = @uniform @ndrange()[2]

  @uniform FT = eltype(SA)

  r1 = @private FT (M,)
  r2 = @private FT (M-2,)
  r3 = @private FT (M-2,)

  @uniform RhoPos = 1
  @uniform wPos = 4
  @uniform ThPos = 5

  if ID <= ND

    b[1,iz,ID,ThPos] = rs[ID,1,iz]
    b[1,iz,ID,wPos] = rs[ID,2,iz]
    b[M,iz,ID,wPos] = rs[ID,3,iz]
    b[M,iz,ID,ThPos] = rs[ID,4,iz]
    @unroll for i = 1 : M
      r1[i] = b[i,iz,ID,RhoPos]
    end  
    @unroll for i = 1 : M - 2
      r2[i] = b[i+1,iz,ID,ThPos]
      r3[i] = b[i+1,iz,ID,wPos]
    end  
    r1[1] -= B1_1[iz,ID] * rs[ID,1,iz]
    r1[M] -= B1_4[iz,ID] * rs[ID,4,iz]
    @unroll for i = 1 : M
      r1[i] -= B1_23[i,1,iz,ID] * rs[ID,2,iz] + B1_23[i,2,iz,ID] * rs[ID,3,iz]
    end
    @unroll for i = 1 : M - 2
      r2[i] -= B2_23[i,1,iz,ID] * rs[ID,2,iz] + B2_23[i,2,iz,ID] * rs[ID,3,iz]
      r3[i] -= B3_14[i,1,iz,ID] * rs[ID,1,iz] + B3_14[i,2,iz,ID] * rs[ID,4,iz]
    end
    if iz > 1
      izm1 = iz - 1  
      r1[1] -= B1m_34[1,iz,ID] * rs[ID,3,izm1] + B1m_34[2,iz,ID] * rs[ID,4,izm1]
    end
    if iz < nz
      izp1 = iz + 1  
      r1[M] -= B1p_12[1,iz,ID] * rs[ID,1,izp1] + B1p_12[2,iz,ID] * rs[ID,2,izp1]
    end
    @unroll for i = 1 : M - 2
      @unroll for j = 1 : M 
        r3[i] += -invfac * A31[i,j,iz,ID] * r1[j]
      end
    end
    @unroll for i = 1 : M - 2
      @unroll for j = 1 : M - 2
        r3[i] += -invfac * A32[i,j,iz,ID] * r2[j]
      end
    end
    ldivFull!(iz,ID,SA,r3,Val(M - 2))
    @unroll for i = 1 : M
      @unroll for j = 1 : M - 2
        r1[i] += -A13[i,j,iz,ID] * r3[j]
      end
      r1[i] *= invfac
    end
    @unroll for i = 1 : M - 2
      @unroll for j = 1 : M - 2
        r2[i] += -A23[i,j,iz,ID] * r3[j]
      end
      r2[i] *= invfac
    end
    @unroll for i = 1 : M
      b[i,iz,ID,RhoPos] = r1[i]
    end
    @unroll for i = 1 : M -2 
      b[i+1,iz,ID,ThPos] = r2[i]
      b[i+1,iz,ID,wPos] = r3[i]
    end  
  end  
end  

function FillJacSplitDGVert!(Jac::JacSplitDGVert,U,DG,dz,fac,Phys)
  
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
  KFillJacSplitDGVertKernel! = FillJacSplitDGVertKernel!(backend,group)
  KFillJacSplitDGVertKernel!(Jac.A13,Jac.A23,Jac.A31,Jac.A32,Jac.B1m_34,Jac.B1_1,Jac.B1_23,Jac.B1_4,
  Jac.B2_23,Jac.B3_14,Jac.B1p_12,Jac.C23_1,Jac.C23_2,Jac.C14_3,Jac.SA,Jac.DL,Jac.DD,Jac.DU,U,dz,
  DWZ,DG.wZ,fac,Phys.cS,Phys,Val(M);ndrange=ndrange)

end   

function SchurBoundary!(Jac::JacSplitDGVert)

  backend = get_backend(Jac.DD)
  FTB = eltype(Jac.DD)

  M = Jac.M
  nz = Jac.nz
  ND = size(Jac.DD,1)

  NDG = 1
  group = (nz, NDG)
  ndrange = (nz, ND)
  KSchurBoundaryKernel! = SchurBoundaryKernel!(backend,group)
  KSchurBoundaryKernel!(Jac.A13,Jac.A23,Jac.A31,Jac.A32,Jac.B1m_34,Jac.B1_1,
    Jac.B1_23,Jac.B1_4, Jac.B2_23,Jac.B3_14,Jac.B1p_12,Jac.C23_1,Jac.C23_2,
    Jac.C14_3,Jac.SA,Jac.DU,Jac.DD,Jac.DL,Jac.fac,Val(M);ndrange=ndrange)
end  

function SchurBoundary1!(Jac::JacSplitDGVert)

  backend = get_backend(Jac.DD)
  FTB = eltype(Jac.DD)

  M = Jac.M
  nz = Jac.nz
  ND = size(Jac.DD,1)

  NDG = 1
  group = (nz, NDG)
  ndrange = (nz, ND)
  KSchurBoundaryKernel! = SchurBoundaryKernel1!(backend,group)
  KSchurBoundaryKernel!(Jac.A13,Jac.A23,Jac.A31,Jac.A32,Jac.B1m_34,Jac.B1_1,
    Jac.B1_23,Jac.B1_4, Jac.B2_23,Jac.B3_14,Jac.B1p_12,Jac.C23_1,Jac.C23_2,
    Jac.C14_3,Jac.SA,Jac.DU,Jac.DD,Jac.DL,Jac.fac,Val(M);ndrange=ndrange)
end  


function Solve!(Jac::JacSplitDGVert,b)

  backend = get_backend(Jac.DD)
  FTB = eltype(Jac.DD)

  M = Jac.M
  nz = Jac.nz
  ND = size(Jac.DD,1)

  invfac = FTB(1) / Jac.fac

  NDG = 32
  group = (nz, NDG)
  ndrange = (nz, ND)
  KldivVerticalFKernel! = ldivSplitVerticalFKernel!(backend,group)
  KldivVerticalFKernel!(Jac.A13,Jac.A23,Jac.A31,Jac.A32,Jac.B1m_34,Jac.B1_1,
    Jac.B1_23,Jac.B1_4, Jac.B2_23,Jac.B3_14,Jac.B1p_12,Jac.C23_1,Jac.C23_2,
    Jac.C14_3,Jac.SA,b,Jac.rs,invfac,Val(M);ndrange=ndrange)

  group = (NDG)
  ndrange = (ND)
  KldivVerticalSKernel! = BlockThomasKernel!(backend,group)
  KldivVerticalSKernel!(Jac.DL,Jac.DD,Jac.DU,Jac.rs,Val(nz);ndrange=ndrange)

  group = (nz, NDG)
  ndrange = (nz, ND)
  KldivVerticalBKernel! = ldivSplitVerticalBKernel!(backend,group)
  KldivVerticalBKernel!(Jac.A13,Jac.A23,Jac.A31,Jac.A32,Jac.B1m_34,Jac.B1_1,
    Jac.B1_23,Jac.B1_4, Jac.B2_23,Jac.B3_14,Jac.B1p_12,
    Jac.SA,b,Jac.rs,invfac,Val(M);ndrange=ndrange)

end

function Jac!(U,fac,DG,Metric,Phys,Cache,JCache::JacSplitDGVert,Global,VelForm)
  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU
  if JCache.grav_do
    @views Geo = Cache.Aux[:,:,:,2]
    precompute_gravity!(Geo,Metric.dz,DG.DWZ,JCache,NumberThreadGPU)
    JCache.grav_do = false
  end  
  Invfac = eltype(U)(1) / fac 
  dz = Metric.dz
  FillJacSplitDGVert!(JCache,U,DG,dz,Invfac,Phys)
  SchurBoundary!(JCache)
end

function Solve!(k,v,Jac::JacSplitDGVert,fac,DG::FiniteElements.DGElement,Metric,Global,VelForm)

  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU
  @. k = v
  @views TendVCart2VSp!(k,DG,Metric,NumberThreadGPU,VelForm)
  Solve!(Jac,k)
  @views @. k[:,:,:,2:3] *= fac
  @views TendVSp2VCart!(k,DG,Metric,NumberThreadGPU,VelForm)
end

function precompute_gravity!(GeoPot,dz, DWZ,cache, NumberThreadGPU)
  (; A31, C23_1) = cache
  backend = get_backend(dz)
  M = size(DWZ,1)
  nz = size(A31,3)
  NumI = size(A31,4)
  NumIG = min(div(NumberThreadGPU,M*nz),NumI)
  group = (M,nz,NumIG)
  ndrange = (M,nz,NumI)
  Kprecompute_gravityKernel! = precompute_gravityKernel!(backend,group)
  Kprecompute_gravityKernel!(A31,C23_1,GeoPot,DWZ,dz,Val(M),ndrange=ndrange)
end

@kernel inbounds = true function precompute_gravityKernel!(A31,C23_1,@Const(Geo),@Const(DS),@Const(dz), ::Val{M}) where {M}

  i,iz,ID = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]
  if ID <= ND
    acc = eltype(dz)(0)
    Geoi = Geo[i, iz, ID]
    inv2dz = eltype(A31)(2) / dz[iz, ID]
    if i == 1
      @unroll for k in 1:M
        dphi = Geo[k, iz, ID] - Geoi
        C23_1[1,k,iz,ID] = inv2dz * 0.5 * DS[i, k] * dphi
        acc += DS[i, k] * dphi
      end
      C23_1[1, 1, iz, ID] += inv2dz * 0.5 * acc
    elseif i == M  
      @unroll for k in 1:M
        dphi = Geo[k, iz, ID] - Geoi
        C23_1[2,k,iz,ID] = inv2dz * 0.5 * DS[i, k] * dphi
        acc += DS[i, k] * dphi
      end
      C23_1[2, M, iz, ID] += inv2dz * 0.5 * acc
    else  
      @unroll for k in 1:M
        dphi = Geo[k, iz, ID] - Geoi
        A31[i-1,k,iz,ID] = inv2dz * 0.5 * DS[i, k] * dphi
        acc += DS[i, k] * dphi
      end
      A31[i-1, i, iz, ID] += inv2dz * 0.5 * acc
    end  
  end
end
