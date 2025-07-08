mutable struct JacDGVert{FT<:AbstractFloat,
                         AT2<:AbstractArray,
                         AT3<:AbstractArray,
                         AT4<:AbstractArray}
  M::Int
  nz::Int
  fac::FT
  FacGrav::FT
  A13::AT4
  A23::AT4
  A32::AT4
  B1m_34::AT3
  B1_1::AT2
  B1_23::AT4
  B1_4::AT2
  B2_23::AT4
  B3_14::AT4
  B1p_12::AT3
  C23_2::AT4
  C14_3::AT4
  SA::AT4
  SchurBand::AT3
  rs::AT2
end  

function LowOrder(DG)

  M = size(DG.xwZ,1)
  D = zeros(M,M)
  D[1,1] = 0.5
  D[1,2] = 0.5
  for i = 2 : M - 1
    # cI = c  
    D[i,i-1] = -0.5
    D[i,i+1] = 0.5
  end
  D[M,M-1] = -0.5
  D[M,M] = -0.5
  DW = inv(diagm(DG.wZ)) * D
end

@inline function LUFull!(A)

  n = size(A,1)
  @inbounds for j = 1 : n - 1
    @inbounds for i = j + 1 : n 
      A[i,j] /= A[j,j]  
      @inbounds for k = j + 1 : n
        A[i,k] -= A[i,j] * A[j,k]
      end  
    end  
  end  
end

@inline function ldivFull2!(A,b)

# Forward loop
  n = size(A,1)
  @inbounds for k = 1 : n - 1
    @inbounds for i = k + 1 : n
      @views @.  b[i,:] -= A[i,k] * b[k,:]
    end
  end
#  Backward loop
  @inbounds for k = n : -1 : 1
    @views @.  b[k,:] /= A[k,k]
    @inbounds for i = 1 : k - 1
      @views @.  b[i,:] -= A[i,k] * b[k,:]
    end
  end
end

@inline function ldivFull!(A,b)

# Forward loop
  n = size(A,1)
  @inbounds for k = 1 : n - 1
    @inbounds for i = k + 1 : n
      @views b[i] -= A[i,k] * b[k]
    end
  end
#  Backward loop
  @inbounds for k = n : -1 : 1
    @views b[k] /= A[k,k]
    @inbounds for i = 1 : k - 1
      @views b[i] -= A[i,k] * b[k]
    end
  end
end

@inline function ldivBand!(A,b,l,u)

# Ly = b
  A = A
  n = size(A,2)
  up1 = u + 1 
  uplp1 = up1 + l

  @inbounds for k = 1 : n - 1
    @inbounds for i = k + 1 : min(k+l,n)
      b[i] -= A[i - k + up1,k] * b[k]  
    end    
  end
# Ux = y
  @inbounds for k = n : -1 : 1
    b[k] /= A[up1,k]
    @inbounds for i = max(k - u, 1) : k - 1
      b[i] -= A[up1 + i - k,k] * b[k]  
    end    
  end
end

@inline function luBand!(A,l,u)
#1  *   *   *  a14  ...    
#2  *   *  a13 a24  ... 
#3  *  a12 a23 a34  ... 
#4 a11 a22 a33 a44  ...  
#5 a21 a32 a43 a54  ...  
#6 a31 a42 a53 a64  ...  
#7 a41 a52 a63 a74  ...  

#1  an-6n-3 an-5n-2 an-4n-1 an-3n
#2  an-5n-3 an-4n-2 an-3n-1 an-2n
#3  an-4n-3 an-3n-2 an-2n-1 an-1n
#4  an-3n-3 an-2n-2 an-1n-1 an  n
#5  an-2n-3 an-1n-2 an  n-1   *
#6  an-1n-3 an  n-2   *       *
#7  an  n-3    *      *       *


# a11 a12 a13 a14 ....
# a21 a22 a23 a24 a25 ...
# a31 a32 a33 a34 a35 a36 ...
# a41 a42 a43 a44 a45 a46 a47 ...
# ... a52 a53 a54 a55 a56 a57 a58 ...

#  an-1n-3 an-1n-2 an-1n-1 an-1n 
#    ann-3   ann-2   ann-1   ann

# a21 => (5,1)
# a31 => (6,1)
# a41 => (7,1)

# a12 => (3,2)
# a13 => (2,3) 
# a14 => (1,4)
  n = size(A,2)
  up1 = u + 1
  up2 = u + 2
  uplp1 = u + l + 1
  @inbounds for k = 1 : n -1 
    @inbounds for i = up2 : uplp1 
      A[i,k] /= A[up1,k]  
    end  
    @inbounds for i = up2 : uplp1 
      @inbounds for j = k + 1 : min(k+u,n)
        A[i+k-j,j] -= A[i,k] * A[k + up1 - j,j]
      end
    end  
  end    
end

function JacDGVert{FT}(backend,M,nz,NumG) where FT<:AbstractFloat
  M2 = M - 2
  fac = 0
  FacGrav = 0
  A13 = KernelAbstractions.zeros(backend,FT,M,M2,nz,NumG)
  A23 = KernelAbstractions.zeros(backend,FT,M2,M2,nz,NumG)
  A32 = KernelAbstractions.zeros(backend,FT,M2,M2,nz,NumG)
  B1m_34 = KernelAbstractions.zeros(backend,FT,2,nz,NumG)
  B1_1 = KernelAbstractions.zeros(backend,FT,nz,NumG)
  B1_23 = KernelAbstractions.zeros(backend,FT,M,2,nz,NumG)
  B1_4 = KernelAbstractions.zeros(backend,FT,nz,NumG)
  B2_23 = KernelAbstractions.zeros(backend,FT,M2,2,nz,NumG)
  B3_14 = KernelAbstractions.zeros(backend,FT,M2,2,nz,NumG)
  B1p_12 = KernelAbstractions.zeros(backend,FT,2,nz,NumG)
  C23_2 = KernelAbstractions.zeros(backend,FT,2,M2,nz,NumG)
  C14_3 = KernelAbstractions.zeros(backend,FT,2,M2,nz,NumG)
  SA = KernelAbstractions.zeros(backend,FT,M2,M2,nz,NumG)
  SchurBand = KernelAbstractions.zeros(backend,FT,7,4*nz,NumG)
  rs = KernelAbstractions.zeros(backend,FT,4*nz,NumG)

  return JacDGVert{FT,
                   typeof(B1_1),
                   typeof(B1m_34),
                   typeof(A13)}(

    M,
    nz,
    fac,
    FacGrav,
    A13,
    A23,
    A32,
    B1m_34,
    B1_1,
    B1_23,
    B1_4,
    B2_23,
    B3_14,
    B1p_12,
    C23_2,
    C14_3,
    SA,
    SchurBand,
    rs
  )
end  

@kernel inbounds = true function ldivVerticalFKernel!(@Const(A13),@Const(A23),@Const(A32),
  @Const(B1m_34),@Const(B1_1),@Const(B1_23),@Const(B1_4), @Const(B2_23),@Const(B3_14),
  @Const(B1p_12),@Const(C23_2),@Const(C14_3),@Const(SA),@Const(SchurBand),
  b,rs,invfac,FacGrav, ::Val{M}) where {M}

  iz,ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  DoF = @uniform @ndrange()[2]

  @uniform FT = eltype(SchurBand)

  r1 = @private FT (M,)
  r2 = @private FT (M-2,)
  r3 = @private FT (M-2,)
  rsS = @private FT (4,)

  @uniform RhoPos = 1
  @uniform wPos = 4
  @uniform ThPos = 5

  if ID <= DoF
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
      r3[i] += -invfac * FacGrav * r1[i+1]
      @unroll for j = 1 : M - 2
        r3[i] += -invfac * A32[i,j,iz,ID] * r2[j]
      end
    end  

    @views ldivFull!(SA[:,:,iz,ID],r3)

    r11 = r1[1]
    r1M = r1[M]
    @unroll for j = 1 : M - 2
      r11 += -A13[1,j,iz,ID] * r3[j]
      r1M += -A13[M,j,iz,ID] * r3[j]
    end
    r11 *= invfac
    r1M *= invfac

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
    rsS[2] += -FacGrav * r11
    rsS[3] += -FacGrav * r1M
    rs[sh + 1,ID] = rsS[1]
    rs[sh + 2,ID] = rsS[2]
    rs[sh + 3,ID] = rsS[3]
    rs[sh + 4,ID] = rsS[4]
  end
end

@kernel inbounds = true function ldivVerticalSKernel!(A,rs)

  ID, = @index(Global, NTuple)

  DoF = @uniform @ndrange()[1]

  if ID <= DoF
    @views ldivBand!(A[:,:,ID],rs[:,ID],3,3)
  end  
end  

@kernel inbounds = true function ldivVerticalBKernel!(@Const(A13),@Const(A23),@Const(A32),
  @Const(B1m_34),@Const(B1_1),@Const(B1_23),@Const(B1_4), @Const(B2_23),@Const(B3_14),
  @Const(B1p_12),@Const(SA),@Const(SchurBand),b,rs,invfac,FacGrav,  ::Val{M}) where {M}

  iz,ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  DoF = @uniform @ndrange()[2]

  @uniform FT = eltype(SchurBand)

  r1 = @private FT (M,)
  r2 = @private FT (M-2,)
  r3 = @private FT (M-2,)

  @uniform RhoPos = 1
  @uniform wPos = 4
  @uniform ThPos = 5

  if ID <= DoF
    sh = (iz - 1) * 4

    b[1,iz,ID,ThPos] = rs[sh + 1,ID]
    b[1,iz,ID,wPos] = rs[sh + 2,ID]
    b[M,iz,ID,wPos] = rs[sh + 3,ID]
    b[M,iz,ID,ThPos] = rs[sh + 4,ID]
    @unroll for i = 1 : M
      r1[i] = b[i,iz,ID,RhoPos]
    end  
    @unroll for i = 1 : M - 2
      r2[i] = b[i+1,iz,ID,ThPos]
      r3[i] = b[i+1,iz,ID,wPos]
    end  
    r1[1] -= B1_1[iz,ID] * rs[sh + 1,ID]
    r1[M] -= B1_4[iz,ID] * rs[sh + 4,ID]
    @unroll for i = 1 : M
      r1[i] -= B1_23[i,1,iz,ID] * rs[sh + 2,ID] + B1_23[i,2,iz,ID] * rs[sh + 3,ID]
    end
    @unroll for i = 1 : M - 2
      r2[i] -= B2_23[i,1,iz,ID] * rs[sh + 2,ID] + B2_23[i,2,iz,ID] * rs[sh + 3,ID]
      r3[i] -= B3_14[i,1,iz,ID] * rs[sh + 1,ID] + B3_14[i,2,iz,ID] * rs[sh + 4,ID]
    end
    if iz > 1
      r1[1] -= B1m_34[1,iz,ID] * rs[sh-1,ID] + B1m_34[2,iz,ID] * rs[sh,ID]
    end
    if iz < nz
      r1[M] -= B1p_12[1,iz,ID] * rs[sh+5,ID] + B1p_12[2,iz,ID] * rs[sh+6,ID]
    end
    @unroll for i = 1 : M - 2
      r3[i] += -invfac * FacGrav * r1[i+1]
      @unroll for j = 1 : M - 2
        r3[i] += -invfac * A32[i,j,iz,ID] * r2[j]
      end
    end
    @views ldivFull!(SA[:,:,iz,ID],r3)
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


@kernel inbounds = true function luBandKernel!(A)
  ID, = @index(Global, NTuple)

  DoF = @uniform @ndrange()[1]

  if ID <= DoF
    @views luBand!(A[:,:,ID],3,3)
  end  
end  

@kernel inbounds = true function SchurBoundaryKernel!(@Const(A13),@Const(A23),@Const(A32),
  @Const(B1m_34),@Const(B1_1),@Const(B1_23),@Const(B1_4), @Const(B2_23),@Const(B3_14),
  @Const(B1p_12),@Const(C23_2),@Const(C14_3),@Const(SA),SchurBand,fac,FacGrav, ::Val{M}) where {M}

  iz,ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  DoF = @uniform @ndrange()[2]

  @uniform FT = eltype(SchurBand)

  r1 = @private FT (M,2,)
  r2 = @private FT (M-2,2,)
  r3 = @private FT (M-2,2,)
  r11 = @private FT (2,)
  r1M = @private FT (2,)
  s = @private FT (4,2,)
  invfac = FT(1) / fac


  if ID <= DoF
    sh = (iz - 1) * 4
#   Column 1 and 4
    @unroll for i = 1 : M - 2
      r3[i,1] = B3_14[i,1,iz,ID]
      r3[i,2] = B3_14[i,2,iz,ID]
    end  

    @views ldivFull2!(SA[:,:,iz,ID],r3)

    r11[1] = -B1_1[iz,ID]
    r11[2] = FT(0)
    r1M[1] = FT(0)
    r1M[2] = -B1_4[iz,ID]
    @unroll for j = 1 : M - 2
      r11[1] += A13[1,j,iz,ID] * r3[j,1]
      r11[2] += A13[1,j,iz,ID] * r3[j,2]
      r1M[1] += A13[M,j,iz,ID] * r3[j,1]
      r1M[2] += A13[M,j,iz,ID] * r3[j,2]
    end  
    r11[1] *= -invfac
    r11[2] *= -invfac
    r1M[1] *= -invfac
    r1M[2] *= -invfac

    @unroll for i = 1 : M - 2
      r2[i,1] = FT(0)
      r2[i,2] = FT(0)
      @unroll for j = 1 : M - 2
        r2[i,1] += A23[i,j,iz,ID] * r3[j,1]
        r2[i,2] += A23[i,j,iz,ID] * r3[j,2]
      end
      r2[i,1] *= -invfac
      r2[i,2] *= -invfac
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
    s[2,1] = s[2,1] - FacGrav * r11[1]
    s[3,1] = s[3,1] - FacGrav * r1M[1]
    s[2,2] = s[2,2] - FacGrav * r11[2]
    s[3,2] = s[3,2] - FacGrav * r1M[2]
    @unroll for i = 1 : 4
      @atomic :monotonic SchurBand[i + 3,sh + 1,ID] += s[i,1]
      @atomic :monotonic SchurBand[i,sh + 4,ID] += s[i,2]
    end

    @unroll for i = 1 : M
      r1[i,1] = B1_23[i,1,iz,ID]
      r1[i,2] = B1_23[i,2,iz,ID]
    end  
    @unroll for i = 1 : M - 2
      r2[i,1] = B2_23[i,1,iz,ID]
      r2[i,2] = B2_23[i,2,iz,ID]
    end  
    @unroll for i = 1 : M - 2
      r3[i,1] = FacGrav * r1[i+1,1]
      r3[i,2] = FacGrav * r1[i+1,2]
      @unroll for j = 1 : M - 2
        r3[i,1] += A32[i,j,iz,ID] * r2[j,1]
        r3[i,2] += A32[i,j,iz,ID] * r2[j,2]
      end
      r3[i,1] *= -invfac
      r3[i,2] *= -invfac
    end

    @views ldivFull2!(SA[:,:,iz,ID],r3)

    #r11 = invfac * (r1[1:1,:] - A13[1:1,:,iz,ID] * r3[:,:])
    #r1M = invfac * (r1[M:M,:] - A13[M:M,:,iz,ID] * r3[:,:])
    r11[1] = r1[1,1]
    r11[2] = r1[1,2]
    r1M[1] = r1[M,1]
    r1M[2] = r1[M,2]
    @unroll for j = 1 : M - 2
      r11[1] += -A13[1,j,iz,ID] * r3[j,1]
      r11[2] += -A13[1,j,iz,ID] * r3[j,2]
      r1M[1] += -A13[M,j,iz,ID] * r3[j,1]
      r1M[2] += -A13[M,j,iz,ID] * r3[j,2]
    end
    r11[1] *= invfac
    r11[2] *= invfac
    r1M[1]*= invfac
    r1M[2]*= invfac

    #@views r2 = invfac * (r2 - A23[:,:,iz,ID] * r3)
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

    s[2,1] = s[2,1] - FacGrav * r11[1]
    s[3,1] = s[3,1] - FacGrav * r1M[1]
    s[2,2] = s[2,2] - FacGrav * r11[2]
    s[3,2] = s[3,2] - FacGrav * r1M[2]
    @unroll for i = 1 : 4
      @atomic :monotonic SchurBand[i + 2,sh + 2,ID] += s[i,1]
      @atomic :monotonic SchurBand[i + 1,sh + 3,ID] += s[i,2]
    end
    if iz > 1
      @atomic :monotonic SchurBand[7,sh-1,ID] -= FacGrav * invfac * B1m_34[1,iz,ID]
      @atomic :monotonic SchurBand[6,sh,ID] -= FacGrav * invfac * B1m_34[2,iz,ID]
    end
    if iz < nz
      @atomic :monotonic SchurBand[2,sh+5,ID] -= FacGrav * invfac * B1p_12[1,iz,ID]
      @atomic :monotonic SchurBand[1,sh+6,ID] -= FacGrav * invfac * B1p_12[2,iz,ID]
    end
  end
end      


@kernel inbounds = true function FillJacDGVertKernel!(A13,A23,A32,B1m_34,B1_1,B1_23,B1_4,
  B2_23,B3_14,B1p_12,C23_2,C14_3,SA,SchurBand,@Const(U),@Const(dz),@Const(DW),@Const(w),fac,
  FacGrav,cS,Phys, ::Val{M}) where {M}

  iz,ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  DoF = @uniform @ndrange()[2]
  Th = @private eltype(U) (M,)
  dpdRhoTh = @private eltype(U) (M,)
  DWS = @localmem eltype(U) (M,M)
  SAL = @localmem eltype(U) (M,M)
  @uniform RhoPos = 1
  @uniform ThPos = 5
  @uniform invcS = eltype(U)(1) / cS
  @uniform wB = w[1]

  if iz == 1
    @. DWS = DW   
  end
  @synchronize 

  if ID <= DoF
    sh = (iz - 1) * 4 + 1
    inv2dz = eltype(U)(2) / dz[iz,ID]
    invdz = eltype(U)(1) / dz[iz,ID]
    @unroll for i = 1 : M
      Th[i] = U[i,iz,ID,ThPos] / U[i,iz,ID,RhoPos]
      dpdRhoTh[i] = eltype(U)(1) / (eltype(U)(1) - Phys.kappa) * Phys.Rd *
          (Phys.Rd * U[i,iz,ID,ThPos] / Phys.p0)^(Phys.kappa / (eltype(U)(1) - Phys.kappa))
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

    @unroll for i = 1 : M - 2
      @unroll for j = 1 : M - 2
        SAL[i,j] = eltype(U)(0)
        @unroll for k = 1 : M - 2
          SAL[i,j] += A32[i,k,iz,ID] * A23[k,j,iz,ID]
        end
        if i == j
          SA[i,j,iz,ID] = fac - (FacGrav / fac) * A13[i+1,j,iz,ID] -
            (eltype(U)(1) / fac) * SAL[i,j]
        else
          SA[i,j,iz,ID] = -(FacGrav / fac) * A13[i+1,j,iz,ID] -
            (eltype(U)(1) / fac) * SAL[i,j]
        end
      end
    end
    @views LUFull!(SA[:,:,iz,ID])

    if iz > 1
      B1m_34[1,iz,ID] = -eltype(U)(1) / (wB * dz[iz,ID])

      dpdRhoThM = eltype(U)(1) / (eltype(U)(1) - Phys.kappa) * Phys.Rd *
        (Phys.Rd * U[M,iz-1,ID,ThPos] ./ Phys.p0).^(Phys.kappa / (eltype(U)(1) - Phys.kappa))
      B1m_34[2,iz,ID] = -dpdRhoThM * invcS / (wB * dz[iz,ID])
    end
    @unroll for i = 1 : M
      B1_23[i,1,iz,ID] = 2.0 * DW[i,1] / dz[iz,ID]
      B1_23[i,2,iz,ID] = 2.0 * DW[i,M] / dz[iz,ID]
    end  

    if iz == 1
      B1_1[iz,ID] = eltype(U)(0)
    else
      B1_23[1,1,iz,ID] = eltype(U)(0)
      B1_1[iz,ID] = dpdRhoTh[1] * invcS * invdz / wB 
    end
    if iz == nz
      B1_4[iz,ID] = eltype(U)(0)
    else
      B1_23[M,2,iz,ID] = eltype(U)(0)
      B1_4[iz,ID] = dpdRhoTh[M] * invcS * invdz / wB 
    end
    if iz < nz
      dpdRhoThP = eltype(U)(1) / (eltype(U)(1) - Phys.kappa) * Phys.Rd *
        (Phys.Rd * U[1,iz+1,ID,ThPos] / Phys.p0)^(Phys.kappa / (eltype(U)(1) - Phys.kappa))
      B1p_12[1,iz,ID] = -dpdRhoThP * invcS / (wB * dz[iz,ID])
      B1p_12[2,iz,ID] = eltype(U)(1) / (wB * dz[iz,ID])
    end

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

    if iz > 1
      ThM = U[M,iz-1,1,ThPos] / U[M,iz-1,1,RhoPos]
      Th1 = U[1,iz-1,1,ThPos] / U[1,iz-1,1,RhoPos]
      dpdRhoThM = eltype(U)(1) / (eltype(U)(1) - Phys.kappa) * Phys.Rd *
        (Phys.Rd * U[M,iz-1,ID,ThPos] / Phys.p0)^(Phys.kappa / (eltype(U)(1) - Phys.kappa))
      # (sh-1,sh-1)  
      @atomic :monotonic SchurBand[4,sh-1,ID] += ThM * dpdRhoThM * invcS / dz[iz-1,ID] / wB
      # (sh,sh-1)
      @atomic :monotonic SchurBand[5,sh-1,ID] += -Th[1] * dpdRhoThM * invcS * invdz / wB
      # (sh-1,sh)
      @atomic :monotonic SchurBand[3,sh,ID] += -ThM * dpdRhoTh[1] * invcS / dz[iz-1,ID] / wB
      # (sh,sh)
      @atomic :monotonic SchurBand[4,sh,ID] += Th[1] * dpdRhoTh[1] * invcS * invdz / wB
      # (sh-2,sh-2)
      @atomic :monotonic SchurBand[4,sh-2,ID] += cS / dz[iz-1,ID] / wB
      # (sh+1,sh-2)
      @atomic :monotonic SchurBand[7,sh-2,ID] += -cS * invdz / wB
      # (sh-2,sh+1)
      @atomic :monotonic SchurBand[1,sh+1,ID] += -cS / dz[iz-1,ID] / wB
      # (sh+1,sh+1)
      @atomic :monotonic SchurBand[4,sh+1,ID] += cS * invdz / wB
      # (sh,sh-2)
      @atomic :monotonic SchurBand[6,sh-2,ID] += -ThM * invdz / wB 
      # (sh+1,sh-1)      
      @atomic :monotonic SchurBand[6,sh-1,ID] += -dpdRhoThM * invdz / wB 
      # (sh-1,sh+1)
      @atomic :monotonic SchurBand[2,sh+1,ID] += Th[1] / wB / dz[iz-1,ID]
      # (sh-2,sh)
      @atomic :monotonic SchurBand[2,sh,ID] += dpdRhoTh[1] / wB / dz[iz-1,ID]
    end


    if iz == 1
      @atomic :monotonic SchurBand[3,sh+1,ID] += inv2dz * DWS[1,1] * Th[1]
      @atomic :monotonic SchurBand[5,sh,ID] += -inv2dz * DWS[1,1] * dpdRhoTh[1] 
      @atomic :monotonic SchurBand[4,sh+1,ID] += inv2dz * cS / wB
    end
    @atomic :monotonic SchurBand[2,sh+2,ID] += inv2dz * DWS[1,M] * Th[M]
    @atomic :monotonic SchurBand[2,sh+3,ID] += inv2dz * DWS[1,M] * dpdRhoTh[M]
    @atomic :monotonic SchurBand[6,sh+1,ID] += inv2dz * DWS[M,1] * Th[1]
    @atomic :monotonic SchurBand[6,sh,ID] += inv2dz * DWS[M,1] * dpdRhoTh[1]
    if iz == nz
      @atomic :monotonic SchurBand[5,sh+2,ID] += inv2dz * DWS[M,M] * Th[M]
      @atomic :monotonic SchurBand[3,sh+3,ID] += -inv2dz * DWS[M,M] * dpdRhoTh[M]
      @atomic :monotonic SchurBand[4,sh+2,ID] += inv2dz * cS / wB
    end
  end  
end

function FillJacDGVert!(JacVert,U,DG,dz,fac,Phys,Param)

  backend = get_backend(U)
  FTB = eltype(U)

  M = JacVert.M
  nz = JacVert.nz
  DoF  = DG.NumI

  JacVert.fac = fac
  JacVert.FacGrav = Phys.Grav

  DWZ = DG.DWZ

  DoFG = 10
  group = (nz, DoFG)
  ndrange = (nz, DoF) 
  @. JacVert.SchurBand = 0.0
  @views @. JacVert.SchurBand[4,:,:] = fac
  KFillJacDGVertKernel! = FillJacDGVertKernel!(backend,group)
  KFillJacDGVertKernel!(JacVert.A13,JacVert.A23,JacVert.A32,JacVert.B1m_34,JacVert.B1_1,
  JacVert.B1_23,JacVert.B1_4, JacVert.B2_23,JacVert.B3_14,JacVert.B1p_12,JacVert.C23_2,
  JacVert.C14_3,JacVert.SA,JacVert.SchurBand,U,dz,DWZ,DG.wZ,fac,Phys.Grav,
  Param.cS,Phys,Val(M);ndrange=ndrange) 

end  

function SchurBoundary!(JacVert)

  backend = get_backend(JacVert.SchurBand)
  FTB = eltype(JacVert.SchurBand)

  M = JacVert.M
  nz = JacVert.nz
  DoF = size(JacVert.SchurBand,3)

  DoFG = 10
  group = (nz, DoFG)
  ndrange = (nz, DoF)
  KSchurBoundaryKernel! = SchurBoundaryKernel!(backend,group)
  KSchurBoundaryKernel!(JacVert.A13,JacVert.A23,JacVert.A32,JacVert.B1m_34,JacVert.B1_1,
    JacVert.B1_23,JacVert.B1_4, JacVert.B2_23,JacVert.B3_14,JacVert.B1p_12,JacVert.C23_2,
    JacVert.C14_3,JacVert.SA,JacVert.SchurBand,JacVert.fac,JacVert.FacGrav,Val(M);ndrange=ndrange)
  group = (DoFG)
  ndrange = (DoF)
  KluBandKernel! = luBandKernel!(backend,group)
  KluBandKernel!(JacVert.SchurBand,ndrange=ndrange)

end  

function ldivVertical!(JacVert,b)

  backend = get_backend(JacVert.SchurBand)
  FTB = eltype(JacVert.SchurBand)

  M = JacVert.M
  nz = JacVert.nz
  DoF = size(JacVert.SchurBand,3)

  invfac = 1.0 / JacVert.fac
  FacGrav = JacVert.FacGrav

  DoFG = 10
  group = (nz, DoFG)
  ndrange = (nz, DoF)
  KldivVerticalFKernel! = ldivVerticalFKernel!(backend,group)
  KldivVerticalFKernel!(JacVert.A13,JacVert.A23,JacVert.A32,JacVert.B1m_34,JacVert.B1_1,
    JacVert.B1_23,JacVert.B1_4, JacVert.B2_23,JacVert.B3_14,JacVert.B1p_12,JacVert.C23_2,
    JacVert.C14_3,JacVert.SA,JacVert.SchurBand,b,JacVert.rs,invfac,FacGrav,Val(M);ndrange=ndrange)

  group = (DoFG)
  ndrange = (DoF)
  KldivVerticalSKernel! = ldivVerticalSKernel!(backend,group)
  KldivVerticalSKernel!(JacVert.SchurBand,JacVert.rs;ndrange=ndrange)

  group = (nz, DoFG)
  ndrange = (nz, DoF)
  KldivVerticalBKernel! = ldivVerticalBKernel!(backend,group)
  KldivVerticalBKernel!(JacVert.A13,JacVert.A23,JacVert.A32,JacVert.B1m_34,JacVert.B1_1,
    JacVert.B1_23,JacVert.B1_4, JacVert.B2_23,JacVert.B3_14,JacVert.B1p_12,
    JacVert.SA,JacVert.SchurBand,b,JacVert.rs,invfac,FacGrav,Val(M);ndrange=ndrange)

end

function Permutation(M,nz)
#Permutation
  N = M * nz
  p = zeros(Int,3*N)
  ii = 0
  @inbounds for iz = 1 : nz
    @inbounds for iv = 1 : 1
      @inbounds for k = 1 : M 
        ii += 1
        p[ii] = k + (iz - 1) * M + (iv - 1) * N
      end
    end
    @inbounds for iv = [3 2] 
      @inbounds for k = 2 : M - 1
        ii += 1
        p[ii] = k + (iz - 1) * M + (iv - 1) * N
      end
    end
  end
  ivw = 2
  ivTh = 3
  @inbounds for iz = 1 : nz
    ii += 1
    p[ii] = 1 + (iz-1) * M  + (ivTh - 1) * N
    ii += 1
    p[ii] = 1 + (iz-1) * M  + (ivw - 1) * N
    ii += 1
    p[ii] = M + (iz - 1) * M + (ivw - 1) * N
    ii += 1
    p[ii] = M + (iz - 1) * M + (ivTh - 1) * N
  end
  return p
end  

function DScalarDMomAc(NZ,DG,cS)
  
  fac = 0.5
  M = DG.OrdPolyZ + 1
  N = NZ * M
  RowInd = Int[]
  ColInd = Int[]
  Val = Float64[]
  D = DG.DWZ
  for iZ = 1 : NZ
    for i = 1 : M
      for j = 1 : M
        push!(RowInd,i+(iZ-1)*M)
        push!(ColInd,j+(iZ-1)*M)
        push!(Val,-D[i,j])
      end
    end  
    if iZ < NZ
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M)  
      push!(Val,-fac/DG.wZ[1])
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M+1)  
      push!(Val,-fac/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M)  
      push!(Val,fac/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M+1)  
      push!(Val,fac/DG.wZ[1])
    end  
  end
  dSdM = sparse(RowInd, ColInd, Val,N,N)

  RowInd = Int[]
  ColInd = Int[]
  Val = Float64[]
  for iZ = 1 : NZ
    if iZ < NZ
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M)  
      push!(Val,-fac/cS/DG.wZ[1])
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M+1)  
      push!(Val,+fac/cS/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M)  
      push!(Val,fac/cS/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M+1)  
      push!(Val,-fac/cS/DG.wZ[1])
    end  
  end
  dSdS = sparse(RowInd, ColInd, Val,N,N)
  return dSdS,dSdM
end

function DMomDScalarAc(NZ,DG,cS)
  
  fac = 0.5
  M = DG.OrdPolyZ + 1
  N = NZ * M
  RowInd = Int[]
  ColInd = Int[]
  Val = Float64[]
  D = DG.DWZ
  for iZ = 1 : NZ
    for i = 1 : M
      for j = 1 : M
        push!(RowInd,i+(iZ-1)*M)
        push!(ColInd,j+(iZ-1)*M)
        push!(Val,-D[i,j])
      end
    end  
    if iZ < NZ
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M)  
      push!(Val,-fac/DG.wZ[1])
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M+1)  
      push!(Val,-fac/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M)  
      push!(Val,fac/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M+1)  
      push!(Val,fac/DG.wZ[1])
    end  
    if iZ == 1
      push!(RowInd,1)  
      push!(ColInd,1)  
      push!(Val,1/DG.wZ[1])
    end    
    if iZ == NZ
      push!(RowInd,M*NZ)  
      push!(ColInd,M*NZ)  
      push!(Val,-1/DG.wZ[1])
    end    
  end
  dMdS = sparse(RowInd, ColInd, Val,N,N)
  RowInd = Int[]
  ColInd = Int[]
  Val = Float64[]
  for iZ = 1 : NZ
    if iZ < NZ
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M)  
      push!(Val,-fac*cS/DG.wZ[1])
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M+1)  
      push!(Val,+fac*cS/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M)  
      push!(Val,fac*cS/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M+1)  
      push!(Val,-fac*cS/DG.wZ[1])
    end  
    if iZ == 1
      push!(RowInd,1)  
      push!(ColInd,1)  
      push!(Val,-cS/DG.wZ[1])
    end    
    if iZ == NZ
      push!(RowInd,M*NZ)  
      push!(ColInd,M*NZ)  
      push!(Val,-cS/DG.wZ[1])
    end    
  end  
  dMdM = sparse(RowInd, ColInd, Val,N,N)
  return dMdS,dMdM
end

function InitJacDG(DG,nz,Param)
  N = (DG.OrdPolyZ + 1) * nz
  dSdS,dSdM = DScalarDMomAc(nz,DG,Param.cS)
  dMdS,dMdM = DMomDScalarAc(nz,DG,Param.cS)
  return dSdS,dSdM,dMdS,dMdM
end  

function JacDG(U,DG,fac,dSdS,dSdM,dMdS,dMdM,z,Phys)
  FTB = eltype(U)
  N = size(dSdM,1)
  RhoPos = 1
  ThPos = 5
  nz = size(U,2)
  M = size(U,1)
  oneM = ones(M)
  NF = size(z,3)
  JacLU = Array{SparseArrays.UMFPACK.UmfpackLU}(undef,size(U,3))
  for ID = 1 : DG.NumI
    @views zCol = z[:,ID]
    diagz = spdiagm(2.0 ./ reshape(vec(oneM*zCol'),N))
    Th = reshape(U[:,:,ID,ThPos]./U[:,:,ID,RhoPos],N)
    dpdRhoTh = reshape( FTB(1) / (FTB(1) - Phys.kappa) * Phys.Rd *
      (Phys.Rd * U[:,:,ID,ThPos] ./ Phys.p0).^(Phys.kappa / (1.0 - Phys.kappa)),N)
    Jac = [sparse(fac*I,N,N) -diagz * dSdM              -diagz* dSdS * diagm(dpdRhoTh)
           sparse((1.0 * Phys.Grav)*I,N,N) sparse(fac*I,N,N) - diagz * dMdM -diagz* dMdS * diagm(dpdRhoTh)
           spzeros(N,N) -diagz * dSdM * diagm(Th)  sparse(fac*I,N,N) - diagz * diagm(Th) * dSdS * diagm(dpdRhoTh)]
    JacLU[ID] = lu(Jac)           
  end
#-----------  
  ID = 1
  @views zCol = z[:,ID]
    diagz = spdiagm(2.0 ./ reshape(vec(oneM*zCol'),N))
    Th = reshape(U[:,:,ID,ThPos]./U[:,:,ID,RhoPos],N)
    dpdRhoTh = reshape( FTB(1) / (FTB(1) - Phys.kappa) * Phys.Rd *
      (Phys.Rd * U[:,:,ID,ThPos] ./ Phys.p0).^(Phys.kappa / (1.0 - Phys.kappa)),N)
    Jac = [sparse(fac*I,N,N) -diagz * dSdM              -diagz* dSdS * diagm(dpdRhoTh)
           sparse((1.0 * Phys.Grav)*I,N,N) sparse(fac*I,N,N) - diagz * dMdM -diagz* dMdS * diagm(dpdRhoTh)
           spzeros(N,N) -diagz * dSdM * diagm(Th)  sparse(fac*I,N,N) - diagz * diagm(Th) * dSdS * diagm(dpdRhoTh)]
#-----------  
  return JacLU,Jac
# return JacLU
end

function JacDGT(U,DG,fac,dSdS,dSdM,dMdS,dMdM,z,Phys)
  FTB = eltype(U)
  N = size(dSdM,1)
  RhoPos = 1
  ThPos = 5
  nz = size(U,2)
  M = size(U,1)
  oneM = ones(M)
  NF = size(z,3)
  JacLU = Array{SparseArrays.UMFPACK.UmfpackLU}(undef,size(U,3))
    ID = 1  
    @views zCol = z[:,ID]
    diagz = spdiagm(2.0 ./ reshape(vec(oneM*zCol'),N))
    Th = reshape(U[:,:,ID,ThPos]./U[:,:,ID,RhoPos],N)
    dpdRhoTh = reshape( FTB(1) / (FTB(1) - Phys.kappa) * Phys.Rd *
      (Phys.Rd * U[:,:,ID,ThPos] ./ Phys.p0).^(Phys.kappa / (1.0 - Phys.kappa)),N)
    Jac = [sparse(fac*I,N,N)  -diagz* dSdS * diagm(dpdRhoTh) -diagz * dSdM
           spzeros(N,N) sparse(fac*I,N,N) - diagz*diagm(Th)*dSdS*diagm(dpdRhoTh)  -diagz*dSdM*diagm(Th)
           sparse((Phys.Grav)*I,N,N) -diagz*dMdS*diagm(dpdRhoTh) sparse(fac*I,N,N)-diagz*dMdM]
    JacLU[ID] = lu(Jac)
  return JacLU,Jac
end

