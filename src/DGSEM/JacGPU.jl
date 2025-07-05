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

@kernel inbounds = true function ldivVerticalFKernel!(JacVert,b,rs,invfac,FacGrav, ::Val{M}) where {M}

  iz,ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  DoF = @uniform @ndrange()[2]

  @uniform FT = eltype(JacVert.SchurBand)

  r1 = @private FT (M,)
  r2 = @private FT (M-2,)
  r3 = @private FT (M-2,)

  @uniform RhoPos = 1
  @uniform wPos = 4
  @uniform ThPos = 5

  if ID <= DoF
    sh = (iz - 1) * 4
    rs[sh + 1,ID] = b[1,iz,ID,ThPos]
    rs[sh + 2,ID] = b[1,iz,ID,wPos]
    rs[sh + 3,ID] = b[M,iz,ID,wPos]
    rs[sh + 4,ID] = b[M,iz,ID,ThPos]
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
        r3[i] += -invfac * JacVert.A32[i,j,iz,ID] * r2[j]
      end
    end  

    @views ldivFull!(JacVert.SA[:,:,iz,ID],r3)

    r11 = r1[1]
    r1M = r1[M]
    @unroll for j = 1 : M - 2
      r11 += -JacVert.A13[1,j] * r3[j]
      r1M += -JacVert.A13[M,j] * r3[j]
    end
    r11 *= invfac
    r1M *= invfac

    @unroll for i = 1 : M - 2
      @unroll for j = 1 : M - 2
        r2[i] += -JacVert.A23[i,j,iz,ID] * r3[j]
      end
      r2[i] *= invfac
    end
    @unroll for j = 1 : M - 2
      rs[sh + 1,ID] += -JacVert.C14_3[1,j,iz,ID] * r3[j]
      rs[sh + 4,ID] += -JacVert.C14_3[2,j,iz,ID] * r3[j]
      rs[sh + 2,ID] += -JacVert.C23_2[1,j,iz,ID] * r2[j]
      rs[sh + 3,ID] += -JacVert.C23_2[2,j,iz,ID] * r2[j]
    end
    rs[sh + 2,ID] += -FacGrav * r11
    rs[sh + 3,ID] += -FacGrav * r1M
  end
end

@kernel inbounds = true function ldivVerticalSKernel!(A,rs)

  ID, = @index(Global, NTuple)

  DoF = @uniform @ndrange()[1]

  if ID <= DoF
    @views ldivBand!(A[:,:,ID],rs,3,3)
  end  
end  

@kernel inbounds = true function ldivVerticalBKernel!(JacVert,b,rs,invfac,FacGrav,  ::Val{M}) where {M}

  iz,ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  DoF = @uniform @ndrange()[2]

  @uniform FT = eltype(JacVert.SchurBand)

  r1 = @private FT (M,)
  r2 = @private FT (M-2,)
  r3 = @private FT (M-2,)

  @uniform RhoPos = 1
  @uniform wPos = 4
  @uniform ThPos = 5

  if ID <= DoF
    sh = (iz - 1) * 4

    b[1,iz,ID,ThPos] = rs[sh + 1]
    b[1,iz,ID,wPos] = rs[sh + 2]
    b[M,iz,ID,wPos] = rs[sh + 3]
    b[M,iz,ID,ThPos] = rs[sh + 4]
    @unroll for i = 1 : M
      r1[i] = b[i,iz,ID,RhoPos]
    end  
    @unroll for i = 1 : M - 2
      r2[i] = b[i+1,iz,ID,ThPos]
      r3[i] = b[i+1,iz,ID,wPos]
    end  
    r1[1] -= JacVert.B1_1[iz,ID] * rs[sh + 1]
    r1[M] -= JacVert.B1_4[iz,ID] * rs[sh + 4]
    @unroll for i = 1 : M
      r1[i] -= JacVert.B1_23[i,1,iz,ID] * rs[sh + 2] + JacVert.B1_23[i,2,iz,ID] * rs[sh + 3]
    end
    @unroll for i = 1 : M - 2
      r2[i] -= JacVert.B2_23[i,1,iz,ID] * rs[sh + 2] + JacVert.B2_23[i,2,iz,ID] * rs[sh + 3]
      r3[i] -= JacVert.B3_14[i,1,iz,ID] * rs[sh + 1] + JacVert.B3_14[i,2,iz,ID] * rs[sh + 4]
    end
    if iz > 1
      r1[1] -= JacVert.B1m_34[1,iz,ID] * rs[sh-1] + JacVert.B1m_34[2,iz,ID] * rs[sh]
    end
    if iz < nz
      r1[M] -= JacVert.B1p_12[1,iz,ID] * rs[sh+5] + JacVert.B1p_12[2,iz,ID] * rs[sh+6]
    end
    @unroll for i = 1 : M - 2
      r3[i] += -invfac * FacGrav * r1[i+1]
      @unroll for j = 1 : M - 2
        r3[i] += -invfac * JacVert.A32[i,j,iz,ID] * r2[j]
      end
    end
    @views ldivFull!(JacVert.SA[:,:,iz,ID],r3)
    @unroll for i = 1 : M
      @unroll for j = 1 : M - 2
        r1[i] += -JacVert.A13[i,j,iz,ID] * r3[j]
      end
      r1[i] *= invfac
    end
    @unroll for i = 1 : M - 2
      @unroll for j = 1 : M - 2
        r2[i] += -JacVert.A23[i,j,iz,ID] * r3[j]
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

@kernel inbounds = true function SchurBoundaryKernel!(JacVert, ::Val{M}) where {M}

  iz,ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  DoF = @uniform @ndrange()[2]

  @uniform FT = eltype(JacVert.SchurBand)

  r1 = @private FT (M,2,)
  r2 = @private FT (M-2,2,)
  r3 = @private FT (M-2,2,)
  r11 = @private FT (2,)
  r1M = @private FT (2,)
  s = @private FT (4,2,)
  invfac = 1.0 / JacVert.fac

  @uniform FacGrav = JacVert.FacGrav

  if ID <= DoF
    sh = (iz - 1) * 4
#   Column 1 and 4
    @unroll for i = 1 : M - 2
      r3[i,1] = JacVert.B3_14[i,1,iz,ID]
      r3[i,2] = JacVert.B3_14[i,2,iz,ID]
    end  

    @views ldivFull2!(JacVert.SA[:,:,iz,ID],r3)

    r11[1] = -JacVert.B1_1[iz,ID]
    r11[2] = FT(0)
    r1M[1] = FT(0)
    r1M[2] = -JacVert.B1_4[iz,ID]
    @unroll for j = 1 : M - 2
      r11[1] += JacVert.A13[1,j,iz,ID] * r3[j,1]
      r11[2] += JacVert.A13[1,j,iz,ID] * r3[j,2]
      r1M[1] += JacVert.A13[M,j,iz,ID] * r3[j,1]
      r1M[2] += JacVert.A13[M,j,iz,ID] * r3[j,2]
    end  
    r11[1] *= -invfac
    r11[2] *= -invfac
    r1M[1] *= -invfac
    r1M[2] *= -invfac

    @unroll for i = 1 : M - 2
      r2[i,1] = FT(0)
      r2[i,2] = FT(0)
      @unroll for j = 1 : M - 2
        r2[i,1] += JacVert.A23[i,j,iz,ID] * r3[j,1]
        r2[i,2] += JacVert.A23[i,j,iz,ID] * r3[j,2]
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
      s[1,1] += -JacVert.C14_3[1,j,iz,ID] * r3[j,1]
      s[1,2] += -JacVert.C14_3[1,j,iz,ID] * r3[j,2]
      s[4,1] += -JacVert.C14_3[2,j,iz,ID] * r3[j,1]
      s[4,2] += -JacVert.C14_3[2,j,iz,ID] * r3[j,2]
      s[2,1] += -JacVert.C23_2[1,j,iz,ID] * r2[j,1]
      s[2,2] += -JacVert.C23_2[1,j,iz,ID] * r2[j,2]
      s[3,1] += -JacVert.C23_2[2,j,iz,ID] * r2[j,1]
      s[3,2] += -JacVert.C23_2[2,j,iz,ID] * r2[j,2]
    end
    s[2,1] = s[2,1] - FacGrav * r11[1]
    s[3,1] = s[3,1] - FacGrav * r1M[1]
    s[2,2] = s[2,2] - FacGrav * r11[2]
    s[3,2] = s[3,2] - FacGrav * r1M[2]
    @unroll for i = 1 : 4
      @atomic :monotonic JacVert.SchurBand[i + 3,sh + 1,ID] += s[i,1]
      @atomic :monotonic JacVert.SchurBand[i,sh + 4,ID] += s[i,2]
    end

    @unroll for i = 1 : M
      r1[i,1] = JacVert.B1_23[i,1,iz,ID]
      r1[i,2] = JacVert.B1_23[i,2,iz,ID]
    end  
    @unroll for i = 1 : M - 2
      r2[i,1] = JacVert.B2_23[i,1,iz,ID]
      r2[i,2] = JacVert.B2_23[i,2,iz,ID]
    end  
    @unroll for i = 1 : M - 2
      r3[i,1] = FacGrav * r1[i+1,1]
      r3[i,2] = FacGrav * r1[i+1,2]
      @unroll for j = 1 : M - 2
        r3[i,1] += JacVert.A32[i,j,iz,ID] * r2[j,1]
        r3[i,2] += JacVert.A32[i,j,iz,ID] * r2[j,2]
      end
      r3[i,1] *= -invfac
      r3[i,2] *= -invfac
    end

    @views ldivFull2!(JacVert.SA[:,:,iz,ID],r3)

    #r11 = invfac * (r1[1:1,:] - JacVert.A13[1:1,:,iz,ID] * r3[:,:])
    #r1M = invfac * (r1[M:M,:] - JacVert.A13[M:M,:,iz,ID] * r3[:,:])
    r11[1] = r1[1,1]
    r11[2] = r1[1,2]
    r1M[1] = r1[M,1]
    r1M[2] = r1[M,2]
    @unroll for j = 1 : M - 2
      r11[1] += -JacVert.A13[1,j,iz,ID] * r3[j,1]
      r11[2] += -JacVert.A13[1,j,iz,ID] * r3[j,2]
      r1M[1] += -JacVert.A13[M,j,iz,ID] * r3[j,1]
      r1M[2] += -JacVert.A13[M,j,iz,ID] * r3[j,2]
    end
    r11[1] *= invfac
    r11[2] *= invfac
    r1M[1]*= invfac
    r1M[2]*= invfac

    #@views r2 = invfac * (r2 - JacVert.A23[:,:,iz,ID] * r3)
    @unroll for i = 1 : M - 2
      @unroll for j = 1 : M - 2
        r2[i,1] += -JacVert.A23[i,j,iz,ID] * r3[j,1]
        r2[i,2] += -JacVert.A23[i,j,iz,ID] * r3[j,2]
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
      s[1,1] += - JacVert.C14_3[1,j,iz,ID] * r3[j,1]
      s[1,2] += - JacVert.C14_3[1,j,iz,ID] * r3[j,2]
      s[4,1] += - JacVert.C14_3[2,j,iz,ID] * r3[j,1]
      s[4,2] += - JacVert.C14_3[2,j,iz,ID] * r3[j,2]
      s[2,1] += - JacVert.C23_2[1,j,iz,ID] * r2[j,1]
      s[2,2] += - JacVert.C23_2[1,j,iz,ID] * r2[j,2]
      s[3,1] += - JacVert.C23_2[2,j,iz,ID] * r2[j,1]
      s[3,2] += - JacVert.C23_2[2,j,iz,ID] * r2[j,2]
    end

    s[2,1] = s[2,1] - FacGrav * r11[1]
    s[3,1] = s[3,1] - FacGrav * r1M[1]
    s[2,2] = s[2,2] - FacGrav * r11[2]
    s[3,2] = s[3,2] - FacGrav * r1M[2]
    @unroll for i = 1 : 4
      @atomic :monotonic JacVert.SchurBand[i + 2,sh + 2,ID] += s[i,1]
      @atomic :monotonic JacVert.SchurBand[i + 1,sh + 3,ID] += s[i,2]
    end
    if iz > 1
      @atomic :monotonic JacVert.SchurBand[7,sh-1,ID] -= FacGrav * invfac * JacVert.B1m_34[1,iz,ID]
      @atomic :monotonic JacVert.SchurBand[6,sh,ID] -= FacGrav * invfac * JacVert.B1m_34[2,iz,ID]
    end
    if iz < nz
      @atomic :monotonic JacVert.SchurBand[2,sh+5,ID] -= FacGrav * invfac * JacVert.B1p_12[1,iz,ID]
      @atomic :monotonic JacVert.SchurBand[1,sh+6,ID] -= FacGrav * invfac * JacVert.B1p_12[2,iz,ID]
    end
  end
end      


@kernel inbounds = true function FillJacDGVertKernel!(JacVert,@Const(U),@Const(dz),@Const(DW),wB,fac,
  cS,Phys, ::Val{M}) where {M}

  iz,ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  DoF = @uniform @ndrange()[2]
  Th = @private eltype(U) (M,)
  dpdRhoTh = @private eltype(U) (M,)
  DWS = @localmem eltype(U) (M,M)
  @uniform RhoPos = 1
  @uniform ThPos = 5
  @uniform invcS = eltype(U)(1) / cS

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
        JacVert.A13[i,j,iz,ID] = inv2dz * DWS[i,j+1]
      end
    end  
    @unroll for i = 1 : M - 2
      @unroll for j = 1 : M - 2
        JacVert.A23[i,j,iz,ID] = inv2dz * DWS[i+1,j+1] * Th[j+1]  
        JacVert.A32[i,j,iz,ID] = inv2dz * DWS[i+1,j+1] * dpdRhoTh[j+1]  
      end
    end  
    @unroll for i = 1 : M - 2
      @unroll for j = 1 : M - 2
        JacVert.SA[i,j,iz,ID] = eltype(U)(0)
        @unroll for k = 1 : M - 2
          JacVert.SA[i,j,iz,ID] += JacVert.A32[i,k,iz,ID] * JacVert.A23[k,j,iz,ID]
        end
        if i == j
          JacVert.SA[i,j,iz,ID] = fac - (JacVert.FacGrav / fac) * JacVert.A13[i+1,j,iz,ID] -
            (eltype(U)(1) / fac) * JacVert.SA[i,j,iz,ID]
        else
          JacVert.SA[i,j,iz,ID] = -(JacVert.FacGrav / fac) * JacVert.A13[i+1,j,iz,ID] -
            (eltype(U)(1) / fac) * JacVert.SA[i,j,iz,ID]
        end
      end
    end
    @views LUFull!(JacVert.SA[:,:,iz,ID])

    if iz > 1
      JacVert.B1m_34[1,iz,ID] = -eltype(U)(1) / (wB * dz[iz-1,ID])

      dpdRhoThM = eltype(U)(1) / (eltype(U)(1) - Phys.kappa) * Phys.Rd *
        (Phys.Rd * U[M,iz-1,ID,ThPos] ./ Phys.p0).^(Phys.kappa / (eltype(U)(1) - Phys.kappa))
      JacVert.B1m_34[2,iz,ID] = -dpdRhoThM * invcS / (wB * dz[iz-1,ID])
    end
    if iz == 1
      JacVert.B1_1[iz,ID] = eltype(U)(0)
    else
      JacVert.B1_23[1,1,iz,ID] = eltype(U)(0)
      JacVert.B1_1[iz,ID] = dpdRhoTh[1] * invcS * invdz / wB 
    end
    if iz == nz
      JacVert.B1_4[iz,ID] = eltype(U)(0)
    else
      JacVert.B1_23[M,2,iz,ID] = eltype(U)(0)
      JacVert.B1_4[iz,ID] = dpdRhoTh[M] * invcS * invdz / wB 
    end
    if iz < nz
      dpdRhoThP = eltype(U)(1) / (eltype(U)(1) - Phys.kappa) * Phys.Rd *
        (Phys.Rd * U[1,iz+1,ID,ThPos] / Phys.p0)^(Phys.kappa / (eltype(U)(1) - Phys.kappa))
      JacVert.B1p_12[1,iz,ID] = -dpdRhoThP * invcS / (wB * dz[iz+1,ID])
      JacVert.B1p_12[2,iz,ID] = eltype(U)(1) / (wB * dz[iz+1,ID])
    end

    @unroll for i = 1 : M - 2
      JacVert.B2_23[i,1,iz,ID] = inv2dz * DWS[i+1,1] * Th[1]
      JacVert.B2_23[i,2,iz,ID] = inv2dz * DWS[i+1,M] * Th[M]
      JacVert.B3_14[i,1,iz,ID] = inv2dz * DWS[i+1,1] * dpdRhoTh[1]
      JacVert.B3_14[i,2,iz,ID] = inv2dz * DWS[i+1,M] * dpdRhoTh[M]
      JacVert.C23_2[1,i,iz,ID] = inv2dz * DWS[1,i+1] * dpdRhoTh[i+1]
      JacVert.C23_2[2,i,iz,ID] = inv2dz * DWS[M,i+1] * dpdRhoTh[i+1]
      JacVert.C14_3[1,i,iz,ID] = inv2dz * DWS[1,i+1] * Th[i+1]
      JacVert.C14_3[2,i,iz,ID] = inv2dz * DWS[M,i+1] * Th[i+1]
    end

    if iz > 1
      ThM = U[M,iz-1,1,ThPos] / U[M,iz-1,1,RhoPos]
      Th1 = U[1,iz-1,1,ThPos] / U[1,iz-1,1,RhoPos]
      dpdRhoThM = eltype(U)(1) / (eltype(U)(1) - Phys.kappa) * Phys.Rd *
        (Phys.Rd * U[M,iz-1,ID,ThPos] / Phys.p0)^(Phys.kappa / (eltype(U)(1) - Phys.kappa))
      @atomic :monotonic JacVert.SchurBand[4,sh-1,ID] += ThM * dpdRhoThM * invcS / dz[iz-1,ID] / wB
      @atomic :monotonic JacVert.SchurBand[5,sh-1,ID] += -Th[1] * dpdRhoThM * invcS / dz[iz-1,ID] / wB
      @atomic :monotonic JacVert.SchurBand[3,sh,ID] += -ThM * dpdRhoTh[1] * invcS * invdz / wB
      @atomic :monotonic JacVert.SchurBand[4,sh,ID] += Th[1] * dpdRhoTh[1] * invcS * invdz  / wB
      @atomic :monotonic JacVert.SchurBand[4,sh-2,ID] += cS / dz[iz-1,ID] / wB
      @atomic :monotonic JacVert.SchurBand[7,sh-2,ID] += -cS / dz[iz-1,ID] / wB
      @atomic :monotonic JacVert.SchurBand[1,sh+1,ID] += -cS * invdz / wB
      @atomic :monotonic JacVert.SchurBand[4,sh+1,ID] += cS * invdz / wB
      @atomic :monotonic JacVert.SchurBand[6,sh-2,ID] += -ThM * invdz / wB 
      @atomic :monotonic JacVert.SchurBand[6,sh-1,ID] += -dpdRhoThM * invdz / wB 
      @atomic :monotonic JacVert.SchurBand[2,sh+1,ID] += Th[1] / wB / dz[iz-1,ID]
      @atomic :monotonic JacVert.SchurBand[2,sh,ID] += dpdRhoTh[1] / wB / dz[iz-1,ID]
    end


    if iz == 1
      @atomic :monotonic JacVert.SchurBand[3,sh+1,ID] += inv2dz * DWS[1,1] * Th[1]
      @atomic :monotonic JacVert.SchurBand[5,sh,ID] += -inv2dz * DWS[1,1] * dpdRhoTh[1] 
      @atomic :monotonic JacVert.SchurBand[4,sh+1,ID] += inv2dz * cS / wB
    end
    @atomic :monotonic JacVert.SchurBand[2,sh+2,ID] += inv2dz * DWS[1,M] * Th[M]
    @atomic :monotonic JacVert.SchurBand[2,sh+3,ID] += inv2dz * DWS[1,M] * dpdRhoTh[M]
    @atomic :monotonic JacVert.SchurBand[6,sh+1,ID] += inv2dz * DWS[M,1] * Th[1]
    @atomic :monotonic JacVert.SchurBand[6,sh,ID] += inv2dz * DWS[M,1] * dpdRhoTh[1]
    if iz == nz
      @atomic :monotonic JacVert.SchurBand[5,sh+2,ID] += inv2dz * DWS[M,M] * Th[M]
      @atomic :monotonic JacVert.SchurBand[3,sh+3,ID] += -inv2dz * DWS[M,M] * dpdRhoTh[M]
      @atomic :monotonic JacVert.SchurBand[4,sh+2,ID] += inv2dz * cS / wB
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

  DoFG = 10
  group = (nz, DoFG)
  ndrange = (nz, DoF) 
  @. JacVert.SchurBand = 0.0
  @views @. JacVert.SchurBand[4,:,:] = fac
  KFillJacDGVertKernel! = FillJacDGVertKernel!(backend,group)
  KFillJacDGVertKernel!(JacVert,U,dz,DG.DWZ,DG.wZ[1],fac,Param.cS,Phys,Val(M);ndrange=ndrange) 

end  

function FillJacDGVertOld!(JacVert,U,DG,dz,fac,Phys,Param)
  backend = get_backend(U)
  FTB = eltype(U)
  M = JacVert.M
  M2 = JacVert.M - 2
  M1 = M - 1
  M23= M2 * 2 + M
  nz = JacVert.nz
  JacVert.fac = fac
  JacVert.FacGrav = Phys.Grav
  DW = DG.DWZ
# DW = LowOrder(DG)
  wZ = DG.wZ
  ThPos = 5
  RhoPos = 1
  Th = zeros(M)
  dpdRhoTh = zeros(M)
  cS = Param.cS
  invcS = 1.0 / cS
  DoF  = DG.NumI
  S = zeros(2,2)
  @inbounds for ID = 1 : DoF
    sh = 1
    @views @. JacVert.SchurBand[:,:,ID] = 0.0
    @views @. JacVert.SchurBand[4,:,ID] = fac
    @inbounds for iz = 1 : nz
      @views @. Th = U[:,iz,ID,ThPos]/U[:,iz,ID,RhoPos]
      @views @. dpdRhoTh = FTB(1) / (FTB(1) - Phys.kappa) * Phys.Rd *
        (Phys.Rd * U[:,iz,ID,ThPos] / Phys.p0)^(Phys.kappa / (1.0 - Phys.kappa))  
      @views @. JacVert.A13[:,:,iz,ID] = (2 / dz[iz,ID]) * DW[:,2:M1]
      @inbounds for j = 2 : M1
        @views @. JacVert.A23[:,j-1,iz,ID] = 2 / dz[iz,ID] * DW[2:M1,j] * Th[j]  
        @views @. JacVert.A32[:,j-1,iz,ID] = 2 / dz[iz,ID] * DW[2:M1,j] * dpdRhoTh[j]  
      end  
      @inbounds for i = 1 : M2
        @inbounds for j = 1 : M2
          JacVert.SA[i,j,iz,ID] = 0.0
          @inbounds for k = 1 : M2
            JacVert.SA[i,j,iz,ID] += JacVert.A32[i,k,iz,ID] * JacVert.A23[k,j,iz,ID]
          end
          if i == j
            JacVert.SA[i,j,iz,ID] = fac - (JacVert.FacGrav / fac) * JacVert.A13[i+1,j,iz,ID] -
              (1.0 / fac) * JacVert.SA[i,j,iz,ID]
          else
            JacVert.SA[i,j,iz,ID] = -(JacVert.FacGrav / fac) * JacVert.A13[i+1,j,iz,ID] -
              (1.0 / fac) * JacVert.SA[i,j,iz,ID]
          end    
        end
      end  
      @views LUFull!(JacVert.SA[:,:,iz,ID])  

      if iz > 1
        JacVert.B1m_34[1,iz,ID] = -1.0 / (wZ[1] * dz[iz-1,ID])

        dpdRhoThM = FTB(1) / (FTB(1) - Phys.kappa) * Phys.Rd *
          (Phys.Rd * U[M,iz-1,ID,ThPos] ./ Phys.p0).^(Phys.kappa / (1.0 - Phys.kappa))  
        JacVert.B1m_34[2,iz,ID] = -dpdRhoThM * invcS / (wZ[1] * dz[iz-1,ID])
      end  
      if iz == 1
        JacVert.B1_1[iz,ID] = 0
      else
        JacVert.B1_23[1,1,iz,ID] = 0.0
        JacVert.B1_1[iz,ID] = dpdRhoTh[1] * invcS / (wZ[1] * dz[iz,ID])
      end  
      if iz == nz
        JacVert.B1_4[iz,ID] = 0
      else
        JacVert.B1_23[M,2,iz,ID] = 0.0
        JacVert.B1_4[iz,ID] = dpdRhoTh[M] * invcS / (wZ[1] * dz[iz,ID])
      end  
      if iz < nz
        dpdRhoThP = FTB(1) / (FTB(1) - Phys.kappa) * Phys.Rd *
          (Phys.Rd * U[1,iz+1,ID,ThPos] / Phys.p0)^(Phys.kappa / (1.0 - Phys.kappa))  
        JacVert.B1p_12[1,iz,ID] = -dpdRhoThP * invcS / (wZ[1] * dz[iz+1,ID])
        JacVert.B1p_12[2,iz,ID] = 1.0 / (wZ[1] * dz[iz+1,ID])
      end  

      @views @. JacVert.B2_23[:,1,iz,ID] = 2.0 * DW[2:M1,1] * Th[1] / dz[iz,ID]
      @views @. JacVert.B2_23[:,2,iz,ID] = 2.0 * DW[2:M1,M] * Th[M] / dz[iz,ID]
      @views @. JacVert.B3_14[:,1,iz,ID] = 2.0 * DW[2:M1,1] * dpdRhoTh[1] / dz[iz,ID]
      @views @. JacVert.B3_14[:,2,iz,ID] = 2.0 * DW[2:M1,M] * dpdRhoTh[M] / dz[iz,ID]

      @views @. JacVert.C23_2[1,:,iz,ID] = 2.0 * DW[1,2:M1] * dpdRhoTh[2:M1] / dz[iz,ID]
      @views @. JacVert.C23_2[2,:,iz,ID] = 2.0 * DW[M,2:M1] * dpdRhoTh[2:M1] / dz[iz,ID]
      @views @. JacVert.C14_3[1,:,iz,ID] = 2.0 * DW[1,2:M1] * Th[2:M1] / dz[iz,ID]
      @views @. JacVert.C14_3[2,:,iz,ID] = 2.0 * DW[M,2:M1] * Th[2:M1] / dz[iz,ID]

      if iz > 1 
        ThM = U[M,iz-1,1,ThPos] / U[M,iz-1,1,RhoPos]  
        Th1 = U[1,iz-1,1,ThPos] / U[1,iz-1,1,RhoPos]  
        dpdRhoThM = FTB(1) / (FTB(1) - Phys.kappa) * Phys.Rd *
          (Phys.Rd * U[M,iz-1,ID,ThPos] / Phys.p0)^(Phys.kappa / (1.0 - Phys.kappa))  
        S[1,1] = ThM * dpdRhoThM * invcS / dz[iz-1,ID] / wZ[1]
        S[2,1] = -Th[1] * dpdRhoThM * invcS / dz[iz-1,ID] / wZ[1]
        S[1,2] = -ThM * dpdRhoTh[1] * invcS / dz[iz,ID] / wZ[1]
        S[2,2] = Th[1] * dpdRhoTh[1] * invcS / dz[iz,ID] / wZ[1]
        JacVert.SchurBand[4,sh-1,ID] += S[1,1]
        JacVert.SchurBand[5,sh-1,ID] += S[2,1]
        JacVert.SchurBand[3,sh,ID] += S[1,2]
        JacVert.SchurBand[4,sh,ID] += S[2,2]
        S[1,1] = cS / dz[iz-1,ID] / wZ[1]
        S[2,1] = -cS / dz[iz-1,ID] / wZ[1]
        S[1,2] = -cS / dz[iz,ID] / wZ[1]
        S[2,2] = cS / dz[iz,ID] / wZ[1]
        JacVert.SchurBand[4,sh-2,ID] += S[1,1]
        JacVert.SchurBand[1,sh+1,ID] += S[1,2]
        JacVert.SchurBand[7,sh-2,ID] += S[2,1]
        JacVert.SchurBand[4,sh+1,ID] += S[2,2]
        JacVert.SchurBand[6,sh-2,ID] += -ThM / wZ[1] / dz[iz,ID]
        JacVert.SchurBand[6,sh-1,ID] += -dpdRhoThM / wZ[1] / dz[iz,ID]
        JacVert.SchurBand[2,sh+1,ID] += Th[1] / wZ[1] / dz[iz-1,ID]
        JacVert.SchurBand[2,sh,ID] += dpdRhoTh[1] / wZ[1] / dz[iz-1,ID]
      end  
      if iz == 1
        JacVert.SchurBand[3,sh+1,ID] += 2.0 * DW[1,1] * Th[1] / dz[iz,ID]
        JacVert.SchurBand[5,sh,ID] += -2.0 * DW[1,1] * dpdRhoTh[1] / dz[iz,ID]
        JacVert.SchurBand[4,sh+1,ID] += 2.0 * cS / dz[iz,ID] / wZ[1]
      end  
      JacVert.SchurBand[2,sh+2,ID] += 2.0 * DW[1,M] * Th[M] / dz[iz,ID]
      JacVert.SchurBand[2,sh+3,ID] += 2.0 * DW[1,M] * dpdRhoTh[M] / dz[iz,ID]
      JacVert.SchurBand[6,sh+1,ID] += 2.0 * DW[M,1] * Th[1] / dz[iz,ID]
      JacVert.SchurBand[6,sh,ID] += 2.0 * DW[M,1] * dpdRhoTh[1] / dz[iz,ID]
      if iz == nz
        JacVert.SchurBand[5,sh+2,ID] += 2.0 * DW[M,M] * Th[M] / dz[iz,ID]
        JacVert.SchurBand[3,sh+3,ID] += -2.0 * DW[M,M] * dpdRhoTh[M] / dz[iz,ID]
        JacVert.SchurBand[4,sh+2,ID] += 2.0 * cS / dz[iz,ID] / wZ[1]
      end  
      sh += 4
    end
  end
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
  KSchurBoundaryKernel!(JacVert,Val(M);ndrange=ndrange)

  group = (DoFG)
  ndrange = (DoF)
  KluBandKernel! = luBandKernel!(backend,group)
  KluBandKernel!(JacVert.SchurBand,ndrange=ndrange)

end  

function SchurBoundaryOld!(JacVert)
  M = JacVert.M
  nz = JacVert.nz
  M2 = M - 2
  invfac = 1.0 / JacVert.fac
  FacGrav = JacVert.FacGrav
  r1 = zeros(M,2)
  r2 = zeros(M2,2)
  r3 = zeros(M2,2)
  s = zeros(4,2)
  r11 = zeros(2)
  r1M = zeros(2)
  DoF = size(JacVert.A13,4)
  @inbounds for ID = 1 : DoF
    @views A22B = JacVert.SchurBand[:,:,ID]
    @inbounds for iz = 1 : nz
      sh = (iz - 1) * 4  
#     Column 1 and 4
      @views @. r3 = JacVert.B3_14[:,:,iz,ID]   

      @views ldivFull2!(JacVert.SA[:,:,iz,ID],r3)

      #r11 = -invfac * (JacVert.A13[1:1,:,iz,ID] * r3[:,:])
      #r1M = -invfac * (JacVert.A13[M:M,:,iz,ID] * r3[:,:])
      #r11[1,1] += invfac * JacVert.B1_1[iz,ID]
      #r1M[1,2] += invfac * JacVert.B1_4[iz,ID]

      r11[1] = -JacVert.B1_1[iz,ID]
      r11[2] = 0.0
      r1M[1] = 0.0
      r1M[2] = -JacVert.B1_4[iz,ID]
      @inbounds for j = 1 : M2
        @views @. r11 += JacVert.A13[1,j,iz,ID] * r3[j,:] 
        @views @. r1M += JacVert.A13[M,j,iz,ID] * r3[j,:] 
      end    
      r11 .*= -invfac
      r1M .*= -invfac

      #@views r2 = -invfac * (JacVert.A23[:,:,iz,ID] * r3)
      @inbounds for i = 1 : M2
        @views @. r2[i,:] = 0.0 
        @inbounds for j = 1 : M2
          @views @. r2[i,:] += JacVert.A23[i,j,iz,ID] * r3[j,:]
        end  
        @views @. r2[i,:] *= -invfac
      end  

      #@views s[[1,4],:] = -JacVert.C14_3[:,:,iz,ID] * r3
      #@views s[[2,3],:] = -JacVert.C23_2[:,:,iz,ID] * r2
      @. s = 0.0
      @inbounds for j = 1 : M2
        @views @. s[1,:] += -JacVert.C14_3[1,j,iz,ID] * r3[j,:]
        @views @. s[4,:] += -JacVert.C14_3[2,j,iz,ID] * r3[j,:]
        @views @. s[2,:] += -JacVert.C23_2[1,j,iz,ID] * r2[j,:]
        @views @. s[3,:] += -JacVert.C23_2[2,j,iz,ID] * r2[j,:]
      end
      s[2,1] = s[2,1] - FacGrav * r11[1]
      s[3,1] = s[3,1] - FacGrav * r1M[1]
      s[2,2] = s[2,2] - FacGrav * r11[2]
      s[3,2] = s[3,2] - FacGrav * r1M[2]
      #@views A22B[sh + 1:sh + 4,[sh + 1,sh + 4]] .+= s
      @inbounds for i = 1 : 4
         #A22B[i + sh - (sh +1) + 4,sh + 1] += s[i,1]
         A22B[i + 3,sh + 1] += s[i,1]
         #A22B[i + sh - (sh +4) + 4,sh + 4] += s[i,2]
         A22B[i,sh + 4] += s[i,2]
      end   

#     Column 2 and 3 
      @views @. r1 = JacVert.B1_23[:,:,iz,ID]
      @views @. r2 = JacVert.B2_23[:,:,iz,ID]
      #@views r3 = -invfac * (JacVert.A32[:,:,iz,ID] * r2 + FacGrav * r1[2:M-1,:]) 
      @inbounds for i = 1 : M2
        @views @. r3[i,:] = FacGrav * r1[i+1,:]
        @inbounds for j = 1 : M2
          @views @. r3[i,:] += JacVert.A32[i,j,iz,ID] * r2[j,:] 
        end  
        @views @. r3[i,:] *= -invfac
      end  

      @views ldivFull2!(JacVert.SA[:,:,iz,ID],r3)

      #r11 = invfac * (r1[1:1,:] - JacVert.A13[1:1,:,iz,ID] * r3[:,:])
      #r1M = invfac * (r1[M:M,:] - JacVert.A13[M:M,:,iz,ID] * r3[:,:])
      @views @. r11 = r1[1,:]
      @views @. r1M = r1[M,:]
      @inbounds for j = 1 : M2
        @views @. r11 += -JacVert.A13[1,j,iz,ID] * r3[j,:]
        @views @. r1M += -JacVert.A13[M,j,iz,ID] * r3[j,:]
      end  
      @views @. r11 *= invfac
      @views @. r1M *= invfac

      #@views r2 = invfac * (r2 - JacVert.A23[:,:,iz,ID] * r3)
      @inbounds for i = 1 : M2
        @inbounds for j = 1 : M2
          @views @. r2[i,:] += -JacVert.A23[i,j,iz,ID] * r3[j,:]
        end
        @views @. r2[i,:] *= invfac
      end  
        

      #@views s[[1,4],:] = -JacVert.C14_3[:,:,iz,ID] * r3
      #@views s[[2,3],:] = -JacVert.C23_2[:,:,iz,ID] * r2
      @. s = 0.0
      @inbounds for j = 1 : M2
        @views @. s[1,:] += - JacVert.C14_3[1,j,iz,ID] * r3[j,:]
        @views @. s[4,:] += - JacVert.C14_3[2,j,iz,ID] * r3[j,:]
        @views @. s[2,:] += - JacVert.C23_2[1,j,iz,ID] * r2[j,:]
        @views @. s[3,:] += - JacVert.C23_2[2,j,iz,ID] * r2[j,:]
      end  
        
      s[2,1] = s[2,1] - FacGrav * r11[1]
      s[3,1] = s[3,1] - FacGrav * r1M[1]
      s[2,2] = s[2,2] - FacGrav * r11[2]
      s[3,2] = s[3,2] - FacGrav * r1M[2]
      #@views A22B[sh + 1:sh + 4,[sh + 2,sh + 3]] .+= s
      @inbounds for i = 1 : 4
        #A22B[i + sh - (sh +2) + 4,sh + 2] += s[i,1]
        A22B[i + 2,sh + 2] += s[i,1]
        #A22B[i + sh - (sh +3) + 4,sh + 3] += s[i,2]
        A22B[i + 1,sh + 3] += s[i,2]
      end   
      if iz > 1
        #Column -2   
        #A22B[sh+2,sh-1] -= FacGrav * invfac * JacVert.B1m_34[1,iz,ID]
        #A22B[sh+2-(sh-1)+4,sh-1] -= FacGrav * invfac * JacVert.B1m_34[1,iz,ID]
        A22B[7,sh-1] -= FacGrav * invfac * JacVert.B1m_34[1,iz,ID]
        # Column -1
        #A22B[sh+2,sh] -= FacGrav * invfac * JacVert.B1m_34[2,iz,ID]
        #A22B[sh+2-sh+4,sh] -= FacGrav * invfac * JacVert.B1m_34[2,iz,ID]
        A22B[6,sh] -= FacGrav * invfac * JacVert.B1m_34[2,iz,ID]
      end
      if iz < nz 
#       Column +1  
        #A22B[sh+3,sh+5] -= FacGrav * invfac * JacVert.B1p_12[1,iz,ID]
        #A22B[sh+3-(sh+5)+4,sh+5] -= FacGrav * invfac * JacVert.B1p_12[1,iz,ID]
        A22B[2,sh+5] -= FacGrav * invfac * JacVert.B1p_12[1,iz,ID]
#       Column +2  
        #A22B[sh+3,sh+6] -= FacGrav * invfac * JacVert.B1p_12[2,iz,ID]
        #A22B[sh+3-(sh+6)+4,sh+6] -= FacGrav * invfac * JacVert.B1p_12[2,iz,ID]
        A22B[1,sh+6] -= FacGrav * invfac * JacVert.B1p_12[2,iz,ID]
      end
    end    
    luBand!(A22B,3,3)
  end
end

function ldivBlockAF(M,invfac,FacGrav,A13,A23,A32,SA,r1,r2,r3,r11,r1M)
   
  @inbounds for i = 1 : M - 2 
    r3[i] += -invfac * FacGrav * r1[i+1]
    @inbounds for j = 1 : M - 2  
      r3[i] += -invfac * A32[i,j] * r2[j]
    end
  end  

  ldivFull!(SA,r3)

  r11[1] = r1[1]
  r1M[1] = r1[M]
  @inbounds for j = 1 : M - 2
    r11[1] += -A13[1,j] * r3[j]  
    r1M[1] += -A13[M,j] * r3[j]  
  end  
  r11[1] *= invfac
  r1M[1] *= invfac

  @inbounds for i = 1 : M - 2
    @inbounds for j = 1 : M - 2  
      r2[i] += -A23[i,j] * r3[j]  
    end
    r2[i] *= invfac
  end  
end    

function ldivBlockAB(M,invfac,FacGrav,A13,A23,A32,SA,r1,r2,r3)

  #@views r3 .+= -invfac * (A32 * r2 + FacGrav * r1[2:end-1])
  @inbounds for i = 1 : M - 2 
    r3[i] += -invfac * FacGrav * r1[i+1]
    @inbounds for j = 1 : M - 2  
      r3[i] += -invfac * A32[i,j] * r2[j]
    end
  end  

  ldivFull!(SA,r3)

  #r1 .= invfac * (r1 - A13 * r3)
  @inbounds for i = 1 : M
    @inbounds for j = 1 : M - 2  
      r1[i] += -A13[i,j] * r3[j]  
    end
    r1[i] *= invfac
  end  
  #r2 .= invfac * (r2 - A23 * r3)
  @inbounds for i = 1 : M - 2
    @inbounds for j = 1 : M - 2  
      r2[i] += -A23[i,j] * r3[j]  
    end
    r2[i] *= invfac
  end  
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
  KldivVerticalFKernel!(JacVert,b,JacVert.rs,invfac,FacGrav,Val(M);ndrange=ndrange)

  group = (DoFG)
  ndrange = (DoF)
  KldivVerticalSKernel! = ldivVerticalSKernel!(backend,group)
  KldivVerticalSKernel!(JacVert.SchurBand,JacVert.rs;ndrange=ndrange)

  group = (nz, DoFG)
  ndrange = (nz, DoF)
  KldivVerticalBKernel! = ldivVerticalBKernel!(backend,group)
  KldivVerticalBKernel!(JacVert,b,JacVert.rs,invfac,FacGrav,Val(M);ndrange=ndrange)

end

function ldivVerticalOld!(JacVert,b)
  M = JacVert.M
  nz = JacVert.nz
  M2 = M - 2
  M1 = M - 1
  invfac = 1.0 / JacVert.fac
  FacGrav = JacVert.FacGrav
  r1 = zeros(M)
  r2 = zeros(M2)
  r3 = zeros(M2)
  rs = zeros(4*nz)
  s = zeros(4)
  r11 = zeros(1)
  r1M = zeros(1)
  RhoPos = 1
  ThPos = 5
  wPos = 4
  DoF = size(b,3)
  @inbounds for ID = 1 : DoF
    @views A22B = JacVert.SchurBand[:,:,ID]
#   Forward substitution  
    @inbounds for iz = 1 : nz
      sh = (iz - 1) * 4  
      rs[sh + 1] = b[1,iz,ID,ThPos]
      rs[sh + 2] = b[1,iz,ID,wPos]
      rs[sh + 3] = b[M,iz,ID,wPos]
      rs[sh + 4] = b[M,iz,ID,ThPos]
      @views @. r1 = b[:,iz,ID,RhoPos]
      @views @. r2 = b[2:M1,iz,ID,ThPos]
      @views @. r3 = b[2:M1,iz,ID,wPos]

      @views ldivBlockAF(M,invfac,FacGrav,JacVert.A13[:,:,iz,ID],JacVert.A23[:,:,iz,ID],JacVert.A32[:,:,iz,ID],
        JacVert.SA[:,:,iz,ID],r1,r2,r3,r11,r1M)
      @inbounds for j = 1 : M2
        rs[sh + 1] += -JacVert.C14_3[1,j,iz,ID] * r3[j]
        rs[sh + 4] += -JacVert.C14_3[2,j,iz,ID] * r3[j]
        rs[sh + 2] += -JacVert.C23_2[1,j,iz,ID] * r2[j]
        rs[sh + 3] += -JacVert.C23_2[2,j,iz,ID] * r2[j]
      end  
      rs[sh + 2] += -FacGrav * r11[1]
      rs[sh + 3] += -FacGrav * r1M[1]
    end    
    ldivBand!(A22B,rs,3,3)
#   Back substitution  
    @inbounds for iz = 1 : nz
      sh = (iz - 1) * 4
      b[1,iz,ID,ThPos] = rs[sh + 1]
      b[1,iz,ID,wPos] = rs[sh + 2]
      b[M,iz,ID,wPos] = rs[sh + 3]
      b[M,iz,ID,ThPos] = rs[sh + 4]
      @views @. r1 = b[:,iz,ID,RhoPos]
      @views @. r2 = b[2:M1,iz,ID,ThPos]
      @views @. r3 = b[2:M1,iz,ID,wPos]
      #@views rsC = rs[sh + 1 : sh + 4] 
      #r1[1] -= JacVert.B1_1[iz,ID] * rsC[1]
      #r1[M] -= JacVert.B1_4[iz,ID] * rsC[4]
      r1[1] -= JacVert.B1_1[iz,ID] * rs[sh + 1]
      r1[M] -= JacVert.B1_4[iz,ID] * rs[sh + 4]
      #@views r1 .-= JacVert.B1_23[:,:,iz,ID] * rsC[2:3]
      @inbounds for j = 1 : M
        r1[j] -= JacVert.B1_23[j,1,iz,ID] * rs[sh + 2] + JacVert.B1_23[j,2,iz,ID] * rs[sh + 3] 
      end  
      #@views r2 .-= JacVert.B2_23[:,:,iz,ID] * rsC[2:3]
      #@views r3 .-= JacVert.B3_14[:,:,iz,ID] * rsC[[1,4]]
      @inbounds for j = 1 : M2
        r2[j] -= JacVert.B2_23[j,1,iz,ID] * rs[sh + 2] + JacVert.B2_23[j,2,iz,ID] * rs[sh + 3] 
        r3[j] -= JacVert.B3_14[j,1,iz,ID] * rs[sh + 1] + JacVert.B3_14[j,2,iz,ID] * rs[sh + 4] 
      end  
      if iz > 1
        @views rsM = rs[sh - 1 : sh] 
        r1[1] -= JacVert.B1m_34[1,iz,ID] * rsM[1] + JacVert.B1m_34[2,iz,ID] * rsM[2]
      end
      if iz < nz
        @views rsP = rs[sh + 5 : sh + 6] 
        r1[M] -= JacVert.B1p_12[1,iz,ID] * rsP[1] + JacVert.B1p_12[2,iz,ID] * rsP[2]
      end
      @views ldivBlockAB(M,invfac,FacGrav,JacVert.A13[:,:,iz,ID],JacVert.A23[:,:,iz,ID],JacVert.A32[:,:,iz,ID],
        JacVert.SA[:,:,iz,ID],r1,r2,r3)
      @views @. b[:,iz,ID,RhoPos] = r1
      @views @. b[2:M1,iz,ID,ThPos] = r2
      @views @. b[2:M1,iz,ID,wPos] = r3
    end  
  end
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
