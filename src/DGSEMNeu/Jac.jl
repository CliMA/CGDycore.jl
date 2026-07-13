mutable struct JacDGVert{FT<:AbstractFloat,
                         AT2<:AbstractArray,
                         AT3<:AbstractArray,
                         AT4<:AbstractArray}
  CompTri::Bool                        
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

@inline function LUFull!(iz,ID,A,::Val{n}) where {n}

  @inbounds for j = 1 : n - 1
    @inbounds for i = j + 1 : n 
      A[i,j,iz,ID] /= A[j,j,iz,ID]  
      @inbounds for k = j + 1 : n
        A[i,k,iz,ID] -= A[i,j,iz,ID] * A[j,k,iz,ID]
      end  
    end  
  end  
end

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

@inline function ldivFull!(iz,ID,A,b,::Val{n}) where {n}

# Forward loop
  @inbounds for k = 1 : n - 1
    @inbounds for i = k + 1 : n
      b[i] -= A[i,k,iz,ID] * b[k]
    end
  end
#  Backward loop
  @inbounds for k = n : -1 : 1
    b[k] /= A[k,k,iz,ID]
    @inbounds for i = 1 : k - 1
      b[i] -= A[i,k,iz,ID] * b[k]
    end
  end
end

@inline function ldivBand!(ID,A,b,l,u,::Val{n}) where {n}

# Ly = b
  up1 = u + 1 
  uplp1 = up1 + l

  @inbounds for k = 1 : n - 1
    @inbounds for i = k + 1 : min(k+l,n)
      b[i,ID] -= A[i - k + up1,k,ID] * b[k,ID]  
    end    
  end
# Ux = y
  @inbounds for k = n : -1 : 1
    b[k,ID] /= A[up1,k,ID]
    @inbounds for i = max(k - u, 1) : k - 1
      b[i,ID] -= A[up1 + i - k,k,ID] * b[k,ID]  
    end    
  end
end

@inline function luBand!(ID,A,l,u,::Val{n}) where {n}
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
  up1 = u + 1
  up2 = u + 2
  uplp1 = u + l + 1
  @inbounds for k = 1 : n -1 
    @inbounds for i = up2 : uplp1 
      A[i,k,ID] /= A[up1,k,ID]  
    end  
    @inbounds for i = up2 : uplp1 
      @inbounds for j = k + 1 : min(k+u,n)
        A[i+k-j,j,ID] -= A[i,k,ID] * A[k + up1 - j,j,ID]
      end
    end  
  end    
end

function JacDGVert(backend,FT,M,nz,DG) 
  CompTri = false
  M2 = M - 2
  fac = 0
  FacGrav = 0
  NumI = DG.NumI
  A13 = KernelAbstractions.zeros(backend,FT,M,M2,nz,NumI)
  A23 = KernelAbstractions.zeros(backend,FT,M2,M2,nz,NumI)
  A32 = KernelAbstractions.zeros(backend,FT,M2,M2,nz,NumI)
  B1m_34 = KernelAbstractions.zeros(backend,FT,2,nz,NumI)
  B1_1 = KernelAbstractions.zeros(backend,FT,nz,NumI)
  B1_23 = KernelAbstractions.zeros(backend,FT,M,2,nz,NumI)
  B1_4 = KernelAbstractions.zeros(backend,FT,nz,NumI)
  B2_23 = KernelAbstractions.zeros(backend,FT,M2,2,nz,NumI)
  B3_14 = KernelAbstractions.zeros(backend,FT,M2,2,nz,NumI)
  B1p_12 = KernelAbstractions.zeros(backend,FT,2,nz,NumI)
  C23_2 = KernelAbstractions.zeros(backend,FT,2,M2,nz,NumI)
  C14_3 = KernelAbstractions.zeros(backend,FT,2,M2,nz,NumI)
  SA = KernelAbstractions.zeros(backend,FT,M2,M2,nz,NumI)
  SchurBand = KernelAbstractions.zeros(backend,FT,7,4*nz,NumI)
  rs = KernelAbstractions.zeros(backend,FT,4*nz,NumI)

  return JacDGVert{FT,
                   typeof(B1_1),
                   typeof(B1m_34),
                   typeof(A13)}(

    CompTri,
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

    ldivFull!(iz,ID,SA,r3,Val(M - 2))

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

@kernel inbounds = true function ldivVerticalSKernel!(A,rs,::Val{n}) where {n}

  ID, = @index(Global, NTuple)

  DoF = @uniform @ndrange()[1]

  if ID <= DoF
    ldivBand!(ID,A,rs,3,3,Val(n))
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


@kernel inbounds = true function luBandKernel!(A,::Val{n}) where {n}
  ID, = @index(Global, NTuple)

  DoF = @uniform @ndrange()[1]

  if ID <= DoF
    luBand!(ID,A,3,3,Val(n))
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

  kappa   = Phys.kappa
  kexp    = kappa / (eltype(U)(1) - kappa)
  kfac    = eltype(U)(1) / (eltype(U)(1) - kappa) * Phys.Rd

  if ID <= DoF
    sh = (iz - 1) * 4 + 1
    inv2dz = eltype(U)(2) / dz[iz,ID]
    invdz = eltype(U)(1) / dz[iz,ID]
    @unroll for i = 1 : M
      Th[i] = U[i,iz,ID,ThPos] / U[i,iz,ID,RhoPos]
      dpdRhoTh[i] = kfac * (Phys.Rd * U[i,iz,ID,ThPos] / Phys.p0)^kexp
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
    LUFull!(iz,ID,SA,Val(M - 2))

    if iz > 1
      B1m_34[1,iz,ID] = -eltype(U)(1) / (wB * dz[iz,ID])

      dpdRhoThM   = kfac * (Phys.Rd * U[M,iz-1,ID,ThPos] / Phys.p0)^kexp
      B1m_34[2,iz,ID] = -dpdRhoThM * invcS / (wB * dz[iz,ID])
    end
    @unroll for i = 1 : M
      B1_23[i,1,iz,ID] = eltype(U)(2) * DW[i,1] / dz[iz,ID]
      B1_23[i,2,iz,ID] = eltype(U)(2) * DW[i,M] / dz[iz,ID]
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
      dpdRhoThP   = kfac * (Phys.Rd * U[1,iz+1,ID,ThPos] / Phys.p0)^kexp
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
      ThM = U[M,iz-1,ID,ThPos] / U[M,iz-1,ID,RhoPos]
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

@kernel inbounds = true function FillJacFrozenDGVertKernel!(A13,A23,A32,B1m_34,B1_1,B1_23,B1_4,
  B2_23,B3_14,B1p_12,C23_2,C14_3,SA,SchurBand,@Const(Aux),@Const(dz),@Const(DW),@Const(w),fac,
  FacGrav,cS,Phys, ::Val{M}) where {M}

  iz,ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  DoF = @uniform @ndrange()[2]
  Th = @private eltype(Aux) (M,)
  dpdRhoTh = @private eltype(Aux) (M,)
  DWS = @localmem eltype(Aux) (M,M)
  SAL = @localmem eltype(Aux) (M,M)
  @uniform dpdThPos = 3
  @uniform ThPos = 4
  @uniform invcS = eltype(Aux)(1) / cS
  @uniform wB = w[1]

  if iz == 1
    @. DWS = DW   
  end
  @synchronize 

  kappa   = Phys.kappa
  kexp    = kappa / (eltype(Aux)(1) - kappa)
  kfac    = eltype(Aux)(1) / (eltype(Aux)(1) - kappa) * Phys.Rd

  if ID <= DoF
    sh = (iz - 1) * 4 + 1
    inv2dz = eltype(Aux)(2) / dz[iz,ID]
    invdz = eltype(Aux)(1) / dz[iz,ID]
    @unroll for i = 1 : M
      Th[i] = Aux[i,iz,ID,ThPos] 
      dpdRhoTh[i] = Aux[i,iz,ID,dpdThPos]
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
        SAL[i,j] = eltype(Aux)(0)
        @unroll for k = 1 : M - 2
          SAL[i,j] += A32[i,k,iz,ID] * A23[k,j,iz,ID]
        end
        if i == j
          SA[i,j,iz,ID] = fac - (FacGrav / fac) * A13[i+1,j,iz,ID] -
            (eltype(Aux)(1) / fac) * SAL[i,j]
        else
          SA[i,j,iz,ID] = -(FacGrav / fac) * A13[i+1,j,iz,ID] -
            (eltype(Aux)(1) / fac) * SAL[i,j]
        end
      end
    end
    LUFull!(iz,ID,SA,Val(M - 2))

    if iz > 1
      B1m_34[1,iz,ID] = -eltype(Aux)(1) / (wB * dz[iz,ID])

      dpdRhoThM   = Aux[M,iz-1,ID,dpdThPos]
      B1m_34[2,iz,ID] = -dpdRhoThM * invcS / (wB * dz[iz,ID])
    end
    @unroll for i = 1 : M
      B1_23[i,1,iz,ID] = eltype(Aux)(2) * DW[i,1] / dz[iz,ID]
      B1_23[i,2,iz,ID] = eltype(Aux)(2) * DW[i,M] / dz[iz,ID]
    end  

    if iz == 1
      B1_1[iz,ID] = eltype(Aux)(0)
    else
      B1_23[1,1,iz,ID] = eltype(Aux)(0)
      B1_1[iz,ID] = dpdRhoTh[1] * invcS * invdz / wB 
    end
    if iz == nz
      B1_4[iz,ID] = eltype(Aux)(0)
    else
      B1_23[M,2,iz,ID] = eltype(Aux)(0)
      B1_4[iz,ID] = dpdRhoTh[M] * invcS * invdz / wB 
    end
    if iz < nz
      dpdRhoThP   = Aux[1,iz+1,ID,dpdThPos]
      B1p_12[1,iz,ID] = -dpdRhoThP * invcS / (wB * dz[iz,ID])
      B1p_12[2,iz,ID] = eltype(Aux)(1) / (wB * dz[iz,ID])
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
      ThM = Aux[M,iz-1,ID,ThPos]
      dpdRhoThM = Aux[M,iz-1,ID,dpdThPos]
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

function FillJacDGVert!(JacVert::JacDGVert,U,DG,dz,fac,Phys)

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
  @. JacVert.SchurBand = FTB(0)
  @views @. JacVert.SchurBand[4,:,:] = fac
  KFillJacDGVertKernel! = FillJacDGVertKernel!(backend,group)
  KFillJacDGVertKernel!(JacVert.A13,JacVert.A23,JacVert.A32,JacVert.B1m_34,JacVert.B1_1,
  JacVert.B1_23,JacVert.B1_4, JacVert.B2_23,JacVert.B3_14,JacVert.B1p_12,JacVert.C23_2,
  JacVert.C14_3,JacVert.SA,JacVert.SchurBand,U,dz,DWZ,DG.wZ,fac,JacVert.FacGrav,
  Phys.cS,Phys,Val(M);ndrange=ndrange) 

end  

function FillJacFrozenDGVert!(JacVert::JacDGVert,Aux,DG,dz,fac,Phys)

  backend = get_backend(Aux)
  FTB = eltype(Aux)

  M = JacVert.M
  nz = JacVert.nz
  DoF  = DG.NumI

  JacVert.fac = fac
  JacVert.FacGrav = Phys.Grav

  DWZ = DG.DWZ

  DoFG = 10
  group = (nz, DoFG)
  ndrange = (nz, DoF) 
  @. JacVert.SchurBand = FTB(0)
  @views @. JacVert.SchurBand[4,:,:] = fac
  KFillJacDGVertKernel! = FillJacFrozenDGVertKernel!(backend,group)
  KFillJacDGVertKernel!(JacVert.A13,JacVert.A23,JacVert.A32,JacVert.B1m_34,JacVert.B1_1,
  JacVert.B1_23,JacVert.B1_4, JacVert.B2_23,JacVert.B3_14,JacVert.B1p_12,JacVert.C23_2,
  JacVert.C14_3,JacVert.SA,JacVert.SchurBand,Aux,dz,DWZ,DG.wZ,fac,JacVert.FacGrav,
  Phys.cS,Phys,Val(M);ndrange=ndrange) 

end  
function SchurBoundary!(JacVert::JacDGVert)

  backend = get_backend(JacVert.SchurBand)
  FTB = eltype(JacVert.SchurBand)

  M = JacVert.M
  nz = JacVert.nz
  DoF = size(JacVert.SchurBand,3)

  DoFG = 1
  group = (nz, DoFG)
  ndrange = (nz, DoF)
  KSchurBoundaryKernel! = SchurBoundaryKernel!(backend,group)
  KSchurBoundaryKernel!(JacVert.A13,JacVert.A23,JacVert.A32,JacVert.B1m_34,JacVert.B1_1,
    JacVert.B1_23,JacVert.B1_4, JacVert.B2_23,JacVert.B3_14,JacVert.B1p_12,JacVert.C23_2,
    JacVert.C14_3,JacVert.SA,JacVert.SchurBand,JacVert.fac,JacVert.FacGrav,Val(M);ndrange=ndrange)
  group = (DoFG)
  ndrange = (DoF)
  KluBandKernel! = luBandKernel!(backend,group)
  KluBandKernel!(JacVert.SchurBand,Val(4*nz),ndrange=ndrange)

end  

function Solve!(JacVert::JacDGVert,b)

  backend = get_backend(JacVert.SchurBand)
  FTB = eltype(JacVert.SchurBand)

  M = JacVert.M
  nz = JacVert.nz
  DoF = size(JacVert.SchurBand,3)

  invfac = FTB(1) / JacVert.fac
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
  KldivVerticalSKernel!(JacVert.SchurBand,JacVert.rs,Val(4*nz);ndrange=ndrange)

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
    @inbounds for iv = [2 3] 
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
    p[ii] = 1 + (iz-1) * M  + (ivw - 1) * N
    ii += 1
    p[ii] = 1 + (iz-1) * M  + (ivTh - 1) * N
    ii += 1
    p[ii] = M + (iz - 1) * M + (ivTh - 1) * N
    ii += 1
    p[ii] = M + (iz - 1) * M + (ivw - 1) * N
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

function Jac!(U,fac,DG,Metric,Phys,Cache,JCache::JacDGVert,Global,VelForm)
  Invfac = eltype(U)(1) / fac
  dz = Metric.dz
  FillJacDGVert!(JCache,U,DG,dz,Invfac,Phys)
  SchurBoundary!(JCache)
end  

function JacFrozen!(Aux,fac,DG,Metric,Phys,Cache,JCache::JacDGVert,Global,VelForm)
  Invfac = eltype(Aux)(1) / fac
  dz = Metric.dz
  FillJacFrozenDGVert!(JCache,Aux,DG,dz,Invfac,Phys)
  SchurBoundary!(JCache)
end  



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

mutable struct JacobianCacheMarcoSplitNS{NZ, M, RealT, AType2 <: AbstractArray,
                             AType3 <: AbstractArray, AType4 <: AbstractArray}
    CompTri::Bool                         
    M::Int
    nz::Int
    NumI::Int
    fac::RealT
    invwB::RealT
    schur_dinv :: RealT      # = inv(fac) = γΔt  (D11 = fac always; scalar)
    A12_11::AType2
    A22_diag::AType3
    A33_diag::AType3
    DWS::AType2
    Th_cache::AType3
    dpdRhoTh_cache::AType3
    phi_cache::AType3        # (M, nz, NumI) geopotential (frozen) — for split gravity
    Gblk::AType4             # (M, M, nz, NumI) gravity block A31 = G[i,k] (constant)
    grav_done::Base.RefValue{Bool}   # lazy precompute flag for Gblk
    SA::AType4
    SchurD  :: AType4
    SchurL  :: AType4
    SchurU  :: AType4
    Schurm      :: AType2    # (nz, NumI) per-element ρ-elimination multiplier D21/fac
    SchurLhat   :: AType4
    SchurDinv   :: AType4
    SchurUhat   :: AType4
    rs      :: AType2
end

@inline nvelem(::JacobianCacheMarcoSplitNS{NZ, M}) where {NZ, M} = NZ
@inline nvnodes(::JacobianCacheMarcoSplitNS{NZ, M}) where {NZ, M} = M

function JacobianCacheMarcoSplitNS(backend,FT,M,nz,DG) 
    CompTri = false
    polydeg = M - 1
    NumI = DG.NumI
    fac = 1.0
    schur_dinv = inv(fac)
    invwb = inv(DG.wZ[1])
    A12_11 = KernelAbstractions.zeros(backend,FT, nz, NumI)
    A22_diag = KernelAbstractions.zeros(backend,FT, polydeg, nz, NumI)
    A33_diag = KernelAbstractions.zeros(backend,FT, polydeg, nz, NumI)
    DWS = KernelAbstractions.zeros(backend,FT,M,M)
    copyto!(DWS,DG.DWZ)
    Th_cache        = KernelAbstractions.zeros(backend,FT, M, nz, NumI)
    dpdRhoTh_cache  = KernelAbstractions.zeros(backend,FT, M, nz, NumI)
    phi_cache       = KernelAbstractions.zeros(backend,FT, M, nz, NumI)
    Gblk            = KernelAbstractions.zeros(backend,FT, M, M, nz, NumI)
    grav_done       = Ref(false)
    SA = KernelAbstractions.zeros(backend,FT, polydeg, polydeg, nz, NumI)
    SchurD = KernelAbstractions.zeros(backend,FT, 3, 3, nz, NumI)
    # Full 3×3 block-tridiagonal (the split gravity makes ϱ_M couple densely,
    # so the pointwise 2×2 reduction is invalid — see memory note).
    SchurL = KernelAbstractions.zeros(backend,FT, 3, 3, nz - 1, NumI)
    SchurU = KernelAbstractions.zeros(backend,FT, 3, 3, nz - 1, NumI)
    Schurm = KernelAbstractions.zeros(backend,FT, nz, NumI)             # unused by the 3×3 Thomas (kept for struct compat)
    Schurm .= P.Grav / fac
    SchurLhat = KernelAbstractions.zeros(backend,FT, 3, 3, nz - 1, NumI)
    SchurDinv = KernelAbstractions.zeros(backend,FT, 3, 3, nz, NumI)
    SchurUhat = KernelAbstractions.zeros(backend,FT, 3, 3, nz - 1, NumI)
    rs = zeros(FT, 3 * nz, NumI)

    invwB = inv(DG.wZ[1])
#   zCol = compute_vertical_step_size(NumI, polydeg, semi, nx, nz, typeof(semi.mesh))
#   wPos, ThPos = compute_variable_position(equations)
#   jacobian_aux = [JacobianAuxMarcoNS(M, FT) for _ in 1:Threads.nthreads()]
#   rotation_matrix = pre_compute_rotation_matrix(semi, equations)
#   ID_map, iz_map = build_vertical_mapping(semi, nx, nz, equations)
#   u_vert = zeros(FT, M, nz, NumI, nvariables(semi))

    return JacobianCacheMarcoSplitNS{nz, M, FT, typeof(A12_11), typeof(Th_cache),
                         typeof(SA)}(CompTri, M, nz, NumI, fac, invwb, schur_dinv, 
                                         A12_11, A22_diag, A33_diag,
                                         DWS, Th_cache, dpdRhoTh_cache, phi_cache,
                                         Gblk, grav_done,
                                         SA, SchurD, SchurL, SchurU,
                                         Schurm, SchurLhat, SchurDinv, SchurUhat,
                                         rs)
end

# Lazy precompute of the (state-independent) gravity block G[i,k] per element.
#   :pointwise  -> g·I
#   :split      -> G_ik = inv2dz[ 0.5 D_ik (φ_k-φ_i) + δ_ik 0.5 Σ_m D_im (φ_m-φ_i) ]
# Requires phi_cache (filled by rotate_wrap_to_vertical!).  Run once.
function precompute_gravity!(phi_cache,dz, cache::JacobianCacheMarcoSplitNS)
  cache.grav_done[] && return
  (; M, nz, NumI, DWS, Gblk) = cache
  @inbounds for ID in 1:NumI
    @inbounds for iz in 1:nz
      inv2dz = 2 / dz[iz, ID]
      @inbounds for i in 1:M
        acc = zero(eltype(Gblk))
        @inbounds for k in 1:M
          dphi = phi_cache[k, iz, ID] - phi_cache[i, iz, ID]
          Gblk[i, k, iz, ID] = inv2dz * 0.5 * DWS[i, k] * dphi
          acc += DWS[i, k] * dphi
        end
        Gblk[i, i, iz, ID] += inv2dz * 0.5 * acc
      end
    end
  end
  cache.grav_done[] = true
  return nothing
end

function precompute_gravity!(phi_cache,dz, cache::JacobianCacheMarcoSplitNS, NumberThreadGPU)
  cache.grav_done[] && return

  (; M, nz, NumI, DWS, Gblk) = cache
  backend = get_backend(Gblk)
  NumIG = min(div(NumberThreadGPU,M*nz),NumI)
  group = (M,nz,NumIG)
  ndrange = (M,nz,NumI) 
  Kprecompute_gravityKernel! = precompute_gravityKernel!(backend,group)
  Kprecompute_gravityKernel!(Gblk,phi_cache,DWS,dz,Val(M),ndrange=ndrange)
  cache.grav_done[] = true
end  
  

@kernel inbounds = true function precompute_gravityKernel!(Gblk,@Const(Geo),@Const(DWS),@Const(dz), ::Val{M}) where {M}

  i,iz,ID = @index(Global, NTuple)

  NumI = @uniform @ndrange()[3]
  if ID <= NumI
    acc = eltype(Gblk)(0)
    Geoi = Geo[i, iz, ID]
    inv2dz = eltype(Gblk)(2) / dz[iz, ID]
    @unroll for k in 1:M
      dphi = Geo[k, iz, ID] - Geoi
      Gblk[i, k, iz, ID] = inv2dz * 0.5 * DWS[i, k] * dphi 
#     Gblk[i, k, iz, ID] = 0
      acc += DWS[i, k] * dphi 
    end  
    Gblk[i, i, iz, ID] += inv2dz * 0.5 * acc
#   Gblk[i, i, iz, ID] = sum(Gblk[i, :, iz, ID])
  end  
end


#
#function FillJacDGVertKernel!(cache_jacobian::JacobianCacheMarcoSplitNS, u, dz, gamdt)
#  (; A12_11, A22_diag, A33_diag, SA, SchurD, SchurL, SchurU, Schurm, M, nz, NumI,
#     fac, DWS, Th_cache, dpdRhoTh_cache, Gblk, invwB) = cache_jacobian
#     
function FillJacDGVert!(cache_jacobian::JacobianCacheMarcoSplitNS, u, dz, fac,NumberThreadGPU)
  (; A12_11, A22_diag, A33_diag, SA, SchurD, SchurL, SchurU, Schurm, M, nz, NumI,
     fac, DWS, Th_cache, dpdRhoTh_cache, Gblk, invwB) = cache_jacobian
  cache_jacobian.fac = fac
  backend = get_backend(u)

  DoFG = 10
  group = (nz,DoFG)
  ndrange = (nz,NumI) 

  KFillJacDGVertMKernel! = FillJacDGVertMKernel!(backend,group)
  KFillJacDGVertMKernel!(u,dz,
    A12_11, A22_diag, A33_diag, SA, SchurD, SchurL, SchurU, Schurm, 
    cache_jacobian.fac, DWS, Th_cache, dpdRhoTh_cache, Gblk, invwB,Val(M),ndrange=ndrange)
end  



@kernel inbounds = true function FillJacDGVertMKernel!(@Const(U),@Const(dz),
A12_11, A22_diag, A33_diag, SA, SchurD, SchurL, SchurU, Schurm, 
     fac, @Const(DW), Th_cache, dpdRhoTh_cache, @Const(Gblk), invwB, ::Val{M}) where {M}
  
  iz,ID = @index(Global, NTuple)
  
  nz = @uniform @ndrange()[1]
  NumI = @uniform @ndrange()[2]

  Th = @private eltype(U) (M,)
  dpdRhoTh = @private eltype(U) (M,)
  DWS = @private eltype(U) (M,M)
  SAL1 = @private eltype(U) (M-1,M-1)
  SAL2 = @private eltype(U) (M-1,M-1)
  @uniform RhoPos = 1
  @uniform ThPos = 5
  

  if ID <= NumI
    @. DWS = DW
     
    FacGrav = P.Grav   
    kappa   = P.kappa
    cS = P.cS
    invcS = eltype(U)(1) / cS
    kexp    = kappa / (eltype(U)(1) - kappa)
    kfac    = eltype(U)(1) / (eltype(U)(1) - kappa) * P.Rd
    pre_fac = kfac * (P.Rd / P.p0)^kexp
    invfac  = inv(fac)
    A22_diag[:,iz,ID] .= invfac
    A33_diag[:,iz,ID] .= fac
    inv2dz = eltype(U)(2) / dz[iz, ID]
    invdz  = eltype(U)(1) / dz[iz, ID]

    @unroll for i in 1:M
      Th[i]       = U[i, iz, ID, ThPos] / U[i, iz, ID, RhoPos]
      dpdRhoTh[i] = pre_fac * U[i, iz, ID, ThPos]^kexp
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

    @unroll for j in 1:(M - 1)
      thj_inv2dz = Th[j] * inv2dz
      @unroll for k in 1:(M - 1)
        SAL1[k, j] = DWS[k, j] * thj_inv2dz * A22_diag[k, iz, ID]
      end
    end
    if !is_iz1
      SAL1[1, 1] = eltype(U)(0)
    end

    # SA gravity term: -invfac·(G_II·A13),  A13[l,j]=inv2dz·DWS[l,j], A13[1,1]=0 (iz>1)
    @unroll for j in 1:(M - 1)
      @unroll for i in 1:(M - 1)
        SAL2[i, j] = eltype(U)(0)
      end
    end
    @unroll for j in 1:(M - 1)
      @unroll for l in 1:(M - 1)
        a13lj = inv2dz * DWS[l, j] * invfac
        @unroll for i in 1:(M - 1)
          SAL2[i, j] -= Gblk[i, l, iz, ID] * a13lj
        end
      end
      SAL2[j, j] += A33_diag[j, iz, ID]
    end
    if !is_iz1
      # A13[1,1]=0 for iz>1: remove the (l=1,j=1) contribution G[i,1]·A13[1,1]
      corr13 = inv2dz * dws11 * invfac
      @unroll for i in 1:(M - 1)
        SAL2[i, 1] += Gblk[i, 1, iz, ID] * corr13
      end
    end

    # SA -= A32·A22⁻¹·A23
    @unroll for j in 1:(M - 1)
      @unroll for k in 1:(M - 1)
        sal_kj      = SAL1[k, j]
        dpdk_inv2dz = dpdRhoTh[k] * inv2dz
        @unroll for i in 1:(M - 1)
          SAL2[i, j] -= DWS[i, k] * dpdk_inv2dz * sal_kj
        end
      end
    end
    _corr32 = (is_iz1 ? 2 : 1) * dws11 * dpdRhoTh[1] * inv2dz
    @unroll for j in 1:(M - 1)
      SAL2[1, j] += _corr32 * SAL1[1, j]
    end
    # A12 cross-term: +invfac·A12_11·G_II[:,1]·SAL[1,:]  (column G[:,1], all rows)
    a12scale = invfac * A12_11[iz, ID]
    @unroll for j in 1:(M - 1)
      sal1j = SAL1[1, j]
      @unroll for i in 1:(M - 1)
        SAL2[i, j] += a12scale * Gblk[i, 1, iz, ID] * sal1j
      end
    end
    LUFull!(SAL2, Val(M - 1))
    @unroll for j in 1:(M - 1)
      @unroll for i in 1:(M - 1)
        SA[i, j, iz, ID] = SAL2[i, j]
      end
    end

    # ── D block ──────────────────────────────────────────────────
    D11 = fac
    D21 = Gblk[M, M, iz, ID]      # boundary gravity ρw_M ← ρ_M (=g pointwise)
    D31 = eltype(U)(0)
    if iz < nz
      D12 = eltype(U)(0)
      D22 = fac + cS * invdz * invwB
      D32 = eltype(U)(0)
      D13 = dpdRhoTh[M] * invcS * invdz * invwB
      D23 = eltype(U)(0)
      D33 = fac + Th[M] * dpdRhoTh[M] * invcS * invdz * invwB
    else
      D12 = eltype(U)(2) * DWS[M, M] / dz[iz, ID]
      D22 = fac + inv2dz * cS * invwB
      D32 = inv2dz * DWS[M, M] * Th[M]
      D13 = eltype(U)(0)
      D23 = -inv2dz * DWS[M, M] * dpdRhoTh[M]
      D33 = fac
    end
    SchurD[1, 1, iz, ID] = D11; SchurD[2, 1, iz, ID] = D21; SchurD[3, 1, iz, ID] = D31
    SchurD[1, 2, iz, ID] = D12; SchurD[2, 2, iz, ID] = D22; SchurD[3, 2, iz, ID] = D32
    SchurD[1, 3, iz, ID] = D13; SchurD[2, 3, iz, ID] = D23; SchurD[3, 3, iz, ID] = D33
    Schurm[iz, ID] = D21 * invfac            # per-element ρ-elim multiplier

    if iz < nz
      f0 = eltype(U)(0)
      @unroll for c in 1:3
        SchurL[1, c, iz, ID] = f0; SchurL[2, c, iz, ID] = f0; SchurL[3, c, iz, ID] = f0
        SchurU[1, c, iz, ID] = f0; SchurU[2, c, iz, ID] = f0; SchurU[3, c, iz, ID] = f0
      end
    end
  end
end

function FillJacDGVert!(cache_jacobian::JacobianCacheMarcoSplitNS, u, dz, gamdt)
  (; A12_11, A22_diag, A33_diag, SA, SchurD, SchurL, SchurU, Schurm, M, nz, NumI,
     fac, DWS, Th_cache, dpdRhoTh_cache, Gblk, invwB) = cache_jacobian
  cache_jacobian.fac = gamdt
  FillJacDGVertKernel!(cache_jacobian, u, dz)
end


function FillJacDGVertKernel!(cache_jacobian::JacobianCacheMarcoSplitNS, u, dz)
  (; A12_11, A22_diag, A33_diag, SA, SchurD, SchurL, SchurU, Schurm, M, nz, NumI,
     fac, DWS, Th_cache, dpdRhoTh_cache, Gblk, invwB) = cache_jacobian
  kappa   = P.kappa
  cS = P.cS
  invcS   = 1 / cS
  kexp    = kappa / (eltype(u)(1) - kappa)
  kfac    = eltype(u)(1) / (eltype(u)(1) - kappa) * P.Rd
  pre_fac = kfac * (P.Rd / P.p0)^kexp
  RhoPos  = 1
  wPos = 4
  ThPos = 5
  invfac  = inv(fac)

  Th = zeros(M)
  dpdRhoTh = zeros(M)
  SAL = zeros(M-1,M-1)

  @inbounds for ID in 1:NumI

    @inbounds for iz in 1:nz
      inv2dz = eltype(u)(2) / dz[iz, ID]
      invdz  = eltype(u)(1) / dz[iz, ID]
      A22_diag[:,iz,ID] .= invfac
      A33_diag[:,iz,ID] .= fac

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
        @inbounds for k in 1:(M - 1)
          SAL[k, j] = DWS[k, j] * thj_inv2dz * A22_diag[k, iz, ID]
        end
      end
      if !is_iz1
        SAL[1, 1] = zero(eltype(u))
      end

      # SA gravity term: -invfac·(G_II·A13),  A13[l,j]=inv2dz·DWS[l,j], A13[1,1]=0 (iz>1)
      @inbounds for j in 1:(M - 1)
        @inbounds for i in 1:(M - 1)
          SA[i, j, iz, ID] = zero(eltype(u))
        end
      end
      @inbounds for j in 1:(M - 1)
        for l in 1:(M - 1)
          a13lj = inv2dz * DWS[l, j] * invfac
          @inbounds for i in 1:(M - 1)
            SA[i, j, iz, ID] -= Gblk[i, l, iz, ID] * a13lj
          end
        end
        SA[j, j, iz, ID] += A33_diag[j, iz, ID]
      end
      if !is_iz1
        # A13[1,1]=0 for iz>1: remove the (l=1,j=1) contribution G[i,1]·A13[1,1]
        corr13 = inv2dz * dws11 * invfac
        @inbounds for i in 1:(M - 1)
          SA[i, 1, iz, ID] += Gblk[i, 1, iz, ID] * corr13
        end
      end

      # SA -= A32·A22⁻¹·A23
      @inbounds for j in 1:(M - 1)
        @inbounds for k in 1:(M - 1)
          sal_kj      = SAL[k, j]
          dpdk_inv2dz = dpdRhoTh[k] * inv2dz
          @inbounds for i in 1:(M - 1)
            SA[i, j, iz, ID] -= DWS[i, k] * dpdk_inv2dz * sal_kj
          end
        end
      end
      _corr32 = (is_iz1 ? 2 : 1) * dws11 * dpdRhoTh[1] * inv2dz
      @inbounds for j in 1:(M - 1)
        SA[1, j, iz, ID] += _corr32 * SAL[1, j]
      end
      # A12 cross-term: +invfac·A12_11·G_II[:,1]·SAL[1,:]  (column G[:,1], all rows)
      a12scale = invfac * A12_11[iz, ID]
      @inbounds for j in 1:(M - 1)
        sal1j = SAL[1, j]
        @inbounds for i in 1:(M - 1)
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
        @inbounds for c in 1:3
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
  # R3 -= invfac · G_II · bρ   (matvec; pointwise = -fgi·bρ)
  for l in 1:polydeg
    brl = R1[l] * invfac
    @inbounds for i in 1:polydeg
      R3[i] -= Gblk[i, l, iz, ID] * brl
    end
  end
  # R3 -= A32 · A22⁻¹ · bθ
  for j in 1:polydeg
    r2js = R2[j] * A22_diag[j, iz, ID]
    dj   = dpdRhoTh_cache[j, iz, ID] * inv2dz
    @inbounds for i in 1:polydeg
      R3[i] -= DWS[i, j] * dj * r2js
    end
  end
  R3[1] += (is_iz1 ? 2 : 1) * DWS[1, 1] * dpdRhoTh_cache[1, iz, ID] * inv2dz *
            R2[1] * A22_diag[1, iz, ID]
  # A12 cross:  R3 += invfac·a12·G[:,1]·(A22⁻¹ bθ)_1
  a12c = invfac * a12 * R2[1] * A22_diag[1, iz, ID]
  @inbounds for i in 1:polydeg
    R3[i] += a12c * Gblk[i, 1, iz, ID]
  end
  return nothing
end

# Phase C (one column): recover ϱθ (R2) and ϱ (R1) after the SA solve has overwritten R3.
#   On entry:  R1 = bρ, R2 = bθ, R3 = SA-solved momentum.  On exit R1,R2 = responses.
@inline function schur_postsa_col!(R1, R2, R3, iz, ID, A22_diag, DWS, Th_cache,
                                   a12, inv2dz, invfac, is_iz1, polydeg)
  # ϱθ recovery:  R2 = A22⁻¹·(bθ - A23·R3)
  for j in 1:polydeg
    r3j = R3[j];  tj = Th_cache[j, iz, ID] * inv2dz
    @inbounds for i in 1:polydeg
      R2[i] -= DWS[i, j] * tj * r3j
    end
  end
  if !is_iz1
    R2[1] += DWS[1, 1] * Th_cache[1, iz, ID] * inv2dz * R3[1]
  end
  @inbounds for i in 1:polydeg
    R2[i] *= A22_diag[i, iz, ID]
  end
  # ϱ recovery:  R1 = invfac·(bρ - a12·R2[1]·e1 - A13·R3)   (R1 still holds bρ)
  for j in 1:polydeg
    r3j = R3[j]
    @inbounds for i in 1:polydeg
      R1[i] -= inv2dz * DWS[i, j] * r3j
    end
  end
  R1[1] -= a12 * R2[1]
  if !is_iz1
    R1[1] += inv2dz * DWS[1, 1] * R3[1]
  end
  @inbounds for i in 1:polydeg
    R1[i] *= invfac
  end
  return nothing
end

# Batched LU triangular solve of the SA factor over the 3 driver columns of r3.
@inline function ldivFull3!(iz, ID, SA, r3, ::Val{n}) where {n}
  @inbounds begin
    @inbounds for k in 1:(n - 1)
      @inbounds for i in (k + 1):n
        sa = SA[i, k, iz, ID]
        r3[i, 1] -= sa * r3[k, 1];  r3[i, 2] -= sa * r3[k, 2];  r3[i, 3] -= sa * r3[k, 3]
      end
    end
    @inbounds for k in n:-1:1
      sak = inv(SA[k, k, iz, ID])
      r3[k, 1] *= sak;  r3[k, 2] *= sak;  r3[k, 3] *= sak
      @inbounds for i in 1:(k - 1)
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
    @inbounds for k in 1:(n - 1)
      @inbounds for i in (k + 1):n
        sa = SA[i, k, iz, ID]
        r3[i, 1] -= sa * r3[k, 1];  r3[i, 2] -= sa * r3[k, 2]
      end
    end
    @inbounds for k in n:-1:1
      sak = inv(SA[k, k, iz, ID])
      r3[k, 1] *= sak;  r3[k, 2] *= sak
      @inbounds for i in 1:(k - 1)
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
function SchurBoundary!(cache::JacobianCacheMarcoSplitNS, dz, Phys)
  (; M, nz, NumI, A12_11, A22_diag, DWS, Th_cache, dpdRhoTh_cache, Gblk,
     SA, SchurD, SchurL, SchurU, fac, invwB) = cache
  FT      = eltype(SA)
  invfac  = inv(fac)
  cS = P.cS
  invcS   = inv(cS)
  polydeg = M - 1

  r1  = zeros(FT, M - 1, 3)   # 3 driver columns (ϱ_M, ϱw_M, ϱθ_M)
  r2  = zeros(FT, M - 1, 3)
  r3  = zeros(FT, M - 1, 3)
  Th  = zeros(FT, M)
  dpdRhoTh = zeros(FT, M)
  SAL = zeros(FT, M - 1, M - 1)
  @inbounds for ID in 1:NumI
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
      @inbounds for i in 1:polydeg
        r1[i,1]=zero(FT); r1[i,2]=zero(FT); r1[i,3]=zero(FT)
        r2[i,1]=zero(FT); r2[i,2]=zero(FT); r2[i,3]=zero(FT)
        r3[i,1]=zero(FT); r3[i,2]=zero(FT); r3[i,3]=zero(FT)
      end
      @inbounds for i in 1:polydeg
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
      ldivFull3!(iz, ID, SA, r3, Val(M - 1))
      @inbounds for col in 1:3
        schur_postsa_col!(view(r1,:,col), view(r2,:,col), view(r3,:,col), iz, ID,
                          A22_diag, DWS, Th_cache, a12, inv2dz, invfac, is_iz1, polydeg)
        # C_self → SchurD[:,col,iz]
        p1 = zero(FT);  p2 = zero(FT);  p3 = zero(FT);  pg = zero(FT)
        @inbounds for j in 1:polydeg
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
        @inbounds for i in 1:polydeg
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
        ldivFull2!(iz, ID, SA, r3, Val(M - 1))
        @inbounds for (slot, col) in ((1, 2), (2, 3))   # slot1→col2 (ϱw_M), slot2→col3 (ϱθ_M)
          schur_postsa_col!(view(r1,:,slot), view(r2,:,slot), view(r3,:,slot), iz, ID,
                            A22_diag, DWS, Th_cache, a12, inv2dz, invfac, is_iz1, polydeg)
          # C_p → SchurD[:,col,iz-1] self correction (node-1 response)
          r2_1 = r2[1,slot];  r3_1 = r3[1,slot]
          SchurD[1, col, iz-1, ID] -= c2p1*r2_1 + c3p1*r3_1
          SchurD[2, col, iz-1, ID] -= c2p2*r2_1 + c3p2*r3_1
          SchurD[3, col, iz-1, ID] -= c2p3*r2_1 + c3p3*r3_1
          # C_self → SchurL[:,col,iz-1]  (full interior(iz) response + Cgrav)
          p1 = zero(FT);  p2 = zero(FT);  p3 = zero(FT);  pg = zero(FT)
          @inbounds for j in 1:polydeg
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

#=
function SchurBoundary!(cache::JacobianCacheMarcoSplitNS, dz, NumberThreadGPU::Int64)
  (; M, nz, NumI, A12_11, A22_diag, DWS, Th_cache, dpdRhoTh_cache, Gblk,
     SA, SchurD, SchurL, SchurU, fac, invwB) = cache
  FT      = eltype(SA)
  invfac  = inv(fac)
  cS = P.cS
  invcS   = inv(cS)
  polydeg = M - 1
end  

@kernel inbounds=true function SchurBoundaryKernel!(cache::JacobianCacheMarcoSplitNS, dz, Phys)
  (; M, nz, NumI, A12_11, A22_diag, DWS, Th_cache, dpdRhoTh_cache, Gblk,
     SA, SchurD, SchurL, SchurU, fac, invwB) = cache
  FT      = eltype(SA)
  iz,ID = @index(Global, NTuple)
  NumI = @uniform @ndrange()[2]

  invfac  = inv(fac)
  cS = P.cS
  invcS   = inv(cS)
  polydeg = M - 1

  r1 = @private FT (M - 1, 3)
  r2 = @private FT (M - 1, 3)
  r3 = @private FT (M - 1, 3)
  Th = @private FT (M - 1)
  dpdRhoTh = @private FT (M - 1)
  SAL = @private FT (M - 1, M -1)

  if ID <= NumI
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
    @unroll for i in 1:polydeg
      r1[i,1]=zero(FT); r1[i,2]=zero(FT); r1[i,3]=zero(FT)
      r2[i,1]=zero(FT); r2[i,2]=zero(FT); r2[i,3]=zero(FT)
      r3[i,1]=zero(FT); r3[i,2]=zero(FT); r3[i,3]=zero(FT)
    end
    @unroll for i in 1:polydeg
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
    ldivFull3!(iz, ID, SA, r3, Val(M - 1))
    @inbounds for col in 1:3
      schur_postsa_col!(view(r1,:,col), view(r2,:,col), view(r3,:,col), iz, ID,
                        A22_diag, DWS, Th_cache, a12, inv2dz, invfac, is_iz1, polydeg)
      # C_self → SchurD[:,col,iz]
      p1 = zero(FT);  p2 = zero(FT);  p3 = zero(FT);  pg = zero(FT)
      @inbounds for j in 1:polydeg
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
      @unroll for i in 1:polydeg
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
      ldivFull2!(iz, ID, SA, r3, Val(M - 1))
      @unroll for (slot, col) in ((1, 2), (2, 3))   # slot1→col2 (ϱw_M), slot2→col3 (ϱθ_M)
        schur_postsa_col!(view(r1,:,slot), view(r2,:,slot), view(r3,:,slot), iz, ID,
                          A22_diag, DWS, Th_cache, a12, inv2dz, invfac, is_iz1, polydeg)
        # C_p → SchurD[:,col,iz-1] self correction (node-1 response)
        r2_1 = r2[1,slot];  r3_1 = r3[1,slot]
        SchurD[1, col, iz-1, ID] -= c2p1*r2_1 + c3p1*r3_1
        SchurD[2, col, iz-1, ID] -= c2p2*r2_1 + c3p2*r3_1
        SchurD[3, col, iz-1, ID] -= c2p3*r2_1 + c3p3*r3_1
        # C_self → SchurL[:,col,iz-1]  (full interior(iz) response + Cgrav)
        p1 = zero(FT);  p2 = zero(FT);  p3 = zero(FT);  pg = zero(FT)
        @unroll for j in 1:polydeg
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
=#

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
  D = F.SchurD;  U = F.SchurU;  L = F.SchurL;  Lh = F.SchurLhat;  Di = F.SchurDinv
  nz = F.nz
  NumI = F.NumI
  @inbounds for ID in 1:NumI
    @inbounds begin
      invert_3x3!(Di, D[1,1,1,ID],D[1,2,1,ID],D[1,3,1,ID],
      D[2,1,1,ID],D[2,2,1,ID],D[2,3,1,ID],
      D[3,1,1,ID],D[3,2,1,ID],D[3,3,1,ID], 1, ID)
      @inbounds for i in 2:nz
        # Lhat = L[i-1] · Dinv[i-1]
        @inbounds for c in 1:3
          a1 = Di[1,c,i-1,ID];  a2 = Di[2,c,i-1,ID];  a3 = Di[3,c,i-1,ID]
          Lh[1,c,i-1,ID] = L[1,1,i-1,ID]*a1 + L[1,2,i-1,ID]*a2 + L[1,3,i-1,ID]*a3
          Lh[2,c,i-1,ID] = L[2,1,i-1,ID]*a1 + L[2,2,i-1,ID]*a2 + L[2,3,i-1,ID]*a3
          Lh[3,c,i-1,ID] = L[3,1,i-1,ID]*a1 + L[3,2,i-1,ID]*a2 + L[3,3,i-1,ID]*a3
        end
        # D[i] -= Lhat · U[i-1]
        @inbounds for c in 1:3
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
  nz = F.nz
  NumI = F.NumI
  U = F.SchurU;  Lh = F.SchurLhat;  Di = F.SchurDinv
  @inbounds for ID in 1:NumI
    @inbounds begin
    # forward:  y[i] = rhs[i] - Lhat[i-1]·y[i-1]
      @inbounds for i in 2:nz
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
      @inbounds for i in nz-1:-1:1
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

function ldivVerticalFKernel!(cache_jacobian::JacobianCacheMarcoSplitNS, b, dz)
  (; M, nz, NumI, A12_11, A22_diag, DWS, Th_cache, dpdRhoTh_cache, Gblk,
     rs, SA, invwB) = cache_jacobian
  FT      = eltype(SA)
  invfac  = inv(cache_jacobian.fac)
  cS = P.cS
  invcS   = 1 / cS
  FacGrav = P.Grav
  polydeg = M - 1
  RhoPos  = 1
  wPos = 4
  ThPos = 5
  r2  = zeros(FT, M - 1)
  r3  = zeros(FT, M - 1)
  dpdRhoTh = zeros(FT, M)
  @inbounds for ID in 1:NumI
    @inbounds for iz in 1:nz
      sh     = (iz - 1) * 3
      inv2dz = 2 / dz[iz, ID]
      is_iz1 = (iz == 1)
      rs[sh+1, ID] = b[M, iz, ID, RhoPos]
      rs[sh+2, ID] = b[M, iz, ID, wPos]
      rs[sh+3, ID] = b[M, iz, ID, ThPos]
      @inbounds for i in 1:polydeg
        r2[i] = b[i, iz, ID, ThPos]
        r3[i] = b[i, iz, ID, wPos]
      end
      # r3 -= invfac · G_II · b_ρ   (matvec; pointwise = -fgi·b_ρ)
      @inbounds for l in 1:polydeg
        brl = b[l, iz, ID, RhoPos] * invfac
        @inbounds for i in 1:polydeg
            r3[i] -= Gblk[i, l, iz, ID] * brl
        end
      end
      @inbounds for j in 1:polydeg
        r2j_sc      = r2[j] * A22_diag[j, iz, ID]
        dpdj_inv2dz = dpdRhoTh_cache[j, iz, ID] * inv2dz
        @inbounds for i in 1:polydeg
          r3[i] -= DWS[i, j] * dpdj_inv2dz * r2j_sc
        end
      end
      r3[1] += (is_iz1 ? 2 : 1) * DWS[1, 1] * dpdRhoTh_cache[1, iz, ID] * inv2dz *
                r2[1] * A22_diag[1, iz, ID]
      # A12 cross: r3[i] += invfac·A12_11·G_II[i,1]·A22⁻¹[1]·r2[1]  (column; pointwise row 1)
      a12c = invfac * A12_11[iz, ID] * r2[1] * A22_diag[1, iz, ID]
      @inbounds for i in 1:polydeg
        r3[i] += a12c * Gblk[i, 1, iz, ID]
      end
      ldivFull!(iz, ID, SA, r3, Val(nvnodes(cache_jacobian) - 1))
      @inbounds for j in 1:polydeg
        r3j        = r3[j]
        thj_inv2dz = Th_cache[j, iz, ID] * inv2dz
        @inbounds for i in 1:polydeg
          r2[i] -= DWS[i, j] * thj_inv2dz * r3j
        end
      end
      if !is_iz1
        r2[1] += DWS[1, 1] * Th_cache[1, iz, ID] * inv2dz * r3[1]
      end
      @inbounds for i in 1:polydeg
        r2[i] *= A22_diag[i, iz, ID]
      end
      rs1 = rs[sh+1, ID];  rs2 = rs[sh+2, ID];  rs3 = rs[sh+3, ID]
      @inbounds for j in 1:polydeg
        c1j = inv2dz * DWS[M, j]
        c3j = c1j * Th_cache[j, iz, ID]
        c2j = c1j * dpdRhoTh_cache[j, iz, ID]
        rs1 -= c1j * r3[j]
        rs3 -= c3j * r3[j]
        rs2 -= c2j * r2[j]
      end
      # Cgrav: rs[ρw_M] -= Σ_j G[M,j]·ρ_int[j]  (zero for pointwise)
      cg = zero(FT)
      @inbounds for j in 1:polydeg
        rhoj = b[j, iz, ID, RhoPos]
        if j == 1
          rhoj -= A12_11[iz, ID] * r2[1]
        end
        acc = zero(FT)
        @inbounds for l in 1:polydeg
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

function ldivVerticalBKernel!(cache_jacobian::JacobianCacheMarcoSplitNS, b, dz)
    (; A12_11, A22_diag, DWS, Th_cache, dpdRhoTh_cache, Gblk,
       SA, rs, M, nz, NumI, invwB) = cache_jacobian
  FT      = eltype(SA)
  invfac  = inv(cache_jacobian.fac)
  cS = P.cS
  invcS   = 1 / cS
  polydeg = M - 1
  RhoPos  = 1
  ThPos = 5
  wPos = 4
  r2 = zeros(M - 1)
  r3 = zeros(M - 1)
  @inbounds for ID in 1:NumI
    @inbounds for iz in 1:nz
      sh     = (iz - 1) * 3
      fgi    = P.Grav * invfac
      inv2dz = 2 / dz[iz, ID]
      is_iz1 = (iz == 1)
      b[M, iz, ID, RhoPos] = rs[sh+1, ID]
      b[M, iz, ID, wPos]   = rs[sh+2, ID]
      b[M, iz, ID, ThPos]  = rs[sh+3, ID]
      ThM_self  = Th_cache[M, iz, ID]
      dpdM_self = dpdRhoTh_cache[M, iz, ID]
      rsh2 = rs[sh+2, ID];  rsh3 = rs[sh+3, ID]
      @inbounds for i in 1:polydeg
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
      @inbounds for l in 1:polydeg
        brl = b[l, iz, ID, RhoPos] * invfac
        @inbounds for i in 1:polydeg
          r3[i] -= Gblk[i, l, iz, ID] * brl
        end
      end
      @inbounds for i in 1:polydeg
        r3[i] -= Gblk[i, M, iz, ID] * rhoM
      end
      @inbounds for j in 1:polydeg
        r2j_sc      = r2[j] * A22_diag[j, iz, ID]
        dpdj_inv2dz = dpdRhoTh_cache[j, iz, ID] * inv2dz
        @inbounds for i in 1:polydeg
          r3[i] -= DWS[i, j] * dpdj_inv2dz * r2j_sc
        end
      end
      r3[1] += (is_iz1 ? 2 : 1) * DWS[1, 1] * dpdRhoTh_cache[1, iz, ID] * inv2dz *
      r2[1] * A22_diag[1, iz, ID]
      # A12 cross-term:  r3 += invfac*A12_11 * G[:,1] * (A22 r2)_1   (column over all rows)
      a12c = invfac * A12_11[iz, ID] * r2[1] * A22_diag[1, iz, ID]
      @inbounds for i in 1:polydeg
        r3[i] += a12c * Gblk[i, 1, iz, ID]
      end
      ldivFull!(iz, ID, SA, r3, Val(nvnodes(cache_jacobian) - 1))
      @inbounds for j in 1:polydeg
        r3j        = r3[j]
        thj_inv2dz = Th_cache[j, iz, ID] * inv2dz
        @inbounds for i in 1:polydeg
          r2[i] -= DWS[i, j] * thj_inv2dz * r3j
        end
      end
      if !is_iz1
        r2[1] += DWS[1, 1] * Th_cache[1, iz, ID] * inv2dz * r3[1]
      end
      @inbounds for i in 1:polydeg
        r2[i] *= A22_diag[i, iz, ID]
      end
      b[1, iz, ID, RhoPos] -= A12_11[iz, ID] * r2[1]
      @inbounds for j in 1:polydeg
        r3j_inv2dz = r3[j] * inv2dz
        @inbounds for i in 1:polydeg
          b[i, iz, ID, RhoPos] -= DWS[i, j] * r3j_inv2dz
        end
      end
      if !is_iz1
        b[1, iz, ID, RhoPos] += DWS[1, 1] * inv2dz * r3[1]
      end
      @inbounds for i in 1:polydeg
        b[i, iz, ID, RhoPos] *= invfac
        b[i, iz, ID, ThPos]   = r2[i]
        b[i, iz, ID, wPos]    = r3[i]
      end
    end
  end
end

function Jac!(U,fac,DG,Metric,Phys,Cache,JCache::JacobianCacheMarcoSplitNS,Global,VelForm)
  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU
  @views Geo = Cache.Aux[:,:,:,2]
  precompute_gravity!(Geo,Metric.dz,JCache,NumberThreadGPU)
  Invfac = eltype(U)(1) / fac
  dz = Metric.dz
  @show "Jac KA",Invfac
# FillJacDGVert!(JCache, U, Metric.dz, Invfac,NumberThreadGPU)
  FillJacDGVert!(JCache, U, Metric.dz, Invfac)
  SchurBoundary!(JCache, Metric.dz, NumberThreadGPU::Int64)
  luBandkernel!(JCache)
end

function JacM!(U,fac,DG,Metric,Phys,Cache,JCache::JacobianCacheMarcoSplitNS,Global,VelForm)
  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU
  @views Geo = Cache.Aux[:,:,:,2]
  precompute_gravity!(Geo,Metric.dz,JCache,NumberThreadGPU)
  Invfac = eltype(U)(1) / fac
  dz = Metric.dz
  FillJacDGVert!(JCache, U, Metric.dz, Invfac)
  SchurBoundary!(JCache, Metric.dz, Phys)
  luBandkernel!(JCache)
end

function Solve!(k,v,Jac::JacDGVert,fac,DG::FiniteElements.DGElement,Metric,Global,VelForm)
  
  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU
  @. k = v
  @views TendVCart2VSp!(k,DG,Metric,NumberThreadGPU,VelForm)
  Solve!(Jac,k)
  @views @. k[:,:,:,2:3] *= fac
  @views TendVSp2VCart!(k,DG,Metric,NumberThreadGPU,VelForm)

end

function Solve!(k,v,Jac::JacobianCacheMarcoSplitNS,fac,DG,Metric,Global,VelForm)

  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU
  @. k = v
  @views TendVCart2VSp!(k,DG,Metric,NumberThreadGPU,VelForm)
  solve_jacobian!(k, Jac, Metric)
  @views @. k[:,:,:,2:3] *= fac
  @views TendVSp2VCart!(k,DG,Metric,NumberThreadGPU,VelForm)

end

function solve_jacobian!(b, cache_jacobian::JacobianCacheMarcoSplitNS, Metric)
  ldivVerticalFKernel!(cache_jacobian, b, Metric.dz)
  ldivVerticalSKernel!(cache_jacobian)
  ldivVerticalBKernel!(cache_jacobian, b, Metric.dz)
end

function InitJacDG(DG,nz,Param)
  N = (DG.OrdPolyZ + 1) * nz
  dSdS,dSdM = DScalarDMomAc(nz,DG,Param.cS)
  dMdS,dMdM = DMomDScalarAc(nz,DG,Param.cS)
  return dSdS,dSdM,dMdS,dMdM
end

function JacDGTNeu(U,Aux,DG,fac,dSdS,dSdM,dMdS,dMdM,dGeo,dz,Phys)
  D = DG.DWZ
  FTB = eltype(U)
  N = size(dSdM,1)
  RhoPos = 1
  ThPos = 5
  nz = size(U,2)
  M = size(U,1)
  oneM = ones(M)
  NF = size(dz,3)
  JacLU = Array{SparseArrays.UMFPACK.UmfpackLU}(undef,size(U,3))
    ID = 1
    @views dzCol = dz[:,ID]
    diagdz = spdiagm(2.0 ./ reshape(vec(oneM*dzCol'),N))
    Th = reshape(U[:,:,ID,ThPos]./U[:,:,ID,RhoPos],N)
    dpdRhoTh = reshape( FTB(1) / (FTB(1) - Phys.kappa) * Phys.Rd *
      (Phys.Rd * U[:,:,ID,ThPos] ./ Phys.p0).^(Phys.kappa / (1.0 - Phys.kappa)),N)
    Geo = reshape(Aux[:,:,ID,2],N)
    Jac = [sparse(fac*I,N,N)  -diagdz* dSdS * diagm(dpdRhoTh) -diagdz * dSdM
           spzeros(N,N) sparse(fac*I,N,N) - diagdz*diagm(Th)*dSdS*diagm(dpdRhoTh)  -diagdz*dSdM*diagm(Th)
           dGeo -diagdz*dMdS*diagm(dpdRhoTh) sparse(fac*I,N,N)-diagdz*dMdM]
    JacLU[ID] = lu(Jac)
  return JacLU,Jac
end

function DerivativeGeo(Geo,DG,dz)

  RowInd = Int[]
  ColInd = Int[]
  Val = Float64[]
  DWS = DG.DWZ

  M = size(DWS,1)
  nz = size(dz,1)
  ND = size(dz,2)
  N = M * nz

  Gblk = zeros(M,M,nz,ND)

  @inbounds for ID in 1 : ND
    for iz in 1:nz
      inv2dz = 2 / dz[iz, ID]
      for i in 1:M
        acc = 0.0
        for k in 1:M
          dGeo = Geo[k,iz,ID] - Geo[i,iz,ID]
          Gblk[i,k,iz,ID] = inv2dz * 0.5 * DWS[i,k] * dGeo
          Gblk[i,k,iz,ID] = 0.0
          acc += DWS[i,k] * dGeo
          if ID == 1
#           push!(RowInd,i+(iz-1)*M)
#           push!(ColInd,k+(iz-1)*M)
#           push!(Val,Gblk[i,k,iz,ID])
          end
        end
        if ID == 1
          push!(RowInd,i+(iz-1)*M)
          push!(ColInd,i+(iz-1)*M)
#         push!(Val,inv2dz * 0.5 * acc)
          push!(Val,P.Grav)
        end
#       Gblk[i,i,iz,ID] += inv2dz * 0.5 * acc
        Gblk[i,i,iz,ID] = P.Grav
      end
    end
  end
  DGeo = sparse(RowInd,ColInd,Val,N,N)
  return Gblk, DGeo
end
