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

@kernel inbounds = true function luBandKernel!(A,::Val{n}) where {n}
  ID, = @index(Global, NTuple)

  DoF = @uniform @ndrange()[1]

  if ID <= DoF
    luBand!(ID,A,3,3,Val(n))
  end
end

@kernel inbounds = true function ldivVerticalSKernel!(A,rs,::Val{n}) where {n}

  ID, = @index(Global, NTuple)

  DoF = @uniform @ndrange()[1]

  if ID <= DoF
    ldivBand!(ID,A,rs,3,3,Val(n))
  end
end
