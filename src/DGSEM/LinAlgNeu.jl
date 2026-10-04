################################################################################
# Dense LU for small, register-resident matrices (MMatrix / SMatrix)
#
# CONVENTION (single one, for both Val and Integer entry points):
#   LUFull!  : in-place LU without pivoting. L (unit diagonal) is stored below
#              the diagonal, U on and above it, and the DIAGONAL OF U IS STORED
#              INVERTED (1/u_kk).
#   ldivFull!: solves with these factors; it MULTIPLIES by the stored inverse
#              diagonal (no division).
# Always use the pair LUFull! / ldivFull! (or ldivFull2!) together.
################################################################################

@inline function _LUFull_inv!(A, n)
  T = eltype(A)
  @unroll for j = 1 : n - 1
    invAjj = one(T) / A[j,j]
    @unroll for i = j + 1 : n
      A[i,j] *= invAjj
      @unroll for k = j + 1 : n
        A[i,k] = muladd(-A[i,j], A[j,k], A[i,k])
      end
    end
    A[j,j] = invAjj
  end
  A[n,n] = one(T) / A[n,n]
end

@inline LUFull!(A, ::Val{n}) where {n} = _LUFull_inv!(A, n)
@inline LUFull!(A, n::Integer)         = _LUFull_inv!(A, n)

@inline function _ldivFull_inv!(A, b, n)
  # Forward loop (L has unit diagonal)
  @unroll for k = 1 : n - 1
    @unroll for i = k + 1 : n
      b[i] = muladd(-A[i,k], b[k], b[i])
    end
  end
  # Backward loop (diagonal of U is stored inverted)
  @unroll for k = n : -1 : 1
    b[k] *= A[k,k]
    @unroll for i = 1 : k - 1
      b[i] = muladd(-A[i,k], b[k], b[i])
    end
  end
end

@inline ldivFull!(A, b, ::Val{n}) where {n} = _ldivFull_inv!(A, b, n)
@inline ldivFull!(A, b, n::Integer)         = _ldivFull_inv!(A, b, n)

# Solve two right-hand sides with the same factors in a single sweep:
# the factors are read once and the two independent chains give more ILP.
@inline function ldivFull2!(A, b1, b2, ::Val{n}) where {n}
  @unroll for k = 1 : n - 1
    @unroll for i = k + 1 : n
      a = A[i,k]
      b1[i] = muladd(-a, b1[k], b1[i])
      b2[i] = muladd(-a, b2[k], b2[i])
    end
  end
  @unroll for k = n : -1 : 1
    d = A[k,k]
    b1[k] *= d
    b2[k] *= d
    @unroll for i = 1 : k - 1
      a = A[i,k]
      b1[i] = muladd(-a, b1[k], b1[i])
      b2[i] = muladd(-a, b2[k], b2[i])
    end
  end
end

################################################################################
# Tridiagonal Schur complement for the vertical HDG solver (JacHDGVert)
#
# Memory layout: the column index ID is the FASTEST index, so that consecutive
# threads (one thread per column) access consecutive addresses (coalesced).
#
#   Blk       (NumI, 4, nz)   per-element 2x2 contributions to the Schur matrix
#                             (values that are ADDED):
#                               1: diag  of face iz-1   (needs iz > 1)
#                               2: sub   entry (iz, iz-1)  (needs 1 < iz < nz)
#                               3: super entry (iz-1, iz)  (needs 1 < iz < nz)
#                               4: diag  of face iz     (needs iz < nz)
#   A         (NumI, 3, N)    factored tridiagonal, N = nz-1 faces:
#                               A[ID,3,i]   = l_i      (multiplier, sub-diagonal)
#                               A[ID,2,i]   = 1/d_i    (inverse pivot)
#                               A[ID,1,i+1]= u_i/d_i  (super-diagonal * inv pivot)
#   RsC       (NumI, 2, nz)   per-element RHS contributions:
#                               1: to face iz-1,  2: to face iz
#   rs        (NumI, N)       right-hand side / solution on the faces
################################################################################

# Assemble the band from the element blocks and factor it (Thomas algorithm).
# The running pivot is carried in a register (no read-after-write on memory).
@inline function LUTriBlk!(ID, A, Blk, cS2, N)
  dinv = inv(cS2 + Blk[ID,4,1] + Blk[ID,1,2])
  A[ID,2,1] = dinv
  for i in 1:N-1
    s = Blk[ID,2,i+1]                 # entry (i+1,i)
    u = Blk[ID,3,i+1]                 # entry (i,i+1)
    l = s * dinv
    A[ID,3,i]   = l
    A[ID,1,i+1] = u * dinv
    dinv = inv(cS2 + Blk[ID,4,i+1] + Blk[ID,1,i+2] - l * u)
    A[ID,2,i+1] = dinv
  end
end

# Gather the face RHS from the element contributions and solve L U x = b.
# Forward/backward sweeps carry the running value in a register; the backward
# step is a single muladd on the dependency chain.
@inline function ldivTriBlk!(ID, A, RsC, rs, N)
  # Forward: L y = b
  y = RsC[ID,2,1] + RsC[ID,1,2]
  rs[ID,1] = y
  for i in 2:N
    bi = RsC[ID,2,i] + RsC[ID,1,i+1]
    y  = muladd(-A[ID,3,i-1], y, bi)
    rs[ID,i] = y
  end
  # Backward: U x = y
  x = y * A[ID,2,N]
  rs[ID,N] = x
  for i in N-1:-1:1
    x = muladd(-A[ID,1,i+1], x, rs[ID,i] * A[ID,2,i])
    rs[ID,i] = x
  end
end

@kernel inbounds = true function luTriBlkKernel!(A, @Const(Blk), cS2, N)
  ID, = @index(Global, NTuple)
  ND = @uniform @ndrange()[1]
  if ID <= ND
    LUTriBlk!(ID, A, Blk, cS2, N)
  end
end

@kernel inbounds = true function ldivTriBlkKernel!(@Const(A), @Const(RsC), rs, N)
  ID, = @index(Global, NTuple)
  ND = @uniform @ndrange()[1]
  if ID <= ND
    ldivTriBlk!(ID, A, RsC, rs, N)
  end
end

################################################################################
# LEGACY routines (kept unchanged so that other code keeps working).
# They use their own layouts / conventions. Delete if nothing else calls them.
################################################################################

# --- dense LU on a 4D global array; diagonal NOT inverted (pair with the
# --- 4-argument ldivFull! below) ---------------------------------------------
@inline function LUFull!(iz,ID,A)

  n = size(A,1)
  for j = 1 : n - 1
    invAjj = eltype(A)(1) / A[j,j,iz,ID]
    for i = j + 1 : n 
      A[i,j,iz,ID] *= invAjj
      for k = j + 1 : n
        A[i,k,iz,ID] -= A[i,j,iz,ID] * A[j,k,iz,ID]
      end  
    end  
  end  
end

@inline function ldivFull!(iz,ID,A,b)

  n = size(A,1)
# Forward loop
  @inbounds for k = 1 : n - 1
    for i = k + 1 : n
      b[i] -= A[i,k,iz,ID] * b[k]
    end
  end
# Backward loop
  for k = n : -1 : 1
    b[k] /= A[k,k,iz,ID]
    for i = 1 : k - 1
      b[i] -= A[i,k,iz,ID] * b[k]
    end
  end
end

@inline function ldivBand1!(ID,A,b,::Val{l},::Val{u},::Val{n}) where {l,u,n}

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

# --- tridiagonal with layout (3, N, ND) -------------------------------------
@inline function LUTri!(ID,A,N) 
    # Handle the very first diagonal element
    A[2, 1, ID] = eltype(A)(1) / A[2, 1, ID]

    @unroll for i in 1:(N-1)
        # 1. Compute multiplier for L and store it in Row 3 (dl)
        # We multiply by the pre-computed inverse diagonal of row 'i'
        A[3, i, ID] = A[3, i, ID] * A[2, i, ID]

        # 2. Update the next diagonal element of U in Row 2 (d)
        d_next = A[2, i+1, ID] - A[3, i, ID] * A[1, i+1, ID]

        # 3. Store the reciprocal of the next diagonal element
        A[2, i+1, ID] = eltype(A)(1) / d_next
    end
end

@inline function ldivTri!(ID,A, b, N)
    # --- Step 1: Forward Substitution (Solve L * y = b) ---
    # We overwrite b with y in-place.
    @unroll for i in 2:N
        b[i, ID] = b[i, ID] - A[3, i-1, ID] * b[i-1, ID]
    end

    # --- Step 2: Backward Substitution (Solve U * x = y) ---
    # We overwrite y (currently in b) with x in-place.
    # Uses the pre-computed inverse diagonal (Row 2) to avoid hardware division.
    b[N, ID] = b[N, ID] * A[2, N, ID]
    @unroll for i in (N-1):-1:1
        b[i, ID] = (b[i, ID] - A[1, i+1, ID] * b[i+1, ID]) * A[2, i, ID]
    end
end

@kernel inbounds = true function luTriKernel!(A,::Val{n}) where {n}
  ID, = @index(Global, NTuple)

  ND = @uniform @ndrange()[1]

  if ID <= ND
    LUTri!(ID,A,n)
  end
end

@kernel inbounds = true function ldivVerticalTriKernel!(@Const(A),rs,::Val{n}) where {n}

  ID, = @index(Global, NTuple)

  ND = @uniform @ndrange()[1]

  if ID <= ND
    ldivTri!(ID,A,rs,n)
  end
end

# --- general band solvers ----------------------------------------------------
@inline function luBand!(ID,A,::Val{kl},::Val{ku},::Val{n}) where {kl,ku,n}
  diag_row = ku + 1 # The row index where the main diagonal lives
  @inbounds for k in 1:n
    # Pivot value is at A[k,k] -> Ab[diag_row, k]
    pivot = A[diag_row,k,ID]

    # 1. Compute multipliers for rows below k within lower bandwidth
    @inbounds for i in (k+1):min(k+kl, n)
      # A[i,k] maps to Ab[diag_row + i - k, k]
      A[diag_row + i - k,k,ID] /= pivot
    end

    # 2. Update remaining submatrix elements within the bands
    @inbounds for j in (k+1):min(k+ku, n)
      # A[k,j] maps to Ab[diag_row + k - j, j]
      u_kj = A[diag_row + k - j,j,ID]
            
      @inbounds for i in (k+1):min(k+kl, n)
        # A[i,j] maps to Ab[diag_row + i - j, j]
        A[diag_row + i - j,j,ID] -= A[diag_row + i - k,k,ID] * u_kj
      end
    end
  end
end

@inline function ldivBand!(ID,A,b,::Val{kl},::Val{ku},::Val{n}) where {kl,ku,n}
  diag_row = ku + 1

  # 1. Forward Substitution (L * y = b)
  for k in 1:n
    for i in (k+1):min(k+kl, n)
      # L_ik is stored at Ab_lu[diag_row + i - k, k]
      b[i,ID] -= A[diag_row + i - k,k,ID] * b[k,ID]
    end
  end

  # 2. Backward Substitution (U * x = y)
  for k in n:-1:1
    # U_kk is stored at Ab_lu[diag_row, k]
    b[k,ID] /= A[diag_row,k,ID]

    for i in max(1, k-ku):(k-1)
      # U_ik is stored at Ab_lu[diag_row + i - k, k]
      b[i,ID] -= A[diag_row + i - k,k,ID] * b[k,ID]
    end
  end
end

@inline function luBand1!(ID,A,::Val{l},::Val{u},::Val{n}) where {l,u,n}
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

@kernel inbounds = true function luBandKernel!(A,::Val{l},::Val{u},::Val{n}) where {l,u,n}
  ID, = @index(Global, NTuple)

  DoF = @uniform @ndrange()[1]

  if ID <= DoF
    luBand!(ID,A,Val(l),Val(u),Val(n))
  end
end

@kernel inbounds = true function ldivVerticalSKernel!(A,rs,::Val{l},::Val{u},::Val{n}) where {l,u,n}

  ID, = @index(Global, NTuple)

  DoF = @uniform @ndrange()[1]

  if ID <= DoF
    ldivBand!(ID,A,rs,Val(l),Val(u),Val(n))
  end
end
