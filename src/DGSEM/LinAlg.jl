################################################################################
# Linear algebra for the vertically implicit solvers of the HEVI scheme
#
#   Dense LU, diagonal of U stored inverted  (JacHDG)
#     LUInv!, ldivInv!, ldivInv2!
#   Dense LU, diagonal of U stored as is     (JacTwoSidedSplit, JacOneSided)
#     LUFull!, ldivFull!
#   Block tridiagonal, NB x NB blocks        (JacOneSided, JacTwoSidedSplit)
#     FactorBlockTriKernel!, ldivBlockTriKernel!
#   Scalar tridiagonal from element blocks   (JacHDG)
#     LUTriBlk!, ldivTriBlk!, luTriBlkKernel!, ldivTriBlkKernel!
#
# All factorizations are without pivoting. Always use the LU and ldiv routines
# of the same family together.
################################################################################

################################################################################
# Dense LU for small, register-resident matrices (MMatrix / SMatrix).
# L (unit diagonal) is stored below the diagonal, U on and above it, and the
# diagonal of U is stored inverted (1/u_kk). ldivInv! multiplies by it.
################################################################################

@inline function LUInv!(A, ::Val{n}) where {n}
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

@inline function ldivInv!(A, b, ::Val{n}) where {n}
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

# Solve two right-hand sides with the same factors in a single sweep:
# the factors are read once and the two independent chains give more ILP.
@inline function ldivInv2!(A, b1, b2, ::Val{n}) where {n}
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
# Dense LU with the diagonal of U stored as is, for a local matrix A[i,j] or
# for the element matrices A[i,j,iz,ID] of a global array.
# JacTwoSidedSplit also reads these factors directly.
################################################################################

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
# Backward loop
  @inbounds for k = n : -1 : 1
    b[k] /= A[k,k,iz,ID]
    @inbounds for i = 1 : k - 1
      b[i] -= A[i,k,iz,ID] * b[k]
    end
  end
end

################################################################################
# Tridiagonal Schur complement for the vertical HDG solver (JacHDG)
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
# Block Thomas algorithm for block tridiagonal systems with NB x NB blocks,
# one system per column ID. The diagonal block of row iz is Dself[iz]
# (plus Dprev[iz+1] if HasPrev, contributions computed by the element above),
# the sub- and super-diagonal blocks are Lo[iz] and Up[iz]. The factorization
# stores the inverted pivot blocks Binv and Cp = Binv * Up. The right hand side
# is rsA[iz] (plus rsB[iz+1] if HasPrev).
################################################################################

@inline LoadBlock(A,ID,iz,::Val{NB}) where {NB} =
  SMatrix{NB,NB,eltype(A)}(ntuple(k -> A[ID,mod1(k,NB),cld(k,NB),iz], Val(NB*NB)))

@kernel inbounds = true function FactorBlockTriKernel!(Binv,Cp,@Const(Dself),@Const(Dprev),
  @Const(Lo),@Const(Up),::Val{NB},::Val{nz},::Val{HasPrev}) where {NB,nz,HasPrev}

  ID, = @index(Global, NTuple)
  ND = @uniform @ndrange()[1]

  if ID <= ND
    FT = eltype(Binv)
    CpOld = zero(SMatrix{NB,NB,FT})
    for iz = 1 : nz
      D = LoadBlock(Dself,ID,iz,Val(NB))
      if HasPrev && iz < nz
        D += LoadBlock(Dprev,ID,iz+1,Val(NB))
      end
      if iz > 1
        D -= LoadBlock(Lo,ID,iz,Val(NB)) * CpOld
      end
      Bi = inv(D)
      if iz < nz
        CpOld = Bi * LoadBlock(Up,ID,iz,Val(NB))
      end
      @unroll for j = 1 : NB
        @unroll for i = 1 : NB
          Binv[ID,i,j,iz] = Bi[i,j]
          Cp[ID,i,j,iz] = CpOld[i,j]
        end
      end
    end
  end
end

@kernel inbounds = true function ldivBlockTriKernel!(rs,@Const(Binv),@Const(Cp),@Const(Lo),
  @Const(rsA),@Const(rsB),::Val{NB},::Val{nz},::Val{HasPrev}) where {NB,nz,HasPrev}

  ID, = @index(Global, NTuple)
  ND = @uniform @ndrange()[1]

  if ID <= ND
    FT = eltype(rs)
    y = zero(SVector{NB,FT})
    for iz = 1 : nz
      d = SVector{NB,FT}(ntuple(i -> rsA[ID,i,iz], Val(NB)))
      if HasPrev && iz < nz
        d += SVector{NB,FT}(ntuple(i -> rsB[ID,i,iz+1], Val(NB)))
      end
      if iz > 1
        d -= LoadBlock(Lo,ID,iz,Val(NB)) * y
      end
      y = LoadBlock(Binv,ID,iz,Val(NB)) * d
      @unroll for i = 1 : NB
        rs[ID,i,iz] = y[i]
      end
    end
    for iz = nz - 1 : -1 : 1
      y = SVector{NB,FT}(ntuple(i -> rs[ID,i,iz], Val(NB))) - LoadBlock(Cp,ID,iz,Val(NB)) * y
      @unroll for i = 1 : NB
        rs[ID,i,iz] = y[i]
      end
    end
  end
end
