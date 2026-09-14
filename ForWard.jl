@kernel inbounds = true function ldivHDGVerticalFKernel!(@Const(A13),@Const(A23),@Const(A31),@Const(A32),
  @Const(C2), @Const(C3),@Const(SA),@Const(b),rs,invfac, ::Val{M}) where {M}

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
    @unroll for i = 1 : M
      for j = 1 : M
        r3[i] -= invfac * (A31[i,j,iz,ID] * r1[j] + A32[i,j,iz,ID] * r2[j])
      end
    end  

    ldivFull!(iz,ID,SA,r3,Val(M))

    @unroll for i = 1 : M 
      @unroll for j = 1 : M
        r2[i] -= A23[i,j,iz,ID] * r3[j]
      end
      r2[i] *= invfac
    end

    if iz == 1
      i = 2 * iz - 1
      c2 = MVector{M}([C2[2,1,iz,ID];zeros(eltype(SA),M-1)])
      c3 = MVector{M}([C3[2,1,iz,ID];zeros(eltype(SA),M-1)])
      t = (r2'*c2 + r3'*c3)
      @atomic rs[i,ID] -= t

      i = 2 * iz
      c2 = MVector{M}([zeros(eltype(SA),M-1);C2[1,2,iz,ID]])
      c3 = MVector{M}([zeros(eltype(SA),M-1);C3[1,2,iz,ID]])
      t = (r2'*c2 + r3'*c3)
      @atomic rs[i,ID] -= t

      i = 2 * iz + 1
      c2 = MVector{M}([zeros(eltype(SA),M-1);C2[2,2,iz,ID]])
      c3 = MVector{M}([zeros(eltype(SA),M-1);C3[2,2,iz,ID]])
      t = (r2'*c2 + r3'*c3)
      @atomic rs[i,ID] -= t
    end  
    if iz > 1 && iz < nz
      i = 2 * iz - 2
      c2 = MVector{M}([C2[1,1,iz,ID];zeros(eltype(SA),M-1)])
      c3 = MVector{M}([C3[1,1,iz,ID];zeros(eltype(SA),M-1)])
      t = (r2'*c2 + r3'*c3)
      @atomic rs[i,ID] -= t

      i = 2 * iz - 1
      c2 = MVector{M}([C2[2,1,iz,ID];zeros(eltype(SA),M-1)])
      c3 = MVector{M}([C3[2,1,iz,ID];zeros(eltype(SA),M-1)])
      t = (r2'*c2 + r3'*c3)
      @atomic rs[i,ID] -= t

      i = 2 * iz
      c2 = MVector{M}([zeros(eltype(SA),M-1);C2[1,2,iz,ID]])
      c3 = MVector{M}([zeros(eltype(SA),M-1);C3[1,2,iz,ID]])
      t = (r2'*c2 + r3'*c3)
      @atomic rs[i,ID] -= t

      i = 2 * iz + 1
      c2 = MVector{M}([zeros(eltype(SA),M-1);C2[2,2,iz,ID]])
      c3 = MVector{M}([zeros(eltype(SA),M-1);C3[2,2,iz,ID]])
      t = (r2'*c2 + r3'*c3)
      @atomic rs[i,ID] -= t
    end  
    if iz == nz
      i = 2 * iz - 2
      c2 = MVector{M}([C2[1,1,iz,ID];zeros(eltype(SA),M-1)])
      c3 = MVector{M}([C3[1,1,iz,ID];zeros(eltype(SA),M-1)])
      t = (r2'*c2 + r3'*c3)
      @atomic rs[i,ID] -= t

      i = 2 * iz - 1
      c2 = MVector{M}([C2[2,1,iz,ID];zeros(eltype(SA),M-1)])
      c3 = MVector{M}([C3[2,1,iz,ID];zeros(eltype(SA),M-1)])
      t = (r2'*c2 + r3'*c3)
      @atomic rs[i,ID] -= t

      i = 2 * iz
      c2 = MVector{M}([zeros(eltype(SA),M-1);C2[2,2,iz,ID]])
      c3 = MVector{M}([zeros(eltype(SA),M-1);C3[2,2,iz,ID]])
      t = (r2'*c2 + r3'*c3)
      @atomic rs[i,ID] -= t
    end    
  end
end
