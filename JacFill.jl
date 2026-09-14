@kernel inbounds = true function FillJacHDGVertKernel!(A13,A23,@Const(A31),A32,
  B1,B2,B3,C2,C3,SA,SchurBand,@Const(U),@Const(dz),
  @Const(DW),@Const(w),fac,cS,Phys, ::Val{M}) where {M}

  iz,ID = @index(Global, NTuple)

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
    sh = (iz - 1) * 4 + 1
    inv2dz = eltype(SA)(2) / dz[iz,ID]
    invdz = eltype(SA)(1) / dz[iz,ID]
    invwBdz = invdz / wB
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
      B1[1,iz,ID] = eltype(SA)(0)
      B1[2,iz,ID] = invdz 
      B2[1,iz,ID] = eltype(SA)(0)
      B2[2,iz,ID] = eltype(SA)(0.5) * (Th[M] + Thp) * invdz 
      B3[1,iz,ID] = inv2dz
      B3[2,iz,ID] = -invdz 
      C2[1,1,iz,ID] = eltype(SA)(0)
      C2[1,2,iz,ID] = -dpdRhoTh[M] * invcSwB 
      C2[2,1,iz,ID] = dpdRhoTh[1] * invwB
      C2[2,2,iz,ID] = dpdRhoTh[M] * invwB 
      C3[1,1,iz,ID] = eltype(SA)(0)
      C3[1,2,iz,ID] = -invwB
      C3[2,1,iz,ID] = -invwB * cS
      C3[2,2,iz,ID] = invwB * cS
    elseif iz == nz
      Thm = U[M,iz - 1,ID,ThPos] / U[M,iz - 1,ID,RhoPos]  
      B1[1,iz,ID] = -invdz  
      B1[2,iz,ID] = eltype(SA)(0)
      B2[1,iz,ID] = -eltype(SA)(0.5) * (Th[1] + Thm) * invdz   
      B2[2,iz,ID] = eltype(SA)(0)
      B3[1,iz,ID] = invdz  
      B3[2,iz,ID] = -inv2dz
      C2[1,1,iz,ID] = dpdRhoTh[1] * invcSwB
      C2[1,2,iz,ID] = eltype(SA)(0)
      C2[2,1,iz,ID] = dpdRhoTh[1] * invwB
      C2[2,2,iz,ID] = dpdRhoTh[M] * invwB
      C3[1,1,iz,ID] = -invwB
      C3[1,2,iz,ID] = eltype(SA)(0)
      C3[2,1,iz,ID] = -invwB * cS
      C3[2,2,iz,ID] = invwB * cS
    else
      Thm = U[M,iz - 1,ID,ThPos] / U[M,iz - 1,ID,RhoPos]  
      Thp = U[1,iz + 1,ID,ThPos] / U[1,iz + 1,ID,RhoPos]  
      B1[1,iz,ID] = -invdz   
      B1[2,iz,ID] = invdz 
      B2[1,iz,ID] = -eltype(SA)(0.5) * (Th[1] + Thm) * invdz   
      B2[2,iz,ID] = eltype(SA)(0.5) * (Th[M] + Thp) * invdz 
      B3[1,iz,ID] = invdz   
      B3[2,iz,ID] = -invdz 
      C2[1,1,iz,ID] = dpdRhoTh[1] * invcSwB
      C2[1,2,iz,ID] = -dpdRhoTh[M] * invcSwB 
      C2[2,1,iz,ID] = dpdRhoTh[1] * invwB
      C2[2,2,iz,ID] = dpdRhoTh[M] * invwB 
      C3[1,1,iz,ID] = -invwB
      C3[1,2,iz,ID] = -invwB
      C3[2,1,iz,ID] = -invwB * cS
      C3[2,2,iz,ID] = invwB * cS
    end


    @unroll for i = 1 : M 
      @unroll for j = 1 : M 
        SAL[i,j] = eltype(SA)(0)
        @unroll for k = 1 : M 
          SAL[i,j] -= A32[i,k,iz,ID] * A23[k,j,iz,ID] + A31[i,k,iz,ID] * A13[k,j,iz,ID]
        end
        SAL[i,j] *= fac 
      end
      SAL[i,i] += invfac
    end

    LUFull!(SAL,Val(M))
    @unroll for i = 1 : M
      @unroll for j = 1 : M
        SA[i,j,iz,ID] = SAL[i,j]
      end
    end  


#   iz = 1
#   Spalte iz             : r1[M]=B1[2,iz]      2*iz
#   Spalte iz             : r2[M]=B2[2,iz]      2*iz

#   Spalte nz - 1 + iz    : r3[1]=B3[1,iz]      2*(iz-1)+1

#   Spalte nz - 1 + iz +1 : r3[M]=B3[2,iz+1]    2*iz+1

#   iz = 2
#   Spalte iz - 1         : r1[1]=B1[1,iz-1]    2*(iz-1)
#   Spalte iz - 1         : r2[1]=B2[1,iz-1     2*(iz-1) 

#   Spalte iz             : r1[M]=B1[2,iz]      2*iz
#   Spalte iz             : r2[M]=B2[2,iz]      2*iz

#   Spalte nz - 1 + iz    : r3[1]=B3[1,iz]      2*(iz-1)+1

#   Spalte nz - 1 + iz +1 : r3[M]=B3[2,iz+1]    2*iz+1
  end

  if iz == 1
    j = 2 * iz - 1
    r3[1] = B3[1,iz,ID]
    @unroll for i = 2 : M
      r3[i] = eltype(SA)(0)
    end  
    ldivFull!(iz,ID,SA,r3,Val(M))
    r21 = -A23[1,1,iz,ID] * r3[1]
    r2M = -A23[M,1,iz,ID] * r3[1]
    @unroll for j = 2 : M
      r21 -= A23[1,j,iz,ID] * r3[j]
      r2M -= A23[M,j,iz,ID] * r3[j]
    end  
    r21 *= fac
    r2M *= fac

    i = 2 * iz - 1
    t = r21 * C2[2,1,iz,ID] + r3[1] * C3[2,1,iz,ID]
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t

    i = 2 * iz
    t = r2M * C2[1,2,iz,ID] + r3[M] * C3[1,2,iz,ID]
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t

    i = 2 * iz + 1
    t = r2M * C2[2,2,iz,ID] + r3[M] * C3[2,2,iz,ID]
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t


    j = 2 * iz 
    r1M = B1[2,iz,ID] 
    r2M = B2[2,iz,ID]
    @unroll for i = 1 : M
      r3[i] = -(A31[i,M,iz,ID] * r1M + A32[i,M,iz,ID] * r2M) * fac
    end
    ldivFull!(iz,ID,SA,r3,Val(M))
    r21 = -A23[1,1,iz,ID] * r3[1]
    r2M = r2M - A23[M,1,iz,ID] * r3[1]
    @unroll for j = 1 : M
      r21 -= A23[1,j,iz,ID] * r3[j]
      r2M -= A23[M,j,iz,ID] * r3[j]
    end  
    r21 *= fac
    r2M *= fac

    i = 2 * iz - 1
    t = r21 * C2[2,1,iz,ID] + r3[1] * C3[2,1,iz,ID]
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t

    i = 2 * iz
    t = r2M * C2[1,2,iz,ID] + r3[M] * C3[1,2,iz,ID]
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t

    i = 2 * iz + 1
    t = r2M * C2[2,2,iz,ID] + r3[M] * C3[2,2,iz,ID]
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t

    j = 2 * iz + 1
    r3[M] = B3[2,iz,ID]
    @unroll for i = 1 : M - 1
      r3[i] = eltype(SA)(0)
    end  
    ldivFull!(iz,ID,SA,r3,Val(M))
    r21 = -A23[1,1,iz,ID] * r3[1]
    r2M = -A23[M,1,iz,ID] * r3[1]
    @unroll for j = 2 : M
      r21 -= A23[1,j,iz,ID] * r3[j]
      r2M -= A23[M,j,iz,ID] * r3[j]
    end
    r21 *= fac
    r2M *= fac

    i = 2 * iz - 1
    t = r21 * C2[2,1,iz,ID] + r3[1] * C3[2,1,iz,ID]
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t

    i = 2 * iz
    t = r2M * C2[1,2,iz,ID] + r3[M] * C3[1,2,iz,ID]
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t

    i = 2 * iz + 1
    t = r2M * C2[2,2,iz,ID] + r3[M] * C3[2,2,iz,ID]
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t
  end    

  if iz > 1 && iz < nz
    j = 2 * iz - 2  
    r1 = MVector{M}([B1[1,iz,ID];zeros(eltype(SA),M-1)])
    r2 = MVector{M}([B2[1,iz,ID];zeros(eltype(SA),M-1)])
    r3 = MVector{M}(zeros(eltype(SA),M))

    r3 -= (A31[:,:,iz,ID] * r1 + A32[:,:,iz,ID] * r2) * fac
    ldivFull!(iz,ID,SA,r3,Val(M))
    r2 = (r2 - A23[:,:,iz,ID] * r3) * fac

    i = 2 * iz - 2 
    c2 = MVector{M}([C2[1,1,iz,ID];zeros(eltype(SA),M-1)])
    c3 = MVector{M}([C3[1,1,iz,ID];zeros(eltype(SA),M-1)])
    t = (r2'*c2 + r3'*c3)
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t

    i = 2 * iz - 1  
    c2 = MVector{M}([C2[2,1,iz,ID];zeros(eltype(SA),M-1)])
    c3 = MVector{M}([C3[2,1,iz,ID];zeros(eltype(SA),M-1)])
    t = (r2'*c2 + r3'*c3)
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t
   
    i = 2 * iz
    c2 = MVector{M}([zeros(eltype(SA),M-1);C2[1,2,iz,ID]])
    c3 = MVector{M}([zeros(eltype(SA),M-1);C3[1,2,iz,ID]])
    t = (r2'*c2 + r3'*c3)
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t

    i = 2 * iz + 1 
    c2 = MVector{M}([zeros(eltype(SA),M-1);C2[2,2,iz,ID]])
    c3 = MVector{M}([zeros(eltype(SA),M-1);C3[2,2,iz,ID]])
    t = (r2'*c2 + r3'*c3)
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t

    j = 2 * iz - 1
    r2 = MVector{M}(zeros(eltype(SA),M))
    r3 = MVector{M}([B3[1,iz,ID];zeros(eltype(SA),M-1)])
    ldivFull!(iz,ID,SA,r3,Val(M))
    r2 = (r2 -A23[:,:,iz,ID] * r3) * fac

    i = 2 * iz - 2
    c2 = MVector{M}([C2[1,1,iz,ID];zeros(eltype(SA),M-1)])
    c3 = MVector{M}([C3[1,1,iz,ID];zeros(eltype(SA),M-1)])
    t = (r2'*c2 + r3'*c3)
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t

    i = 2 * iz - 1
    c2 = MVector{M}([C2[2,1,iz,ID];zeros(eltype(SA),M-1)])
    c3 = MVector{M}([C3[2,1,iz,ID];zeros(eltype(SA),M-1)])
    t = (r2'*c2 + r3'*c3)
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t

    i = 2 * iz
    c2 = MVector{M}([zeros(eltype(SA),M-1);C2[1,2,iz,ID]])
    c3 = MVector{M}([zeros(eltype(SA),M-1);C3[1,2,iz,ID]])
    t = (r2'*c2 + r3'*c3)
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t

    i = 2 * iz + 1
    c2 = MVector{M}([zeros(eltype(SA),M-1);C2[2,2,iz,ID]])
    c3 = MVector{M}([zeros(eltype(SA),M-1);C3[2,2,iz,ID]])
    t = (r2'*c2 + r3'*c3)
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t


    j = 2 * iz 
    r1 = MVector{M}([zeros(eltype(SA),M-1);B1[2,iz,ID]])
    r2 = MVector{M}([zeros(eltype(SA),M-1);B2[2,iz,ID]])
    r3 = MVector{M}(zeros(eltype(SA),M))
    r3 -= (A31[:,:,iz,ID] * r1 + A32[:,:,iz,ID] * r2) * fac
    ldivFull!(iz,ID,SA,r3,Val(M))
    r2 = (r2 -A23[:,:,iz,ID] * r3) * fac

    i = 2 * iz - 2
    c2 = MVector{M}([C2[1,1,iz,ID];zeros(eltype(SA),M-1)])
    c3 = MVector{M}([C3[1,1,iz,ID];zeros(eltype(SA),M-1)])
    t = (r2'*c2 + r3'*c3)
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t

    i = 2 * iz - 1
    c2 = MVector{M}([C2[2,1,iz,ID];zeros(eltype(SA),M-1)])
    c3 = MVector{M}([C3[2,1,iz,ID];zeros(eltype(SA),M-1)])
    t = (r2'*c2 + r3'*c3)
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t

    i = 2 * iz
    c2 = MVector{M}([zeros(eltype(SA),M-1);C2[1,2,iz,ID]])
    c3 = MVector{M}([zeros(eltype(SA),M-1);C3[1,2,iz,ID]])
    t = (r2'*c2 + r3'*c3)
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t

    i = 2 * iz + 1
    c2 = MVector{M}([zeros(eltype(SA),M-1);C2[2,2,iz,ID]])
    c3 = MVector{M}([zeros(eltype(SA),M-1);C3[2,2,iz,ID]])
    t = (r2'*c2 + r3'*c3)
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t

    j = 2 * iz + 1
    r2 = MVector{M}(zeros(eltype(SA),M))
    r3 = MVector{M}([zeros(eltype(SA),M-1);B3[2,iz,ID]])
    ldivFull!(iz,ID,SA,r3,Val(M))
    r2 -= (A23[:,:,iz,ID] * r3) * fac

    i = 2 * iz - 2
    c2 = MVector{M}([C2[1,1,iz,ID];zeros(eltype(SA),M-1)])
    c3 = MVector{M}([C3[1,1,iz,ID];zeros(eltype(SA),M-1)])
    t = (r2'*c2 + r3'*c3)
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t

    i = 2 * iz - 1
    c2 = MVector{M}([C2[2,1,iz,ID];zeros(eltype(SA),M-1)])
    c3 = MVector{M}([C3[2,1,iz,ID];zeros(eltype(SA),M-1)])
    t = (r2'*c2 + r3'*c3)
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t

    i = 2 * iz
    c2 = MVector{M}([zeros(eltype(SA),M-1);C2[1,2,iz,ID]])
    c3 = MVector{M}([zeros(eltype(SA),M-1);C3[1,2,iz,ID]])
    t = (r2'*c2 + r3'*c3)
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t

    i = 2 * iz + 1
    c2 = MVector{M}([zeros(eltype(SA),M-1);C2[2,2,iz,ID]])
    c3 = MVector{M}([zeros(eltype(SA),M-1);C3[2,2,iz,ID]])
    t = (r2'*c2 + r3'*c3)
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t
  end    

  if iz == nz
    j = 2 * iz - 2  
    r1 = MVector{M}([B1[1,iz,ID];zeros(eltype(SA),M-1)])
    r2 = MVector{M}([B2[1,iz,ID];zeros(eltype(SA),M-1)])
    r3 = MVector{M}(zeros(eltype(SA),M))
        
    r3 -= (A31[:,:,iz,ID] * r1 + A32[:,:,iz,ID] * r2) * fac
    ldivFull!(iz,ID,SA,r3,Val(M))
    r2 = (r2 - A23[:,:,iz,ID] * r3) * fac

    i = 2 * iz - 2 
    c2 = MVector{M}([C2[1,1,iz,ID];zeros(eltype(SA),M-1)])
    c3 = MVector{M}([C3[1,1,iz,ID];zeros(eltype(SA),M-1)])
    t = (r2'*c2 + r3'*c3)
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t

    i = 2 * iz - 1  
    c2 = MVector{M}([C2[2,1,iz,ID];zeros(eltype(SA),M-1)])
    c3 = MVector{M}([C3[2,1,iz,ID];zeros(eltype(SA),M-1)])
    t = (r2'*c2 + r3'*c3)
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t
   
    i = 2 * iz
    c2 = MVector{M}([zeros(eltype(SA),M-1);C2[2,2,iz,ID]])
    c3 = MVector{M}([zeros(eltype(SA),M-1);C3[2,2,iz,ID]])
    t = (r2'*c2 + r3'*c3)
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t

    j = 2 * iz - 1
    r2 = MVector{M}(zeros(eltype(SA),M))
    r3 = MVector{M}([B3[1,iz,ID];zeros(eltype(SA),M-1)])
    ldivFull!(iz,ID,SA,r3,Val(M))
    r2 -= (A23[:,:,iz,ID] * r3) * fac

    i = 2 * iz - 2
    c2 = MVector{M}([C2[1,1,iz,ID];zeros(eltype(SA),M-1)])
    c3 = MVector{M}([C3[1,1,iz,ID];zeros(eltype(SA),M-1)])
    t = (r2'*c2 + r3'*c3)
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t

    i = 2 * iz - 1
    c2 = MVector{M}([C2[2,1,iz,ID];zeros(eltype(SA),M-1)])
    c3 = MVector{M}([C3[2,1,iz,ID];zeros(eltype(SA),M-1)])
    t = (r2'*c2 + r3'*c3)
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t

    i = 2 * iz
    c2 = MVector{M}([zeros(eltype(SA),M-1);C2[2,2,iz,ID]])
    c3 = MVector{M}([zeros(eltype(SA),M-1);C3[2,2,iz,ID]])
    t = (r2'*c2 + r3'*c3)
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t

    j = 2 * iz 
    r2 = MVector{M}(zeros(eltype(SA),M))
    r3 = MVector{M}([zeros(eltype(SA),M-1);B3[2,iz,ID]])
    ldivFull!(iz,ID,SA,r3,Val(M))
    r2 -= (A23[:,:,iz,ID] * r3) * fac

    i = 2 * iz - 2
    c2 = MVector{M}([C2[1,1,iz,ID];zeros(eltype(SA),M-1)])
    c3 = MVector{M}([C3[1,1,iz,ID];zeros(eltype(SA),M-1)])
    t = (r2'*c2 + r3'*c3)
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t

    i = 2 * iz - 1
    c2 = MVector{M}([C2[2,1,iz,ID];zeros(eltype(SA),M-1)])
    c3 = MVector{M}([C3[2,1,iz,ID];zeros(eltype(SA),M-1)])
    t = (r2'*c2 + r3'*c3)
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t

    i = 2 * iz
    c2 = MVector{M}([zeros(eltype(SA),M-1);C2[2,2,iz,ID]])
    c3 = MVector{M}([zeros(eltype(SA),M-1);C3[2,2,iz,ID]])
    t = (r2'*c2 + r3'*c3)
    iB = 3 + 1 + i - j
    @atomic SchurBand[iB,j,ID] -= t
  end    
end  
