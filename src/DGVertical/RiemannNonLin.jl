@kernel inbounds = true function RiemanNonLinVertKernel!(RiemannSolver!,F,
  @Const(U),@Const(Aux),@Const(w), ::Val{M}, ::Val{NUMV}, ::Val{NAUX}) where {M, NUMV, NAUX}

  Iz = @index(Global, NTuple)

  Nz = @uniform @ndrange()[1]

  VLL = @private eltype(F) (NUMV,)
  VRR = @private eltype(F) (NUMV,)
  AuxL = @private eltype(F) (NAUX,)
  AuxR = @private eltype(F) (NAUX,)
  FLoc = @private eltype(F) (NUMV,)

  RhoPos = @uniform 1
  wPos = @uniform 2
  ThPos = @uniform 3

  if Iz <= Nz
    if Iz > 1    
      @unroll for iAux = 1 : NAUX  
        AuxL[iAux] = Aux[Iz-1,M,iAux]
      end  
      @unroll for iv = 1 : NUMV  
        VLL[iv] = U[Iz-1,M,iv]
      end  
    else
      @unroll for iAux = 1 : NAUX  
        AuxL[iAux] = Aux[Iz,1,iAux]
      end  
      @unroll for iv = 1 : NUMV  
        VLL[iv] = U[Iz,1,iv]
      end  
      VLL[wPos] = -VLL[wPos]
    end  
    if Iz < Nz
      @unroll for iAux = 1 : NAUX  
        AuxR[iAux] = Aux[Iz,1,ind,iAux]
      end  
      @unroll for iv = 1 : NUMV  
        VRR[iv] = U[Iz,1,ind,iv]
      end  
    else  
      @unroll for iAux = 1 : NAUX  
        AuxR[iAux] = Aux[Iz-1,M,ind,iAux]
      end  
      @unroll for iv = 1 : NUMV  
        VRR[iv] = U[Iz-1,M,ind,iv]
      end  
      VRR[wPos] = -VRR[wPos]
    end

    @views RiemannSolver!(FLoc,VLL,VRR,AuxL,AuxR)

    Surf = eltype(U)(1) / w[1]  
    FLoc[RhoPos] *= Surf
    FLoc[wPos] *= Surf
    FLoc[ThPos] *= Surf
    if Iz > 1 
      @atomic :monotonic F[Iz-1,M,RhoPos] -= FLoc[RhoPos]
      @atomic :monotonic F[Iz-1,M,wPos] -= FLoc[wPos]
      @atomic :monotonic F[Iz-1,M,ThPos] -= FLoc[ThPos]
    end  
    if Iz < Nz
      @atomic :monotonic F[Iz,1,RhoPos] += FLoc[RhoPos]
      @atomic :monotonic F[Iz,1,wPos] += FLoc[wPos]
      @atomic :monotonic F[Iz,1,ThPos] += FLoc[ThPos]
    end  
  end  
end

