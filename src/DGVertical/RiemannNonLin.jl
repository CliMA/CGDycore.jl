@kernel inbounds = true function RiemanNonLinVertKernel!(RiemannSolver!,NonConservativeFlux,F,
  @Const(U),@Const(Aux),@Const(w), ::Val{M}, ::Val{NUMV}, ::Val{NAUX}) where {M, NUMV, NAUX}

  Iz = @index(Global, NTuple)

  Nz = @uniform @ndrange()[1]

  VLL = @private eltype(F) (NUMV,)
  VRR = @private eltype(F) (NUMV,)
  AuxL = @private eltype(F) (NAUX,)
  AuxR = @private eltype(F) (NAUX,)
  FLoc = @private eltype(F) (NUMV,)

  RhoPos = @uniform 1
  uPos = @uniform 2
  vPos = @uniform 3
  wPos = @uniform 4
  ThPos = @uniform 5

  if Iz <= Nz
    ind = Glob[ID,IF]
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
      VLL[uPos] = -VLL[uPos]
      VLL[vPos] = -VLL[vPos]
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
      VRR[uPos] = -VRR[uPos]
      VRR[vPos] = -VRR[vPos]
      VRR[wPos] = -VRR[wPos]
    end

    @views RiemannSolver!(FLoc,VLL,VRR,AuxL,AuxR,NV[:,ID,Iz,IF])

    Surf = eltype(U)(1) / w[1]  
    FLoc[RhoPos] *= Surf
    FLoc[uPos] *= Surf
    FLoc[vPos] *= Surf
    FLoc[wPos] *= Surf
    FLoc[ThPos] *= Surf
    if Iz > 1 
      @atomic :monotonic F[Iz-1,M,RhoPos] -= FLoc[RhoPos]
      @atomic :monotonic F[Iz-1,M,uPos] -= FLoc[uPos] 
      @atomic :monotonic F[Iz-1,M,vPos] -= FLoc[vPos]
      @atomic :monotonic F[Iz-1,M,wPos] -= FLoc[wPos]
      @atomic :monotonic F[Iz-1,M,ThPos] -= FLoc[ThPos]
    end  
    if Iz < Nz
      @atomic :monotonic F[Iz,1,RhoPos] += FLoc[RhoPos]
      @atomic :monotonic F[Iz,1,uPos] += FLoc[uPos]
      @atomic :monotonic F[Iz,1,vPos] += FLoc[vPos]
      @atomic :monotonic F[Iz,1,wPos] += FLoc[wPos]
      @atomic :monotonic F[Iz,1,ThPos] += FLoc[ThPos]
    end  
  end  
end

