@kernel inbounds = true function RiemanNonLinV3Kernel!(RiemannSolver!,NonConservativeFlux,F,@Const(U),@Const(Aux),@Const(Glob),
  @Const(NV),@Const(VolSurfV),
  @Const(w), ::Val{M}, ::Val{NUMV}, ::Val{NAUX}) where {M, NUMV, NAUX}

  Iz,iD,_  = @index(Local, NTuple)
  Iz,ID,IF = @index(Global, NTuple)


  Nz = @uniform @ndrange()[1]
  NQ = @uniform @ndrange()[2]

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

  if ID <= NQ
    ind = Glob[ID,IF]
    if Iz > 1    
      for iAux = 1 : NAUX  
        AuxL[iAux] = Aux[Iz-1,M,ind,iAux]
      end  
      @unroll for iv = 1 : NUMV  
        VLL[iv] = U[Iz-1,M,ind,iv]
      end  
    else
      for iAux = 1 : NAUX  
        AuxL[iAux] = Aux[Iz,1,ind,iAux]
      end  
      @unroll for iv = 1 : NUMV  
        VLL[iv] = U[Iz,1,ind,iv]
      end  
      VLL[uPos] = -VLL[uPos]
      VLL[vPos] = -VLL[vPos]
      VLL[wPos] = -VLL[wPos]
    end  
    if Iz < Nz
      for iAux = 1 : NAUX  
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

    Surf = VolSurfV[ID,Iz,IF] / w[1]  
    FLoc[RhoPos] *= Surf
    FLoc[uPos] *= Surf
    FLoc[vPos] *= Surf
    FLoc[wPos] *= Surf
    FLoc[ThPos] *= Surf
    if Iz > 1 
      @atomic :monotonic F[Iz-1,M,ind,RhoPos] -= FLoc[RhoPos]
      @atomic :monotonic F[Iz-1,M,ind,uPos] -= FLoc[uPos] 
      @atomic :monotonic F[Iz-1,M,ind,vPos] -= FLoc[vPos]
      @atomic :monotonic F[Iz-1,M,ind,wPos] -= FLoc[wPos]
      @atomic :monotonic F[Iz-1,M,ind,ThPos] -= FLoc[ThPos]
    end  
    if Iz < Nz
      @atomic :monotonic F[Iz,1,ind,RhoPos] += FLoc[RhoPos]
      @atomic :monotonic F[Iz,1,ind,uPos] += FLoc[uPos]
      @atomic :monotonic F[Iz,1,ind,vPos] += FLoc[vPos]
      @atomic :monotonic F[Iz,1,ind,wPos] += FLoc[wPos]
      @atomic :monotonic F[Iz,1,ind,ThPos] += FLoc[ThPos]
    end  
  end  
end

@kernel inbounds = true function RiemanNonLinH3Kernel!(RiemannSolver!,F,@Const(U),@Const(Aux),@Const(GlobE),
  @Const(EF),@Const(FTE),@Const(NH),@Const(VolSurfH),
  @Const(w), NF, ::Val{NUMV}, ::Val{NAUX}) where {NUMV, NAUX}

  _,_,iz, = @index(Local, NTuple)
  I,K,Iz,IE = @index(Global, NTuple)

  N = @uniform @groupsize()[1]
  M = @uniform @groupsize()[2]
  TilesDim = @uniform @groupsize()[3]

  Nz = @uniform @ndrange()[3]

  FLoc = @private eltype(F) (NUMV,)

  RhoPos = @uniform 1
  uPos = @uniform 2
  vPos = @uniform 3
  wPos = @uniform 4
  ThPos = @uniform 5


  if Iz <= Nz
    iFL = EF[1,IE]
    iFR = EF[2,IE]
    indL = GlobE[1,I,IE]
    indR = GlobE[2,I,IE]
    RiemannSolver!(FLoc,view(U[Iz,K,indL,:],1:NUMV),view(U[Iz,K,indR,:],1:NUMV),
      view(Aux[Iz,K,indL,:],1:NAUX),view(Aux[Iz,K,indR,:],1:NAUX),
      view(NH[:,K,I,Iz,IE],1:3))
    Surf = VolSurfH[K,I,Iz,IE] / w[1]  
    FLoc[RhoPos] *= Surf
    FLoc[uPos] *= Surf
    FLoc[vPos] *= Surf
    FLoc[wPos] *= Surf
    FLoc[ThPos] *= Surf
    if iFL <= NF
      @atomic :monotonic F[Iz,K,indL,RhoPos] -= FLoc[RhoPos]
      @atomic :monotonic F[Iz,K,indL,uPos] -= FLoc[uPos]
      @atomic :monotonic F[Iz,K,indL,vPos] -= FLoc[vPos]
      @atomic :monotonic F[Iz,K,indL,wPos] -= FLoc[wPos]
      @atomic :monotonic F[Iz,K,indL,ThPos] -= FLoc[ThPos]
    end  
    if iFR <= NF
      @atomic :monotonic F[Iz,K,indR,RhoPos] += FLoc[RhoPos]
      @atomic :monotonic F[Iz,K,indR,uPos] += FLoc[uPos]
      @atomic :monotonic F[Iz,K,indR,vPos] += FLoc[vPos]
      @atomic :monotonic F[Iz,K,indR,wPos] += FLoc[wPos]
      @atomic :monotonic F[Iz,K,indR,ThPos] += FLoc[ThPos]
    end  
  end  
end

