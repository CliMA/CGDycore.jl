function RiemannNonLinH(RiemannSolver,F,U,Aux,DG,Metric,Grid,NumberThreadGPU,NV,NAUX)
  backend = get_backend(F)
  DoF = DG.DoF
  DoFE = DG.DoFE
  M = DG.OrdPolyZ + 1
  NE = Grid.NumEdges
  Nz = Grid.nz
  NEG = min(div(NumberThreadGPU,DoFE*M),Nz)
  group = (DoFE,M,NEG,1)
  ndrange = (DoFE,M,Nz,NE)
  KRiemannNonLinH3Kernel! = RiemannNonLinH3Kernel!(backend,group)
  KRiemannNonLinH3Kernel!(RiemannSolver,F,U,Aux,DG.GlobE,Grid.EF,Grid.FE,Metric.NH,
    Metric.VolSurfH,DG.wF,Grid.NumFaces,Val(NV),Val(NAUX);ndrange=ndrange)
end


function RiemannNonLinV(RiemannSolver,NonConservativeFlux,F,U,Aux,DG,Metric,Grid,NumberThreadGPU,NV,NAUX)
  backend = get_backend(F)
  DoF = DG.DoF
  M = DG.OrdPolyZ + 1
  Nz = Grid.nz
  NF = Grid.NumFaces
  DoFG = min(div(NumberThreadGPU,Nz+1),DoF)
  group = (Nz+1,DoFG,1)
  ndrange = (Nz+1,DoF,NF)
  KRiemannNonLinV3Kernel! = RiemannNonLinV3Kernel!(backend,group)
  KRiemannNonLinV3Kernel!(RiemannSolver,NonConservativeFlux,F,U,Aux,DG.Glob,Metric.NV,
    Metric.VolSurfV,DG.wZ,Val(M),Val(NV),Val(NAUX);ndrange=ndrange)
end  

@kernel inbounds = true function RiemannNonLinV3Kernel!(RiemannSolver!,NonConservativeFlux,F,@Const(U),@Const(Aux),@Const(Glob),
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
      @unroll for iAux = 1 : NAUX  
        AuxL[iAux] = Aux[M,Iz-1,ind,iAux]
      end  
      @unroll for iv = 1 : NUMV  
        VLL[iv] = U[M,Iz-1,ind,iv]
      end  
    else
      @unroll for iAux = 1 : NAUX  
        AuxL[iAux] = Aux[1,Iz,ind,iAux]
      end  
      @unroll for iv = 1 : NUMV  
        VLL[iv] = U[1,Iz,ind,iv]
      end  
      VLL[uPos] = -VLL[uPos]
      VLL[vPos] = -VLL[vPos]
      VLL[wPos] = -VLL[wPos]
    end  
    if Iz < Nz
      @unroll for iAux = 1 : NAUX  
        AuxR[iAux] = Aux[1,Iz,ind,iAux]
      end  
      @unroll for iv = 1 : NUMV  
        VRR[iv] = U[1,Iz,ind,iv]
      end  
    else  
      @unroll for iAux = 1 : NAUX  
        AuxR[iAux] = Aux[M,Iz-1,ind,iAux]
      end  
      @unroll for iv = 1 : NUMV  
        VRR[iv] = U[M,Iz-1,ind,iv]
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
      @atomic :monotonic F[M,Iz-1,ind,RhoPos] -= FLoc[RhoPos]
      @atomic :monotonic F[M,Iz-1,ind,uPos] -= FLoc[uPos] 
      @atomic :monotonic F[M,Iz-1,ind,vPos] -= FLoc[vPos]
      @atomic :monotonic F[M,Iz-1,ind,wPos] -= FLoc[wPos]
      @atomic :monotonic F[M,Iz-1,ind,ThPos] -= FLoc[ThPos]
    end  
    if Iz < Nz
      @atomic :monotonic F[1,Iz,ind,RhoPos] += FLoc[RhoPos]
      @atomic :monotonic F[1,Iz,ind,uPos] += FLoc[uPos]
      @atomic :monotonic F[1,Iz,ind,vPos] += FLoc[vPos]
      @atomic :monotonic F[1,Iz,ind,wPos] += FLoc[wPos]
      @atomic :monotonic F[1,Iz,ind,ThPos] += FLoc[ThPos]
    end  
  end  
end

@kernel inbounds = true function RiemannNonLinH3Kernel!(RiemannSolver!,F,@Const(U),@Const(Aux),@Const(GlobE),
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
    RiemannSolver!(FLoc,view(U,K,Iz,indL,1:NUMV),view(U,K,Iz,indR,1:NUMV),
      view(Aux,K,Iz,indL,1:NAUX),view(Aux,K,Iz,indR,1:NAUX),
      view(NH,1:3,K,I,Iz,IE))
    Surf = VolSurfH[K,I,Iz,IE] / w[I]  
    FLoc[RhoPos] = FLoc[RhoPos] * Surf
    FLoc[uPos] = FLoc[uPos] * Surf
    FLoc[vPos] = FLoc[vPos] * Surf
    FLoc[wPos] = FLoc[wPos] * Surf
    FLoc[ThPos] = FLoc[ThPos] * Surf
    if iFL <= NF
      @atomic :monotonic F[K,Iz,indL,RhoPos] += -FLoc[RhoPos]
      @atomic :monotonic F[K,Iz,indL,uPos] += -FLoc[uPos]
      @atomic :monotonic F[K,Iz,indL,vPos] += -FLoc[vPos]
      @atomic :monotonic F[K,Iz,indL,wPos] += -FLoc[wPos]
      @atomic :monotonic F[K,Iz,indL,ThPos] += -FLoc[ThPos]
    end  
    if iFR <= NF
      @atomic :monotonic F[K,Iz,indR,RhoPos] += FLoc[RhoPos]
      @atomic :monotonic F[K,Iz,indR,uPos] += FLoc[uPos]
      @atomic :monotonic F[K,Iz,indR,vPos] += FLoc[vPos]
      @atomic :monotonic F[K,Iz,indR,wPos] += FLoc[wPos]
      @atomic :monotonic F[K,Iz,indR,ThPos] += FLoc[ThPos]
    end  
  end  
end


@kernel inbounds = true function RiemannNonLinV3NonConservativeKernel!(RiemannSolver!,NonConservativeFlux,F,@Const(U),@Const(Aux),@Const(Glob),
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
  FLocL = @private eltype(F) (NUMV,)
  FLocR = @private eltype(F) (NUMV,)

  RhoPos = @uniform 1
  uPos = @uniform 2
  vPos = @uniform 3
  wPos = @uniform 4
  ThPos = @uniform 5

  if ID <= NQ
    ind = Glob[ID,IF]
    if Iz > 1    
      @unroll for iAux = 1 : NAUX  
        AuxL[iAux] = Aux[M,Iz-1,ind,iAux]
      end  
      @unroll for iv = 1 : NUMV  
        VLL[iv] = U[M,Iz-1,ind,iv]
      end  
    else
      @unroll for iAux = 1 : NAUX  
        AuxL[iAux] = Aux[1,Iz,ind,iAux]
      end  
      @unroll for iv = 1 : NUMV  
        VLL[iv] = U[1,Iz,ind,iv]
      end  
      VLL[uPos] = -VLL[uPos]
      VLL[vPos] = -VLL[vPos]
      VLL[wPos] = -VLL[wPos]
    end  
    if Iz < Nz
      @unroll for iAux = 1 : NAUX  
        AuxR[iAux] = Aux[1,Iz,ind,iAux]
      end  
      @unroll for iv = 1 : NUMV  
        VRR[iv] = U[1,Iz,ind,iv]
      end  
    else  
      @unroll for iAux = 1 : NAUX  
        AuxR[iAux] = Aux[M,Iz-1,ind,iAux]
      end  
      @unroll for iv = 1 : NUMV  
        VRR[iv] = U[M,Iz-1,ind,iv]
      end  
      VRR[uPos] = -VRR[uPos]
      VRR[vPos] = -VRR[vPos]
      VRR[wPos] = -VRR[wPos]
    end

    @views RiemannSolver!(FLocL,FLocR,VLL,VRR,AuxL,AuxR,NV[:,ID,Iz,IF])

    Surf = VolSurfV[ID,Iz,IF] / w[1]  
    FLocL[RhoPos] *= Surf
    FLocL[uPos] *= Surf
    FLocL[vPos] *= Surf
    FLocL[wPos] *= Surf
    FLocL[ThPos] *= Surf
    
    FLocR[RhoPos] *= Surf
    FLocR[uPos] *= Surf
    FLocR[vPos] *= Surf
    FLocR[wPos] *= Surf
    FLocR[ThPos] *= Surf
    if Iz > 1 
      @atomic :monotonic F[M,Iz-1,ind,RhoPos] -= FLocR[RhoPos]
      @atomic :monotonic F[M,Iz-1,ind,uPos] -= FLocR[uPos] 
      @atomic :monotonic F[M,Iz-1,ind,vPos] -= FLocR[vPos]
      @atomic :monotonic F[M,Iz-1,ind,wPos] -= FLocR[wPos]
      @atomic :monotonic F[M,Iz-1,ind,ThPos] -= FLocR[ThPos]
    end  
    if Iz < Nz
      @atomic :monotonic F[1,Iz,ind,RhoPos] += FLocL[RhoPos]
      @atomic :monotonic F[1,Iz,ind,uPos] += FLocL[uPos]
      @atomic :monotonic F[1,Iz,ind,vPos] += FLocL[vPos]
      @atomic :monotonic F[1,Iz,ind,wPos] += FLocL[wPos]
      @atomic :monotonic F[1,Iz,ind,ThPos] += FLocL[ThPos]
    end  
  end  
end

@kernel inbounds = true function RiemannNonLinH3NonConservativeKernel!(RiemannSolver!,F,@Const(U),@Const(Aux),@Const(GlobE),
  @Const(EF),@Const(FTE),@Const(NH),@Const(VolSurfH),
  @Const(w), NF, ::Val{NUMV}, ::Val{NAUX}) where {NUMV, NAUX}

  _,_,iz, = @index(Local, NTuple)
  I,K,Iz,IE = @index(Global, NTuple)

  N = @uniform @groupsize()[1]
  M = @uniform @groupsize()[2]
  TilesDim = @uniform @groupsize()[3]

  Nz = @uniform @ndrange()[3]

  FLocL = @private eltype(F) (NUMV,)
  FLocR = @private eltype(F) (NUMV,)
# FLoc = @localmem eltype(F) (N,M,TilesDim,NUMV)

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
    RiemannSolver!(FLocL,FLocR,view(U,K,Iz,indL,1:NUMV),view(U,K,Iz,indR,1:NUMV),
      view(Aux,K,Iz,indL,1:NAUX),view(Aux,K,Iz,indR,1:NAUX),
      view(NH,1:3,K,I,Iz,IE))
    Surf = VolSurfH[K,I,Iz,IE] / w[I]  
    FLocL[RhoPos] = FLocL[RhoPos] * Surf
    FLocL[uPos] = FLocL[uPos] * Surf
    FLocL[vPos] = FLocL[vPos] * Surf
    FLocL[wPos] = FLocL[wPos] * Surf
    FLocL[ThPos] = FLocL[ThPos] * Surf
    FLocR[RhoPos] = FLocR[RhoPos] * Surf
    FLocR[uPos] = FLocR[uPos] * Surf
    FLocR[vPos] = FLocR[vPos] * Surf
    FLocR[wPos] = FLocR[wPos] * Surf
    FLocR[ThPos] = FLocR[ThPos] * Surf
    if iFL <= NF
      @atomic :monotonic F[K,Iz,indL,RhoPos] += -FLocR[RhoPos]
      @atomic :monotonic F[K,Iz,indL,uPos] += -FLocR[uPos]
      @atomic :monotonic F[K,Iz,indL,vPos] += -FLocR[vPos]
      @atomic :monotonic F[K,Iz,indL,wPos] += -FLocR[wPos]
      @atomic :monotonic F[K,Iz,indL,ThPos] += -FLocR[ThPos]
    end  
    if iFR <= NF
      @atomic :monotonic F[K,Iz,indR,RhoPos] += FLocL[RhoPos]
      @atomic :monotonic F[K,Iz,indR,uPos] += FLocL[uPos]
      @atomic :monotonic F[K,Iz,indR,vPos] += FLocL[vPos]
      @atomic :monotonic F[K,Iz,indR,wPos] += FLocL[wPos]
      @atomic :monotonic F[K,Iz,indR,ThPos] += FLocL[ThPos]
    end  
  end  
end

