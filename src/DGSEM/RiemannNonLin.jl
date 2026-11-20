# The interface solvers are two-sided, RiemannSolver!(FL,FR,VLL,VRR,AuxL,AuxR,n1,n2,n3):
# the lower/left element gets -FL, the upper/right element +FR. Conservative
# fluxes have FL = FR (see TwoSided), nonconservative surface terms differ.
# Without atomics: the edge nodes of element nodes shared by several edges
# (corners) write into CacheAux.FC, RiemannGatherKernel! adds them in a fixed
# order, the result does not depend on the thread schedule.
function RiemannNonLinH(RiemannSolver,F,U,Aux,CacheAux,DG,Metric,Grid,NumberThreadGPU,NUMV,NAUX)
  backend = get_backend(F)
  DoF = DG.DoF
  DoFE = DG.DoFE
  M = DG.V.OrdPolyZ + 1
  NE = Grid.NumEdges
  Nz = Grid.nz
  # several edge nodes per workgroup, as in FluxVolumeNonLinV
  DoFEG = min(max(div(2*NumberThreadGPU,M*Nz),1),DoFE)
  group = (M,Nz,DoFEG,1)
  ndrange = (M,Nz,DoFE,NE)
  KRiemannNonLinH3Kernel! = RiemannNonLinH3Kernel!(backend,group)
  KRiemannNonLinH3Kernel!(RiemannSolver,F,CacheAux.FC,U,Aux,DG.GlobE,CacheAux.SlotE,Grid.EF,Metric.NH,
    Metric.VolSurfH,DG.wF,Grid.NumFaces,Val(NUMV),Val(NAUX);ndrange=ndrange)
  NS = length(CacheAux.IndS)
  if NS > 0
    NSG = min(max(div(NumberThreadGPU,M*Nz),1),NS)
    KRiemannGatherKernel! = RiemannGatherKernel!(backend,(M,Nz,NSG))
    KRiemannGatherKernel!(F,CacheAux.FC,CacheAux.IndS,CacheAux.PtrS,Val(NUMV);ndrange=(M,Nz,NS))
  end
end

@kernel inbounds = true function RiemannGatherKernel!(F,@Const(FC),@Const(IndS),@Const(PtrS),
  ::Val{NUMV}) where {NUMV}

  K, Iz, i = @index(Global, NTuple)

  NS = @uniform @ndrange()[3]

  if i <= NS
    ind = IndS[i]
    @unroll for iv = 1 : NUMV
      FSum = F[K,Iz,ind,iv]
      for l = PtrS[i] : PtrS[i+1] - 1
        FSum += FC[K,Iz,l,iv]
      end
      F[K,Iz,ind,iv] = FSum
    end
  end
end

function RiemannNonLinV(RiemannSolver,F,U,Aux,DG,Metric,Grid,NumberThreadGPU,NUMV,NAUX)
  backend = get_backend(F)
  NumI = DG.NumI
  M = DG.V.OrdPolyZ + 1
  Nz = Grid.nz
  NF = Grid.NumFaces
  NumIG = min(div(NumberThreadGPU,Nz+1),NumI)
  group = (Nz+1,NumIG)
  ndrange = (Nz+1,NumI)
  KRiemannNonLinV3Kernel! = RiemannNonLinV3Kernel!(backend,group)
  KRiemannNonLinV3Kernel!(RiemannSolver,F,U,Aux,Metric.NV,
    Metric.VolSurfV,Val(NUMV),Val(NAUX);ndrange=ndrange)
end  

@kernel inbounds = true function RiemannNonLinV3Kernel!(RiemannSolver!,F,@Const(U),
  @Const(Aux),@Const(NV),@Const(VolSurfV),
  ::Val{NUMV}, ::Val{NAUX}) where {NUMV, NAUX}

  Iz,ind = @index(Global, NTuple)


  Nz = @uniform @ndrange()[1]
  NumI = @uniform @ndrange()[2]

  VLL = @private eltype(F) (NUMV,)
  VRR = @private eltype(F) (NUMV,)
  AuxL = @private eltype(F) (NAUX,)
  AuxR = @private eltype(F) (NAUX,)
  FLoc = @private eltype(F) (NUMV,)
  FRLoc = @private eltype(F) (NUMV,)

  if ind <= NumI
    n1 = NV[1,Iz,ind]
    n2 = NV[2,Iz,ind]
    n3 = NV[3,Iz,ind]
    if Iz == 1
      @unroll for iAux = 1 : NAUX
        AuxL[iAux] = Aux[1,Iz,ind,iAux]
        AuxR[iAux] = AuxL[iAux]
      end
      @unroll for iv = 1 : NUMV
        VLL[iv] = U[1,Iz,ind,iv]
        VRR[iv] = VLL[iv]
      end
      t = eltype(F)(2) * (n1 * VLL[2] +
        n2 * VLL[3] +
        n3 * VLL[4])
      VLL[2] -= n1 * t
      VLL[3] -= n2 * t
      VLL[4] -= n3 * t
      RiemannSolver!(FLoc,FRLoc,VLL,VRR,AuxL,AuxR,n1,n2,n3)
      Surf = VolSurfV[Iz,ind]  
      @unroll for iv = 1 : NUMV  
        F[1,Iz,ind,iv] += FRLoc[iv] * Surf
      end  
    elseif Iz == Nz 
      @unroll for iAux = 1 : NAUX  
        AuxL[iAux] = Aux[end,Iz-1,ind,iAux]
        AuxR[iAux] = AuxL[iAux]
      end  
      @unroll for iv = 1 : NUMV  
        VLL[iv] = U[end,Iz-1,ind,iv]
        VRR[iv] = VLL[iv]
      end  
      t = eltype(F)(2) * (n1 * VRR[2] +
        n2 * VRR[3] +
        n3 * VRR[4])
      VRR[2] -= n1 * t
      VRR[3] -= n2 * t
      VRR[4] -= n3 * t
      RiemannSolver!(FLoc,FRLoc,VLL,VRR,AuxL,AuxR,n1,n2,n3)
      Surf = VolSurfV[Iz,ind]   
      @unroll for iv = 1 : NUMV  
        F[end,Iz-1,ind,iv] -= FLoc[iv] * Surf
      end  
    else
      @unroll for iAux = 1 : NAUX  
        AuxL[iAux] = Aux[end,Iz-1,ind,iAux]
        AuxR[iAux] = Aux[1,Iz,ind,iAux]
      end  
      @unroll for iv = 1 : NUMV  
        VLL[iv] = U[end,Iz-1,ind,iv]
        VRR[iv] = U[1,Iz,ind,iv]
      end  
      RiemannSolver!(FLoc,FRLoc,VLL,VRR,AuxL,AuxR,n1,n2,n3)
      Surf = VolSurfV[Iz,ind]  
      @unroll for iv = 1 : NUMV  
        F[end,Iz-1,ind,iv] -= FLoc[iv] * Surf
        F[1,Iz,ind,iv] += FRLoc[iv] * Surf
      end  
    end    
  end  
end

@kernel inbounds = true function RiemannNonLinH3Kernel!(RiemannSolver!,F,FC,@Const(U),@Const(Aux),@Const(GlobE),
  @Const(SlotE),@Const(EF),@Const(NH),@Const(VolSurfH),
  @Const(w), NF, ::Val{NUMV}, ::Val{NAUX}) where {NUMV, NAUX}

  _,iz, = @index(Local, NTuple)
  K, Iz, I, IE = @index(Global, NTuple)

  DoFE = @uniform @ndrange()[3]

  FLoc = @private eltype(F) (NUMV,)
  FRLoc = @private eltype(F) (NUMV,)
  UL = @private eltype(F) (NUMV,)
  UR = @private eltype(F) (NUMV,)
  AuxL = @private eltype(F) (NAUX,)
  AuxR = @private eltype(F) (NAUX,)

  uPos = @uniform 2
  vPos = @uniform 3
  wPos = @uniform 4


  if I <= DoFE
    iFL = EF[1,IE]
    iFR = EF[2,IE]
    indL = GlobE[1,I,IE]
    indR = GlobE[2,I,IE]
    n1 = NH[K,Iz,I,IE,1]
    n2 = NH[K,Iz,I,IE,2]
    n3 = NH[K,Iz,I,IE,3]
    if iFL > 0
      @unroll for iAux = 1 : NAUX  
        AuxL[iAux] = Aux[K,Iz,indL,iAux]
      end  
      @unroll for iv = 1 : NUMV  
        UL[iv] = U[K,Iz,indL,iv]
      end  
    else  
      @unroll for iAux = 1 : NAUX  
        AuxL[iAux] = Aux[K,Iz,indR,iAux]
      end  
      @unroll for iv = 1 : NUMV  
        UL[iv] = U[K,Iz,indR,iv]
      end  
      t = eltype(F)(2) * (n1 * UL[uPos] +
        n2 * UL[vPos] + n3 * UL[wPos]) 
      UL[uPos] -= n1 * t  
      UL[vPos] -= n2 * t  
      UL[wPos] -= n3 * t  
    end  
    if iFR > 0
      @unroll for iAux = 1 : NAUX  
        AuxR[iAux] = Aux[K,Iz,indR,iAux]
      end  
      @unroll for iv = 1 : NUMV  
        UR[iv] = U[K,Iz,indR,iv]
      end  
    else  
      @unroll for iAux = 1 : NAUX  
        AuxR[iAux] = Aux[K,Iz,indL,iAux]
      end  
      @unroll for iv = 1 : NUMV  
        UR[iv] = U[K,Iz,indL,iv]
      end  
      t = eltype(F)(2) * (n1 * UR[uPos] +
        n2 * UR[vPos] +
        n3 * UR[wPos]) 
      UR[uPos] -= n1 * t  
      UR[vPos] -= n2 * t  
      UR[wPos] -= n3 * t  
    end  
    RiemannSolver!(FLoc,FRLoc,UL,UR,AuxL,AuxR,n1,n2,n3)
    Surf = VolSurfH[K,Iz,I,IE] / w[I]  
    if 0 < iFL <= NF
      slot = SlotE[1,I,IE]
      if slot == 0
        @unroll for iv = 1 : NUMV  
          F[K,Iz,indL,iv] += -FLoc[iv] * Surf
        end  
      else
        @unroll for iv = 1 : NUMV  
          FC[K,Iz,slot,iv] = -FLoc[iv] * Surf
        end  
      end
    end  
    if 0 < iFR <= NF
      slot = SlotE[2,I,IE]
      if slot == 0
        @unroll for iv = 1 : NUMV  
          F[K,Iz,indR,iv] += FRLoc[iv] * Surf
        end  
      else
        @unroll for iv = 1 : NUMV  
          FC[K,Iz,slot,iv] = FRLoc[iv] * Surf
        end  
      end
    end  
  end  
end  
