abstract type RiemannSolver end

Base.@kwdef struct RiemannLMARS <: RiemannSolver end

function (::RiemannLMARS)(Param,Phys,hPos,uPos,vPos,wPos,pPos)
  @inline function RiemannByLMARSNonLin!(F,VLL,VRR,AuxL,AuxR)

    FT = eltype(F)
    cS = Param.cS
    pLL = AuxL[pPos]
    pRR = AuxR[pPos]
    hM = FT(0.5) * (VLL[hPos] + VRR[hPos])
    vLL = VLL[uPos] / VLL[hPos]
    vRR = VRR[uPos] / VRR[hPos]
    pM = FT(0.5) * (pLL + pRR) - FT(0.5) * cS * hM * (vRR - vLL)
    vM = FT(0.5) * (vRR + vLL) - FT(1.0) /(FT(2.0) * cS) * (pRR - pLL) / hM
    if vM > FT(0)
      F[hPos] = vM * VLL[hPos]
      F[uPos] = vM * VLL[uPos] + pM
      F[vPos] = vM * VLL[vPos]
      F[wPos] = vM * VLL[wPos]
    else
      F[hPos] = vM * VRR[hPos]
      F[uPos] = vM * VRR[uPos] + pM
      F[vPos] = vM * VRR[vPos]
      F[wPos] = vM * VRR[wPos]
    end
  end
  return RiemannByLMARSNonLin!
end    

function (::RiemannLMARS)(Param,Phys,hPos,uPos,vPos,wPos,ThPos,pPos)
  @inline function RiemannByLMARSNonLin!(F,VLL,VRR,AuxL,AuxR)

    FT = eltype(F)
    cS = Param.cS
    pLL = AuxL[pPos]
    pRR = AuxR[pPos]
    hM = FT(0.5) * (VLL[hPos] + VRR[hPos])
    vLL = VLL[uPos] / VLL[hPos]
    vRR = VRR[uPos] / VRR[hPos]
    pM = FT(0.5) * (pLL + pRR) - FT(0.5) * cS * hM * (vRR - vLL)
    vM = FT(0.5) * (vRR + vLL) - FT(1.0) /(FT(2.0) * cS) * (pRR - pLL) / hM
    if vM > FT(0)
      F[hPos] = vM * VLL[hPos]
      F[uPos] = vM * VLL[uPos] + pM
      F[vPos] = vM * VLL[vPos]
      F[wPos] = vM * VLL[wPos]
      F[ThPos] = vM * VLL[ThPos]
    else
      F[hPos] = vM * VRR[hPos]
      F[uPos] = vM * VRR[uPos] + pM
      F[vPos] = vM * VRR[vPos]
      F[wPos] = vM * VRR[wPos]
      F[ThPos] = vM * VRR[ThPos]
    end
  end
  return RiemannByLMARSNonLin!
end

@kernel inbounds = true function RiemanNonLinKernel!(RiemannSolver!,F,@Const(U),@Const(Aux),@Const(GlobE),
  @Const(EF),@Const(FTE),@Const(NH),@Const(T1H),@Const(T2H),@Const(VolSurfH),
  @Const(w), NF, ::Val{NV}, ::Val{NAUX}) where {NV, NAUX}

  _,_,I,iE = @index(Local, NTuple)
  Iz,K,_,IE = @index(Global, NTuple)

  N = @uniform @groupsize()[1]
  TilesDim = @uniform @groupsize()[2]

  NE = @uniform @ndrange()[2]

  VLL = @localmem eltype(F) (N,TilesDim,NV)
  VRR = @localmem eltype(F) (N,TilesDim,NV)
  FLoc = @private eltype(F) (NV,)
  FE = @private eltype(F) (NV,)

  hPos = @uniform 1
  uPos = @uniform 2
  vPos = @uniform 3
  wPos = @uniform 4


  if IE <= NE
    iFL = EF[1,IE]
    iFR = EF[2,IE]
    indL = GlobE[1,I,IE]
    indR = GlobE[2,I,IE]
    VLL[I,iE,hPos] = U[1,1,indL,hPos]
    VLL[I,iE,uPos] = NH[1,1,I,1,IE] * U[1,1,indL,uPos] +
      NH[2,1,I,1,IE] * U[1,1,indL,vPos] +
      NH[3,1,I,1,IE] * U[1,1,indL,wPos]
    VLL[I,iE,vPos] = T1H[1,1,I,1,IE] * U[1,1,indL,uPos] +
      T1H[2,1,I,1,IE] * U[1,1,indL,vPos] +
      T1H[3,1,I,1,IE] * U[1,1,indL,wPos]
    VLL[I,iE,wPos] = T2H[1,1,I,1,IE] * U[1,1,indL,uPos] +
      T2H[2,1,I,1,IE] * U[1,1,indL,vPos] +
      T2H[3,1,I,1,IE] * U[1,1,indL,wPos]
    VRR[I,iE,hPos] = U[1,1,indR,hPos]
    VRR[I,iE,uPos] = NH[1,1,I,1,IE] * U[1,1,indR,uPos] +
      NH[2,1,I,1,IE] * U[1,1,indR,vPos] +
      NH[3,1,I,1,IE] * U[1,1,indR,wPos]
    VRR[I,iE,vPos] = T1H[1,1,I,1,IE] * U[1,1,indR,uPos] +
      T1H[2,1,I,1,IE] * U[1,1,indR,vPos] +
      T1H[3,1,I,1,IE] * U[1,1,indR,wPos]
    VRR[I,iE,wPos] = T2H[1,1,I,1,IE] * U[1,1,indR,uPos] +
      T2H[2,1,I,1,IE] * U[1,1,indR,vPos] +
      T2H[3,1,I,1,IE] * U[1,1,indR,wPos]
    @views RiemannSolver!(FLoc,VLL[I,iE,:],VRR[I,iE,:],
      Aux[1,1,indL,:],Aux[1,1,indR,:])

    FE[hPos] =  FLoc[hPos]
    FE[uPos] =  NH[1,1,I,1,IE] * FLoc[uPos] +
      T1H[1,1,I,1,IE] * FLoc[vPos] + T2H[1,1,I,1,IE] * FLoc[wPos]
    FE[vPos] =  NH[2,1,I,1,IE] * FLoc[uPos] +
      T1H[2,1,I,1,IE] * FLoc[vPos] + T2H[2,1,I,1,IE] * FLoc[wPos]
    FE[wPos] =  NH[3,1,I,1,IE] * FLoc[uPos] +
      T1H[3,1,I,1,IE] * FLoc[vPos] + T2H[3,1,I,1,IE] * FLoc[wPos]
    Surf = VolSurfH[1,I,1,IE] / w  
    FE[hPos] *= Surf
    FE[uPos] *= Surf
    FE[vPos] *= Surf
    FE[wPos] *= Surf
    if iFL <= NF
      @atomic :monotonic F[1,1,indL,hPos] -= FE[hPos]
      @atomic :monotonic F[1,1,indL,uPos] -= FE[uPos]
      @atomic :monotonic F[1,1,indL,vPos] -= FE[vPos]
      @atomic :monotonic F[1,1,indL,wPos] -= FE[wPos]
    end  
    if iFR <= NF
      @atomic :monotonic F[1,1,indR,hPos] += FE[hPos]
      @atomic :monotonic F[1,1,indR,uPos] += FE[uPos]
      @atomic :monotonic F[1,1,indR,vPos] += FE[vPos]
      @atomic :monotonic F[1,1,indR,wPos] += FE[wPos]
    end  
  end  
end

@kernel inbounds = true function VSp2VCartKernel!(VCart,@Const(VSp),@Const(Rotate))

  K,iQ,  = @index(Local, NTuple)
  _,Iz,IF = @index(Global, NTuple)

  M = @uniform @ndrange()[1]
  NQ = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  if IF <= NF
    iD = iQ + (IF - 1) * NQ
    VCart[1,1,iD,1] = VSp[1,1,iD,1]
    VCart[1,1,iD,2] = Rotate[1,1,1,iQ,1,IF] * VSp[1,1,iD,2] +
      Rotate[2,1,1,iQ,1,IF] * VSp[1,1,iD,3]
    VCart[1,1,iD,3] = Rotate[1,2,1,iQ,1,IF] * VSp[1,1,iD,2] +
      Rotate[2,2,1,iQ,1,IF] * VSp[1,1,iD,3]
    VCart[1,1,iD,4] = Rotate[1,3,1,iQ,1,IF] * VSp[1,1,iD,2] +
      Rotate[2,3,1,iQ,1,IF] * VSp[1,1,iD,3]
  end  
end


@kernel inbounds = true function VCart2VSpKernel!(VSp,@Const(VCart),@Const(Rotate),@Const(J))

  iQ,  = @index(Local, NTuple)
  _,IF = @index(Global, NTuple)

  NQ = @uniform @ndrange()[1]
  NF = @uniform @ndrange()[2]

  if IF <= NF
    iD = iQ + (IF - 1) * NQ
    VSp[1,1,iD,1] = VCart[1,1,iD,1] / J[iQ,1,1,IF]
    VSp[1,1,iD,2] = (Rotate[1,1,1,iQ,1,IF] * VCart[1,1,iD,2] +
      Rotate[1,2,1,iQ,1,IF] * VCart[1,1,iD,3] +
      Rotate[1,3,1,iQ,1,IF] * VCart[1,1,iD,4]) / J[iQ,1,1,IF]
    VSp[1,1,iD,3] = (Rotate[2,1,1,iQ,1,IF] * VCart[1,1,iD,2] +
      Rotate[2,2,1,iQ,1,IF] * VCart[1,1,iD,3] +
      Rotate[2,3,1,iQ,1,IF] * VCart[1,1,iD,4]) / J[iQ,1,1,IF]
  end  
end
