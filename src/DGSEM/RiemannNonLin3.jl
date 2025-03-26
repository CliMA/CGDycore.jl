@kernel inbounds = true function RiemanNonLinV3Kernel!(RiemannSolver!,F,@Const(U),@Const(Aux),@Const(Glob),
  @Const(NV),@Const(T1V),@Const(T2V),@Const(VolSurfV),
  @Const(w), M, ::Val{NUMV}, ::Val{NAUX}) where {NUMV, NAUX}

  Iz,iD,_  = @index(Local, NTuple)
  Iz,ID,IF = @index(Global, NTuple)

  N = @uniform @groupsize()[1]
  TilesDim = @uniform @groupsize()[2]

  Nz = @uniform @ndrange()[1]
  NQ = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  VLL = @localmem eltype(F) (Nz,TilesDim,NUMV)
  VRR = @localmem eltype(F) (Nz,TilesDim,NUMV)
  AuxL = @localmem eltype(F) (Nz,TilesDim,NAUX)
  VRR = @localmem eltype(F) (Nz,TilesDim,NUMV)
  AuxR = @localmem eltype(F) (Nz,TilesDim,NAUX)
  FLoc = @private eltype(F) (NUMV,)
  FE = @private eltype(F) (NUMV,)

  RhoPos = @uniform 1
  uPos = @uniform 2
  vPos = @uniform 3
  wPos = @uniform 4
  ThPos = @uniform 5


  if ID <= NQ
    ind = Glob[ID,IF]
    if Iz > 1    
      for iAux = 1 : NAUX  
        AuxL[Iz,iD,iAux] = Aux[Iz-1,M,ind,iAux]
      end  
      VLL[Iz,iD,RhoPos] = U[Iz-1,M,ind,RhoPos]
      VLL[Iz,iD,uPos] = NV[1,ID,Iz,IF] * U[Iz-1,M,ind,uPos] +
        NV[2,ID,Iz,IF] * U[Iz-1,M,ind,vPos] +
        NV[3,ID,Iz,IF] * U[Iz-1,M,ind,wPos]
      VLL[Iz,iD,vPos] = T1V[1,ID,Iz,IF] * U[Iz-1,M,ind,uPos] +
        T1V[2,ID,Iz,IF] * U[Iz-1,M,ind,vPos] +
        T1V[3,ID,Iz,IF] * U[Iz-1,M,ind,wPos]
      VLL[Iz,iD,wPos] = T2V[1,ID,Iz,IF] * U[Iz-1,M,ind,uPos] +
        T2V[2,ID,Iz,IF] * U[Iz-1,M,ind,vPos] +
        T2V[3,ID,Iz,IF] * U[Iz-1,M,ind,wPos]
      VLL[Iz,iD,ThPos] = U[Iz-1,M,ind,ThPos]
    else
      for iAux = 1 : NAUX  
        AuxL[Iz,iD,iAux] = Aux[Iz,1,ind,iAux]
      end  
      VLL[Iz,iD,RhoPos] = U[Iz,1,ind,RhoPos]
      VLL[Iz,iD,uPos] = -(NV[1,ID,Iz,IF] * U[Iz,1,ind,uPos] +
        NV[2,ID,Iz,IF] * U[Iz,1,ind,vPos] +
        NV[3,ID,Iz,IF] * U[Iz,1,ind,wPos])
      VLL[Iz,iD,vPos] = T1V[1,ID,Iz,IF] * U[Iz,1,ind,uPos] +
        T1V[2,ID,Iz,IF] * U[Iz,1,ind,vPos] +
        T1V[3,ID,Iz,IF] * U[Iz,1,ind,wPos]
      VLL[Iz,iD,wPos] = T2V[1,ID,Iz,IF] * U[Iz,1,ind,uPos] +
        T2V[2,ID,Iz,IF] * U[Iz,1,ind,vPos] +
        T2V[3,ID,Iz,IF] * U[Iz,1,ind,wPos]
    end  
    if Iz < Nz
      for iAux = 1 : NAUX  
        AuxR[Iz,iD,iAux] = Aux[Iz,1,ind,iAux]
      end  
      VRR[Iz,iD,RhoPos] = U[Iz,1,ind,RhoPos]
      VRR[Iz,iD,uPos] = NV[1,ID,Iz,IF] * U[Iz,1,ind,uPos] +
        NV[2,ID,Iz,IF] * U[Iz,1,ind,vPos] +
        NV[3,ID,Iz,IF] * U[Iz,1,ind,wPos]
      VRR[Iz,iD,vPos] = T1V[1,ID,Iz,IF] * U[Iz,1,ind,uPos] +
        T1V[2,ID,Iz,IF] * U[Iz,1,ind,vPos] +
        T1V[3,ID,Iz,IF] * U[Iz,1,ind,wPos]
      VRR[Iz,iD,wPos] = T2V[1,ID,Iz,IF] * U[Iz,1,ind,uPos] +
        T2V[2,ID,Iz,IF] * U[Iz,1,ind,vPos] +
        T2V[3,ID,Iz,IF] * U[Iz,1,ind,wPos]
      VRR[Iz,iD,ThPos] = U[Iz,1,ind,ThPos]
    else  
      for iAux = 1 : NAUX  
        AuxR[Iz,iD,iAux] = Aux[Iz-1,M,ind,iAux]
      end  
      VRR[Iz,iD,RhoPos] = U[Iz-1,M,ind,RhoPos]
      VRR[Iz,iD,uPos] = -(NV[1,ID,Iz,IF] * U[Iz-1,M,ind,uPos] +
        NV[2,ID,Iz,IF] * U[Iz-1,M,ind,vPos] +
        NV[3,ID,Iz,IF] * U[Iz-1,M,ind,wPos])
      VRR[Iz,iD,vPos] = T1V[1,ID,Iz,IF] * U[Iz-1,M,ind,uPos] +
        T1V[2,ID,Iz,IF] * U[Iz-1,M,ind,vPos] +
        T1V[3,ID,Iz,IF] * U[Iz-1,M,ind,wPos]
      VRR[Iz,iD,wPos] = T2V[1,ID,Iz,IF] * U[Iz-1,M,ind,uPos] +
        T2V[2,ID,Iz,IF] * U[Iz-1,M,ind,vPos] +
        T2V[3,ID,Iz,IF] * U[Iz-1,M,ind,wPos]
      VRR[Iz,iD,ThPos] = U[Iz-1,M,ind,ThPos]
    end

    @views RiemannSolver!(FLoc,VLL[Iz,iD,:],VRR[Iz,iD,:],
      AuxL[Iz,iD,:],AuxR[Iz,iD,:])

    FE[RhoPos] =  FLoc[RhoPos]
    FE[uPos] = NV[1,ID,Iz,IF] * FLoc[uPos] +
      T1V[1,ID,Iz,IF] * FLoc[vPos] +
      T2V[1,ID,Iz,IF] * FLoc[wPos] 
    FE[vPos] = NV[2,ID,Iz,IF] * FLoc[uPos] +
      T1V[2,ID,Iz,IF] * FLoc[vPos] +
      T2V[2,ID,Iz,IF] * FLoc[wPos] 
    FE[wPos] = NV[3,ID,Iz,IF] * FLoc[uPos] +
      T1V[3,ID,Iz,IF] * FLoc[vPos] +
      T2V[3,ID,Iz,IF] * FLoc[wPos] 
    FE[ThPos] =  FLoc[ThPos]

    Surf = VolSurfV[ID,Iz,IF] / w[1]  
    FE[RhoPos] *= Surf
    FE[uPos] *= Surf
    FE[vPos] *= Surf
    FE[wPos] *= Surf
    FE[ThPos] *= Surf
    if Iz > 1 
      @atomic :monotonic F[Iz-1,M,ind,RhoPos] -= FE[RhoPos]
      @atomic :monotonic F[Iz-1,M,ind,uPos] -= FE[uPos]
      @atomic :monotonic F[Iz-1,M,ind,vPos] -= FE[vPos]
      @atomic :monotonic F[Iz-1,M,ind,wPos] -= FE[wPos]
      @atomic :monotonic F[Iz-1,M,ind,ThPos] -= FE[ThPos]
    end  
    if Iz < Nz
      @atomic :monotonic F[Iz,1,ind,RhoPos] += FE[RhoPos]
      @atomic :monotonic F[Iz,1,ind,uPos] += FE[uPos]
      @atomic :monotonic F[Iz,1,ind,vPos] += FE[vPos]
      @atomic :monotonic F[Iz,1,ind,wPos] += FE[wPos]
      @atomic :monotonic F[Iz,1,ind,ThPos] += FE[ThPos]
    end  
  end  
end

@kernel inbounds = true function RiemanNonLinH3Kernel!(RiemannSolver!,F,@Const(U),@Const(Aux),@Const(GlobE),
  @Const(EF),@Const(FTE),@Const(NH),@Const(T1H),@Const(T2H),@Const(VolSurfH),
  @Const(w), NF, ::Val{NUMV}, ::Val{NAUX}) where {NUMV, NAUX}

  _,_,iz, = @index(Local, NTuple)
  I,K,Iz,IE = @index(Global, NTuple)

  N = @uniform @groupsize()[1]
  M = @uniform @groupsize()[2]
  TilesDim = @uniform @groupsize()[3]

  Nz = @uniform @ndrange()[3]

  VLL = @localmem eltype(F) (N,M,TilesDim,NUMV)
  VRR = @localmem eltype(F) (N,M,TilesDim,NUMV)
  AuxL = @localmem eltype(F) (N,M,TilesDim,NAUX)
  AuxR = @localmem eltype(F) (N,M,TilesDim,NAUX)
  FLoc = @private eltype(F) (NUMV,)
  FE = @private eltype(F) (NUMV,)

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
    for iAux = 1 : NAUX
      AuxL[I,K,iz,iAux] = Aux[Iz,K,indL,iAux]
      AuxR[I,K,iz,iAux] = Aux[Iz,K,indR,iAux]
    end  
    VLL[I,K,iz,RhoPos] = U[Iz,K,indL,RhoPos]
    VLL[I,K,iz,uPos] = NH[1,K,I,Iz,IE] * U[Iz,K,indL,uPos] +
      NH[2,K,I,Iz,IE] * U[Iz,K,indL,vPos] +
      NH[3,K,I,Iz,IE] * U[Iz,K,indL,wPos]
    VLL[I,K,iz,vPos] = T1H[1,K,I,Iz,IE] * U[Iz,K,indL,uPos] +
      T1H[2,K,I,Iz,IE] * U[Iz,K,indL,vPos] +
      T1H[3,K,I,Iz,IE] * U[Iz,K,indL,wPos]
    VLL[I,K,iz,wPos] = T2H[1,K,I,Iz,IE] * U[Iz,K,indL,uPos] +
      T2H[2,K,I,Iz,IE] * U[Iz,K,indL,vPos] +
      T2H[3,K,I,Iz,IE] * U[Iz,K,indL,wPos]
    VLL[I,K,iz,ThPos] = U[Iz,K,indL,ThPos]
    VRR[I,K,iz,RhoPos] = U[Iz,K,indR,RhoPos]
    VRR[I,K,iz,uPos] = NH[1,K,I,Iz,IE] * U[Iz,K,indR,uPos] +
      NH[2,K,I,Iz,IE] * U[Iz,K,indR,vPos] +
      NH[3,K,I,Iz,IE] * U[Iz,K,indR,wPos]
    VRR[I,K,iz,vPos] = T1H[1,K,I,Iz,IE] * U[Iz,K,indR,uPos] +
      T1H[2,K,I,Iz,IE] * U[Iz,K,indR,vPos] +
      T1H[3,K,I,Iz,IE] * U[Iz,K,indR,wPos]
    VRR[I,K,iz,wPos] = T2H[1,K,I,Iz,IE] * U[Iz,K,indR,uPos] +
      T2H[2,K,I,Iz,IE] * U[Iz,K,indR,vPos] +
      T2H[3,K,I,Iz,IE] * U[Iz,K,indR,wPos]
    VRR[I,K,iz,ThPos] = U[Iz,K,indL,ThPos]
    @views RiemannSolver!(FLoc,VLL[I,K,iz,:],VRR[I,K,iz,:],
      AuxL[I,K,iz,:],AuxR[I,K,iz,:])

    FE[RhoPos] =  FLoc[RhoPos]
    FE[uPos] =  NH[1,K,I,Iz,IE] * FLoc[uPos] +
      T1H[1,K,I,Iz,IE] * FLoc[vPos] + T2H[1,K,I,Iz,IE] * FLoc[wPos]
    FE[vPos] =  NH[2,K,I,Iz,IE] * FLoc[uPos] +
      T1H[2,K,I,Iz,IE] * FLoc[vPos] + T2H[2,K,I,Iz,IE] * FLoc[wPos]
    FE[wPos] =  NH[3,K,I,Iz,IE] * FLoc[uPos] +
      T1H[3,K,I,Iz,IE] * FLoc[vPos] + T2H[3,K,I,Iz,IE] * FLoc[wPos]
    FE[ThPos] =  FLoc[ThPos]
    Surf = VolSurfH[K,I,Iz,IE] / w[1]  
    FE[RhoPos] *= Surf
    FE[uPos] *= Surf
    FE[vPos] *= Surf
    FE[wPos] *= Surf
    FE[ThPos] *= Surf
    if iFL <= NF
      @atomic :monotonic F[Iz,K,indL,RhoPos] -= FE[RhoPos]
      @atomic :monotonic F[Iz,K,indL,uPos] -= FE[uPos]
      @atomic :monotonic F[Iz,K,indL,vPos] -= FE[vPos]
      @atomic :monotonic F[Iz,K,indL,wPos] -= FE[wPos]
      @atomic :monotonic F[Iz,K,indL,ThPos] -= FE[ThPos]
    end  
    if iFR <= NF
      @atomic :monotonic F[Iz,K,indR,RhoPos] += FE[RhoPos]
      @atomic :monotonic F[Iz,K,indR,uPos] += FE[uPos]
      @atomic :monotonic F[Iz,K,indR,vPos] += FE[vPos]
      @atomic :monotonic F[Iz,K,indR,wPos] += FE[wPos]
      @atomic :monotonic F[Iz,K,indR,ThPos] += FE[ThPos]
    end  
  end  
end

@kernel inbounds = true function VSp2VCart3Kernel!(VCart,@Const(VSp),@Const(Rotate),@Const(Glob))

  _,_,iD,  = @index(Local, NTuple)
  Iz,K,ID,IF = @index(Global, NTuple)

  NQ = @uniform @ndrange()[3]

  if ID <= NQ
    ind = Glob[ID,IF]
    VCart[Iz,K,ind,1] = VSp[Iz,K,ind,1]
    VCart[Iz,K,ind,2] = Rotate[1,1,K,ID,Iz,IF] * VSp[Iz,K,ind,2] +
      Rotate[2,1,K,ID,Iz,IF] * VSp[Iz,K,ind,3] +
      Rotate[3,1,K,ID,Iz,IF] * VSp[Iz,K,ind,4]
    VCart[Iz,K,ind,3] = Rotate[1,2,1,ID,Iz,IF] * VSp[Iz,K,ind,2] +
      Rotate[2,2,K,ID,Iz,IF] * VSp[Iz,K,ind,3] +
      Rotate[3,2,K,ID,Iz,IF] * VSp[Iz,K,ind,4] 
    VCart[Iz,K,ind,4] = Rotate[1,3,K,ID,Iz,IF] * VSp[Iz,K,ind,2] +
      Rotate[2,3,K,ID,Iz,IF] * VSp[Iz,K,ind,3] +
      Rotate[3,3,K,ID,Iz,IF] * VSp[Iz,K,ind,4]
    VCart[Iz,K,ind,5] = VSp[Iz,K,ind,5]
  end  
end


@kernel inbounds = true function VCart2VSp3Kernel!(VSp,@Const(VCart),@Const(Rotate),@Const(J),@Const(Glob))

  _,_,iD,  = @index(Local, NTuple)
  Iz,K,ID,IF = @index(Global, NTuple)

  NQ = @uniform @ndrange()[3]

  if ID <= NQ
    ind = Glob[ID,IF]
    VSp[Iz,K,ind,1] = VCart[Iz,K,ind,1] / J[ID,K,Iz,IF]
    VSp[Iz,K,ind,2] = (Rotate[1,1,K,ID,Iz,IF] * VCart[Iz,K,ind,2] +
      Rotate[1,2,K,ID,Iz,IF] * VCart[Iz,K,ind,3] +
      Rotate[1,3,K,ID,Iz,IF] * VCart[Iz,K,ind,4]) / J[ID,K,Iz,IF]
    VSp[Iz,K,ind,3] = (Rotate[2,1,K,ID,Iz,IF] * VCart[Iz,K,ind,2] +
      Rotate[2,2,K,ID,Iz,IF] * VCart[Iz,K,ind,3] +
      Rotate[2,3,K,ID,Iz,IF] * VCart[Iz,K,ind,4]) / J[ID,K,Iz,IF]
    VSp[Iz,K,ind,4] = (Rotate[3,1,K,ID,Iz,IF] * VCart[Iz,K,ind,2] +
      Rotate[3,2,K,ID,Iz,IF] * VCart[Iz,K,ind,3] +
      Rotate[3,3,K,ID,Iz,IF] * VCart[Iz,K,ind,4]) / J[ID,K,Iz,IF]
    VSp[Iz,K,ind,5] = VCart[Iz,K,ind,5] / J[ID,K,Iz,IF]
  end  
end
