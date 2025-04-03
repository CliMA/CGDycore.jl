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
