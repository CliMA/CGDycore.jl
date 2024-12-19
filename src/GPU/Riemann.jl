@kernel inbounds = true function RiemannNonLinH!(F,@Const(U),@Const(N),@Const(T1),@Const(T2),
  @Const(VolSurf),@Const(EF),@Const(FE),@Const(Glob))

# FE (2,Grid.NumEdges)
# Group M,N,Nz,

  K, I, iz   = @index(Local, NTuple)
  _,_,Iz,IE = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  M = @uniform @groupsize()[1]
  N = @uniform @groupsize()[2]
  Nz = @uniform @ndrange()[3]
  NE = @uniform @ndrange()[4]

  UColL = @localmem eltype(F) (M,N,ColumnTilesDim,5)
  UColR = @localmem eltype(F) (M,N,ColumnTilesDim,5)
  FCol = @localmem eltype(F) (M,N,ColumnTilesDim,5)

  if Iz <= Nz && IE <= NE
  SL = FE[1,IE]
  FSL = EF[1,IE]
  if SL == 1
    iDL = I   
  elseif SL == 2
    iDL = N + (I - 1) * N
  elseif SL == 3
    iDL = N * (N - 1) + I
  elseif SL == 4
    iDL = 1 + (I - 1) * N
  end  
  indL = Glob[iDL,FSL]

  SR = FE[2,IE]
  FSR = EF[2,IE]
  if SR == 1
    iDR = I
  elseif SR == 2
    iDR = N + (I - 1) * N
  elseif SR == 3
    iDR = N * (N - 1) + I
  elseif SR == 4
    iDR = 1 + (I - 1) * N
  end
  indR = Glob[iDR,FSR]

  UColL[K,I,iz,1] = U[Iz,K,indL,1]  
  u1 = U[Iz,K,indL,2] 
  u2 = U[Iz,K,indL,3] 
  u3 = U[Iz,K,indL,4] 
  UColL[K,I,iz,2] = N[1,K,I,Iz,IE] * u1 +
    N[2,K,I,Iz,IE] * u2 + N[2,K,I,Iz,IE] * u3
  UColL[K,I,iz,3] = T1[1,K,I,Iz,IE] * u1 +
    T1[2,K,I,Iz,IE] * u2 + T1[2,K,I,Iz,IE] * u3
  UColL[K,I,iz,4] = T2[1,K,I,Iz,IE] * u1 +
    T2[2,K,I,Iz,IE] * u2 + T2[2,K,I,Iz,IE] * u3
  UColL[K,I,iz,5] = U[Iz,K,indL,5]  

  UColR[K,I,iz,1] = U[Iz,K,indR,1]
  u1 = U[Iz,K,indR,2]
  u2 = U[Iz,K,indR,3]
  u3 = U[Iz,K,indR,4]
  UColR[K,I,iz,2] = N[1,K,I,Iz,IE] * u1 +
    N[2,K,I,Iz,IE] * u2 + N[2,K,I,Iz,IE] * u3
  UColR[K,I,iz,3] = T1[1,K,I,Iz,IE] * u1 +
    T1[2,K,I,Iz,IE] * u2 + T1[2,K,I,Iz,IE] * u3
  UColR[K,I,iz,4] = T2[1,K,I,Iz,IE] * u1 +
    T2[2,K,I,Iz,IE] * u2 + T2[2,K,I,Iz,IE] * u3
  UColR[K,I,iz,5] = U[Iz,K,indR,5]

  RiemannByLMARSNonLin!(FCol[K,I,iz,:],UColL[K,I,iz,:],UColR[K,I,iz,:])
  F1 = FCol[K,I,iz,2]
  F2 = FCol[K,I,iz,3]
  F3 = FCol[K,I,iz,4]
  FCol[K,I,iz,2] = N[1,K,I,Iz,IE] * F1 +
    T1[1,K,I,Iz,IE] * F2 + T2[1,K,I,Iz,IE] * F3
  FCol[K,I,iz,3] = N[2,K,I,Iz,IE] * F1 +
    T1[2,K,I,Iz,IE] * F2 + T2[2,K,I,Iz,IE] * F3
  FCol[K,I,iz,4] = N[3,K,I,Iz,IE] * F1 +
    T1[3,K,I,Iz,IE] * F2 + T2[3,K,I,Iz,IE] * F3
  @views @. FCol[K,I,iz,:] *= VolSurf[K,I,Iz,IE]  
  @views @. F[Iz,K,indL,:] -= FCol[K,I,iz,:]
  @views @.  F[Iz,K,indR,:] += FCol[K,I,iz,:]
  end
end

@inline function RiemannByLMARSNonLin!(F,UL,UR)
  cS = 360.0 
  pL = 1.e5 #PresShB(UL,Param)
  pR = 1.e5 #PresShB(UR,Param)
  RhoM = 0.5 * (UL[1] + UR[1])
  uL = UL[2] / UL[1]
  uR = UR[2] / UR[1]

  pM = 0.5 * (pL + pR) - 0.5 * cS * RhoM * (uR - uL)
  uM = 0.5 * (uR + uL) - 1.0 / (2.0 * cS) * (pR - pL) / RhoM

  if uM > 0.0
    F[1]= uM * UL[1] 
    F[2]= pM + uM * UL[2] 
    F[3]= uM * UL[3] 
    F[4]= uM * UL[4] 
    F[5]= uM * UL[5] 
  else  
    F[1]= uM * UR[1] 
    F[2]= pM + uM * UR[2] 
    F[3]= uM * UR[3] 
    F[4]= uM * UR[4] 
    F[5]= uM * UR[5] 
  end  
end
