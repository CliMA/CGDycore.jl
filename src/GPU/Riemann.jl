@kernel inbounds = true function RiemannNonLinH!(F,@Const(U),@Const(p),@Const(NN),@Const(T1),@Const(T2),
  @Const(VolSurf),w,@Const(JJ),@Const(EF),@Const(FE),@Const(Glob))

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
  pColL = @localmem eltype(F) (M,N,ColumnTilesDim)
  UColR = @localmem eltype(F) (M,N,ColumnTilesDim,5)
  FCol = @localmem eltype(F) (M,N,ColumnTilesDim,5)
  pColR = @localmem eltype(F) (M,N,ColumnTilesDim)

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
    UColL[K,I,iz,2] = NN[1,K,I,Iz,IE] * u1 +
      NN[2,K,I,Iz,IE] * u2 + NN[3,K,I,Iz,IE] * u3
    UColL[K,I,iz,3] = T1[1,K,I,Iz,IE] * u1 +
      T1[2,K,I,Iz,IE] * u2 + T1[3,K,I,Iz,IE] * u3
    UColL[K,I,iz,4] = T2[1,K,I,Iz,IE] * u1 +
      T2[2,K,I,Iz,IE] * u2 + T2[3,K,I,Iz,IE] * u3
    UColL[K,I,iz,5] = U[Iz,K,indL,5]  
    pColL[K,I,iz] = p[Iz,K,indL]  

    UColR[K,I,iz,1] = U[Iz,K,indR,1]
    u1 = U[Iz,K,indR,2]
    u2 = U[Iz,K,indR,3]
    u3 = U[Iz,K,indR,4]
    UColR[K,I,iz,2] = NN[1,K,I,Iz,IE] * u1 +
      NN[2,K,I,Iz,IE] * u2 + NN[3,K,I,Iz,IE] * u3
    UColR[K,I,iz,3] = T1[1,K,I,Iz,IE] * u1 +
      T1[2,K,I,Iz,IE] * u2 + T1[3,K,I,Iz,IE] * u3
    UColR[K,I,iz,4] = T2[1,K,I,Iz,IE] * u1 +
      T2[2,K,I,Iz,IE] * u2 + T2[3,K,I,Iz,IE] * u3
    UColR[K,I,iz,5] = U[Iz,K,indR,5]
    pColR[K,I,iz] = p[Iz,K,indR]  

    @views RiemannByLMARSNonLin!(FCol[K,I,iz,:],UColL[K,I,iz,:],pColL[K,I,iz],
      UColR[K,I,iz,:],pColR[K,I,iz])
    F1 = FCol[K,I,iz,2]
    F2 = FCol[K,I,iz,3]
    F3 = FCol[K,I,iz,4]
    FCol[K,I,iz,2] = NN[1,K,I,Iz,IE] * F1 +
      T1[1,K,I,Iz,IE] * F2 + T2[1,K,I,Iz,IE] * F3
    FCol[K,I,iz,3] = NN[2,K,I,Iz,IE] * F1 +
      T1[2,K,I,Iz,IE] * F2 + T2[2,K,I,Iz,IE] * F3
    FCol[K,I,iz,4] = NN[3,K,I,Iz,IE] * F1 +
      T1[3,K,I,Iz,IE] * F2 + T2[3,K,I,Iz,IE] * F3
    @views @. FCol[K,I,iz,:] *= VolSurf[K,I,Iz,IE] / w  
    for iv = 1 : 5
      if SL == 1 || SL == 4  
        @atomic F[Iz,K,indL,iv] += FCol[K,I,iz,iv] / JJ[iDL,K,Iz,FSL]
        @atomic F[Iz,K,indR,iv] += -FCol[K,I,iz,iv] / JJ[iDR,K,Iz,FSR]
      else  
        @atomic F[Iz,K,indL,iv] += -FCol[K,I,iz,iv] / JJ[iDL,K,Iz,FSL]
        @atomic F[Iz,K,indR,iv] += FCol[K,I,iz,iv] / JJ[iDR,K,Iz,FSR]
      end  
    end
  end
end

@kernel inbounds = true function RiemannNonLinV!(F,@Const(U),@Const(p),M,@Const(NN),@Const(T1),@Const(T2),
  @Const(VolSurf),w,@Const(JJ),@Const(Glob))

# FE (2,Grid.NumEdges)
# Group M,N,Nz,

  iD, iz   = @index(Local, NTuple)
  _,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform @ndrange()[1]
  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  UColL = @localmem eltype(F) (N,ColumnTilesDim,5)
  pColL = @localmem eltype(F) (N,ColumnTilesDim)
  UColR = @localmem eltype(F) (N,ColumnTilesDim,5)
  pColR = @localmem eltype(F) (N,ColumnTilesDim)
  FCol = @localmem eltype(F) (N,ColumnTilesDim,5)

  if Iz <= Nz && IF <= NF

    ind = Glob[iD,IF]

    if Iz > 1
      UColL[iD,iz,1] = U[Iz-1,M,ind,1]  
      u1 = U[Iz-1,M,ind,2] 
      u2 = U[Iz-1,M,ind,3] 
      u3 = U[Iz-1,M,ind,4] 
      UColL[iD,iz,2] = NN[1,iD,Iz,IF] * u1 +
      NN[2,iD,Iz,IF] * u2 + NN[3,iD,Iz,IF] * u3
      UColL[iD,iz,3] = T1[1,iD,Iz,IF] * u1 +
        T1[2,iD,Iz,IF] * u2 + T1[3,iD,Iz,IF] * u3
      UColL[iD,iz,4] = T2[1,iD,Iz,IF] * u1 +
        T2[2,iD,Iz,IF] * u2 + T2[3,iD,Iz,IF] * u3
      UColL[iD,iz,5] = U[Iz-1,M,ind,5]  
      pColL[iD,iz] = p[Iz-1,M,ind]  
    end  
    if Iz < Nz 
      UColR[iD,iz,1] = U[Iz,1,ind,1]  
      u1 = U[Iz,1,ind,2] 
      u2 = U[Iz,1,ind,3] 
      u3 = U[Iz,1,ind,4] 
      UColR[iD,iz,2] = NN[1,iD,Iz,IF] * u1 +
        NN[2,iD,Iz,IF] * u2 + NN[3,iD,Iz,IF] * u3
      UColR[iD,iz,3] = T1[1,iD,Iz,IF] * u1 +
        T1[2,iD,Iz,IF] * u2 + T1[3,iD,Iz,IF] * u3
      UColR[iD,iz,4] = T2[1,iD,Iz,IF] * u1 +
        T2[2,iD,Iz,IF] * u2 + T2[3,iD,Iz,IF] * u3
      UColR[iD,iz,5] = U[Iz,1,ind,5]  
      pColR[iD,iz] = p[Iz,1,ind]  
    end  
    if Iz == 1
      UColL[iD,iz,1] = UColR[iD,iz,1]
      UColL[iD,iz,2] = -UColR[iD,iz,2]
      UColL[iD,iz,3] = UColR[iD,iz,3]
      UColL[iD,iz,4] = UColR[iD,iz,4]
      UColL[iD,iz,5] = UColR[iD,iz,5]
      pColL[iD,iz] = pColR[iD,iz]
    end  
    if Iz == Nz
      UColR[iD,iz,1] = UColL[iD,iz,1]
      UColR[iD,iz,2] = -UColL[iD,iz,2]
      UColR[iD,iz,3] = UColL[iD,iz,3]
      UColR[iD,iz,4] = UColL[iD,iz,4]
      UColR[iD,iz,5] = UColL[iD,iz,5]
      pColR[iD,iz] = pColL[iD,iz]
    end  

    @views RiemannByLMARSNonLin!(FCol[iD,iz,:],UColL[iD,iz,:],pColL[iD,iz],
      UColR[iD,iz,:],pColR[iD,iz])
    F1 = FCol[iD,iz,2]
    F2 = FCol[iD,iz,3]
    F3 = FCol[iD,iz,4]
    FCol[iD,iz,2] = NN[1,iD,Iz,IF] * F1 +
      T1[1,iD,Iz,IF] * F2 + T2[1,iD,Iz,IF] * F3
    FCol[iD,iz,3] = NN[2,iD,Iz,IF] * F1 +
      T1[2,iD,Iz,IF] * F2 + T2[2,iD,Iz,IF] * F3
    FCol[iD,iz,4] = NN[3,iD,Iz,IF] * F1 +
      T1[3,iD,Iz,IF] * F2 + T2[3,iD,Iz,IF] * F3
    @views @. FCol[iD,iz,:] *= VolSurf[iD,Iz,IF] / w  
    if Iz > 1
      for iv = 1 : 5
        @atomic F[Iz-1,M,ind,iv] -= FCol[iD,iz,iv] / JJ[iD,M,Iz-1,IF]
      end
    end 
    if Iz < Nz
      for iv = 1 : 5
        @atomic F[Iz,1,ind,iv] += FCol[iD,iz,iv] / JJ[iD,1,Iz,IF]
      end
    end  
  end
end

@inline function RiemannByLMARSNonLin!(F,UL,pL,UR,pR)
  cS = sqrt(10*1.e5) 
  cS = 1000.0
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
