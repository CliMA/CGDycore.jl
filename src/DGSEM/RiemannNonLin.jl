function VSp2VCart!(VCart,VSp,Rotate)
  @. VCart[:,:,:,1] = VSp[:,:,:,1]
  NF = size(Rotate,6)
  nQuad = size(Rotate,4)
  iD = 0
  for iF = 1 : NF
    for iQ = 1 : nQuad  
      iD += 1  
      VCart[1,1,iD,2] = Rotate[1,1,1,iQ,1,iF] * VSp[1,1,iD,2] +
        Rotate[2,1,1,iQ,1,iF] * VSp[1,1,iD,3]
      VCart[1,1,iD,3] = Rotate[1,2,1,iQ,1,iF] * VSp[1,1,iD,2] +
        Rotate[2,2,1,iQ,1,iF] * VSp[1,1,iD,3]
      VCart[1,1,iD,4] = Rotate[1,2,1,iQ,1,iF] * VSp[1,1,iD,2] +
        Rotate[2,3,1,iQ,1,iF] * VSp[1,1,iD,3]
    end
  end  
end

#function VSp=VCart2VSp(VCart,rot)
#NX=size(VCart,2);
#NV=size(VCart,1);
#NF=size(VCart,4);
#VSp=zeros(NV-1,NX,NX,NF);
#VSp(1,:,:,:)=VCart(1,:,:,:);
#for iF=1:NF
#  VSp(2,:,:,iF)=rot(:,:,1,1,iF).*reshape(VCart(2,:,:,iF),NX,NX)...
#    +rot(:,:,1,2,iF).*reshape(VCart(3,:,:,iF),NX,NX)...
#    +rot(:,:,1,3,iF).*reshape(VCart(4,:,:,iF),NX,NX);
#  VSp(3,:,:,iF)=rot(:,:,2,1,iF).*reshape(VCart(2,:,:,iF),NX,NX)...
#    +rot(:,:,2,2,iF).*reshape(VCart(3,:,:,iF),NX,NX)...
#    +rot(:,:,2,3,iF).*reshape(VCart(4,:,:,iF),NX,NX);
#end
#end
function RiemanNonLinKernel(F,U,DG,Metric,Grid)

  OrdPoly = DG.OrdPoly
  Glob = DG.Glob

  NH = Metric.NH
  T1H = Metric.T1H
  T2H = Metric.T2H
  VolSurfH = Metric.VolSurfH

  hPos = 1
  uPos = 2
  vPos = 3
  wPos = 4

  IndE = zeros(Int,4,OrdPoly+1)
  for i = 1 : OrdPoly + 1
    IndE[1,i] = i
    IndE[2,i] = i * (OrdPoly + 1)
    IndE[3,i] = i  + OrdPoly * (OrdPoly + 1)
    IndE[4,i] = 1  + (i-1) * (OrdPoly + 1)
  end  

  VLL = zeros(4)
  VRR = zeros(4)
  FLoc = zeros(4)
  FE = zeros(OrdPoly+1,4)

  for iE = 1 : Grid.NumEdges
    iFL = Grid.Edges[iE].F[1]  
    iFR = Grid.Edges[iE].F[2]  
    EL = Grid.Edges[iE].FE[1]
    ER = Grid.Edges[iE].FE[2]

    for i = 1 : OrdPoly + 1
      indL = Glob[IndE[EL,i],iFL]  
      indR = Glob[IndE[ER,i],iFR]  
      VLL[hPos] = U[1,1,indL,hPos]
      VLL[uPos] = NH[1,1,i,1,iE] * U[1,1,indL,uPos] +  
        NH[2,1,i,1,iE] * U[1,1,indL,vPos]  
        NH[3,1,i,1,iE] * U[1,1,indL,wPos]  
      VLL[vPos] = T1H[1,1,i,1,iE] * U[1,1,indL,uPos] +  
        T1H[2,1,i,1,iE] * U[1,1,indL,vPos]  
        T1H[3,1,i,1,iE] * U[1,1,indL,wPos]  
      VLL[wPos] = T2H[1,1,i,1,iE] * U[1,1,indL,uPos] +  
        T2H[2,1,i,1,iE] * U[1,1,indL,vPos]  
        T2H[3,1,i,1,iE] * U[1,1,indL,wPos]  
      VRR[hPos] = U[1,1,indR,hPos]
      VRR[uPos] = NH[1,1,i,1,iE] * U[1,1,indR,uPos] +  
        NH[2,1,i,1,iE] * U[1,1,indR,vPos]  
        NH[3,1,i,1,iE] * U[1,1,indR,wPos]  
      VRR[vPos] = T1H[1,1,i,1,iE] * U[1,1,indR,uPos] +  
        T1H[2,1,i,1,iE] * U[1,1,indR,vPos]  
        T1H[3,1,i,1,iE] * U[1,1,indR,wPos]  
      VRR[wPos] = T2H[1,1,i,1,iE] * U[1,1,indR,uPos] +  
        T2H[2,1,i,1,iE] * U[1,1,indR,vPos]  
        T2H[3,1,i,1,iE] * U[1,1,indR,wPos]  
      RiemannByLMARSNonLin!(FLoc,VLL,VRR)  
      FE[i,hPos] =  FLoc[hPos]
      FE[i,uPos] =  NH[1,1,i,1,iE] * FLoc[uPos] + 
        T1H[1,1,i,1,iE] * FLoc[vPos] + T2H[1,1,i,1,iE] * FLoc[wPos]
      FE[i,vPos] =  NH[2,1,i,1,iE] * FLoc[uPos] + 
        T1H[2,1,i,1,iE] * FLoc[vPos] + T2H[2,1,i,1,iE] * FLoc[wPos]
      FE[i,wPos] =  NH[3,1,i,1,iE] * FLoc[uPos] + 
        T1H[3,1,i,1,iE] * FLoc[vPos] + T2H[3,1,i,1,iE] * FLoc[wPos]
      FE[i,hPos] *= VolSurfH[1,i,1,iE]  
      FE[i,uPos] *= VolSurfH[1,i,1,iE]  
      FE[i,vPos] *= VolSurfH[1,i,1,iE]  
      FE[i,wPos] *= VolSurfH[1,i,1,iE]  
      F[1,1,indL,hPos] += - FE[i,hPos] / DG.w[1]
      F[1,1,indL,uPos] += - FE[i,uPos] / DG.w[1]
      F[1,1,indL,vPos] += - FE[i,vPos] / DG.w[1]
      F[1,1,indL,wPos] += - FE[i,wPos] / DG.w[1]
      F[1,1,indR,hPos] += FE[i,hPos] / DG.w[1]
      F[1,1,indR,uPos] += FE[i,uPos] / DG.w[1]
      F[1,1,indR,vPos] += FE[i,vPos] / DG.w[1]
      F[1,1,indR,wPos] += FE[i,wPos] / DG.w[1]
    end
  end
end
