function VSp2VCart!(VCart,VSp,Rotate)
  @. VCart[:,:,:,1] = VSp[:,:,:,1]
  NF = size(Rotate,6)
  nQuad = size(Rotate,4)
  iD = 0
  @inbounds for iF = 1 : NF
    @inbounds for iQ = 1 : nQuad  
      iD += 1  
      VCart[1,1,iD,2] = Rotate[1,1,1,iQ,1,iF] * VSp[1,1,iD,2] +
        Rotate[2,1,1,iQ,1,iF] * VSp[1,1,iD,3]
      VCart[1,1,iD,3] = Rotate[1,2,1,iQ,1,iF] * VSp[1,1,iD,2] +
        Rotate[2,2,1,iQ,1,iF] * VSp[1,1,iD,3]
      VCart[1,1,iD,4] = Rotate[1,3,1,iQ,1,iF] * VSp[1,1,iD,2] +
        Rotate[2,3,1,iQ,1,iF] * VSp[1,1,iD,3]
    end
  end  
end

function VCart2VSp!(VSp,VCart,Rotate)
  @. VSp[:,:,:,1] = VCart[:,:,:,1]
  NF = size(Rotate,6)
  nQuad = size(Rotate,4)
  iD = 0
  @inbounds for iF = 1 : NF
    @inbounds for iQ = 1 : nQuad
      iD += 1
      VSp[1,1,iD,2] = Rotate[1,1,1,iQ,1,iF] * VCart[1,1,iD,2] +
        Rotate[1,2,1,iQ,1,iF] * VCart[1,1,iD,3] +
        Rotate[1,3,1,iQ,1,iF] * VCart[1,1,iD,4]
      VSp[1,1,iD,3] = Rotate[2,1,1,iQ,1,iF] * VCart[1,1,iD,2] +
        Rotate[2,2,1,iQ,1,iF] * VCart[1,1,iD,3] +
        Rotate[2,3,1,iQ,1,iF] * VCart[1,1,iD,4]
    end
  end
end

function RiemanNonLinKernel(F,U,DG,Metric,Grid,Phys)

  OrdPoly = DG.OrdPoly
  Glob = DG.Glob

  NH = Metric.NH
  T1H = Metric.T1H
  T2H = Metric.T2H
  VolSurfH = Metric.VolSurfH
  JJ = Metric.J

  hPos = 1
  uPos = 2
  vPos = 3
  wPos = 4

  IndE = zeros(Int,4,OrdPoly+1)
  @inbounds for i = 1 : OrdPoly + 1
    IndE[1,i] = i                               #  1  2  3  4  5
    IndE[2,i] = i * (OrdPoly + 1)               #  5 10 15 20 25
    IndE[3,i] = i  + OrdPoly * (OrdPoly + 1)    # 21 22 23 24 25 
    IndE[4,i] = 1  + (i-1) * (OrdPoly + 1)      #  1  6 11 16 21 
  end  

  VLL = zeros(4)
  VRR = zeros(4)
  FLoc = zeros(4)
  FE = zeros(OrdPoly+1,4)
  VZ = zeros(4)
  VZ[1] = -1.0
  VZ[2] = -1.0
  VZ[3] = 1.0
  VZ[4] = 1.0

  @inbounds for iE = 1 : Grid.NumEdges
    iFL = Grid.Edges[iE].F[1]  
    iFR = Grid.Edges[iE].F[2]  
    EL = Grid.Edges[iE].FE[1]
    ER = Grid.Edges[iE].FE[2]
    VZL = VZ[EL]

    @inbounds for i = 1 : OrdPoly + 1
      indL = Glob[IndE[EL,i],iFL]  
      indR = Glob[IndE[ER,i],iFR]  
      VLL[hPos] = U[1,1,indL,hPos]
      VLL[uPos] = NH[1,1,i,1,iE] * U[1,1,indL,uPos] +  
        NH[2,1,i,1,iE] * U[1,1,indL,vPos] +  
        NH[3,1,i,1,iE] * U[1,1,indL,wPos]  
      VLL[vPos] = T1H[1,1,i,1,iE] * U[1,1,indL,uPos] +  
        T1H[2,1,i,1,iE] * U[1,1,indL,vPos] +  
        T1H[3,1,i,1,iE] * U[1,1,indL,wPos]  
      VLL[wPos] = T2H[1,1,i,1,iE] * U[1,1,indL,uPos] +  
        T2H[2,1,i,1,iE] * U[1,1,indL,vPos] +  
        T2H[3,1,i,1,iE] * U[1,1,indL,wPos]  
      VRR[hPos] = U[1,1,indR,hPos]
      VRR[uPos] = NH[1,1,i,1,iE] * U[1,1,indR,uPos] +  
        NH[2,1,i,1,iE] * U[1,1,indR,vPos] +  
        NH[3,1,i,1,iE] * U[1,1,indR,wPos]  
      VRR[vPos] = T1H[1,1,i,1,iE] * U[1,1,indR,uPos] +  
        T1H[2,1,i,1,iE] * U[1,1,indR,vPos] +  
        T1H[3,1,i,1,iE] * U[1,1,indR,wPos]  
      VRR[wPos] = T2H[1,1,i,1,iE] * U[1,1,indR,uPos] +  
        T2H[2,1,i,1,iE] * U[1,1,indR,vPos] +  
        T2H[3,1,i,1,iE] * U[1,1,indR,wPos]  
      RiemannByLMARSNonLin!(FLoc,VLL,VRR,Phys)  
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
      F[1,1,indL,hPos] += VZL * FE[i,hPos] / DG.w[1]
      F[1,1,indL,uPos] += VZL * FE[i,uPos] / DG.w[1]
      F[1,1,indL,vPos] += VZL * FE[i,vPos] / DG.w[1]
      F[1,1,indL,wPos] += VZL * FE[i,wPos] / DG.w[1]
      F[1,1,indR,hPos] -= VZL * FE[i,hPos] / DG.w[1]
      F[1,1,indR,uPos] -= VZL * FE[i,uPos] / DG.w[1]
      F[1,1,indR,vPos] -= VZL * FE[i,vPos] / DG.w[1]
      F[1,1,indR,wPos] -= VZL * FE[i,wPos] / DG.w[1]
    end
  end
  iDG = 0
  @inbounds for iF = 1 : Grid.NumFaces
     @inbounds for iD = 1 : (OrdPoly + 1) * (OrdPoly + 1) 
       iDG += 1  
       @views F[1,1,iDG,:] /= JJ[iD,1,1,iF]
     end
  end   
end
