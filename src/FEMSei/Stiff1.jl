function CurlVel1(q,FeT,u,uFe::HDivElement,QuadOrd,ElemType,Grid)
#
#
# int q*v dx = int Curl u * v dx = - int u * rot v dx + int_E (u*t)*q ds
#
#
  tBar = [1 0 1 0
          0 1  0 1]
  tBar = [1 0 -1 0
          0 1  0 -1]
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  if ElemType == Grids.Quad()
    uFRef  = zeros(uFe.Comp,uFe.DoF,NumQuadL,4)
    qRef  = zeros(FeT.DoF,NumQuadL,4)
    PointsE = zeros(2,NumQuadL,4)
    for iQ = 1 : NumQuadL
      PointsE[1,iQ,1] = PointsL[iQ]
      PointsE[2,iQ,1] = -1.0
      PointsE[1,iQ,2] = 1.0
      PointsE[2,iQ,2] = PointsL[iQ]
      PointsE[1,iQ,3] = PointsL[iQ]
      PointsE[2,iQ,3] = 1.0
      PointsE[1,iQ,4] = -1.0
      PointsE[2,iQ,4] = PointsL[iQ]
      for iDoF = 1 : uFe.DoF
        uFRef[1,iDoF,iQ,1] = uFe.phi[iDoF,1](PointsL[iQ],-1.0)  
        uFRef[2,iDoF,iQ,1] = uFe.phi[iDoF,2](PointsL[iQ],-1.0)  
        uFRef[1,iDoF,iQ,2] = uFe.phi[iDoF,1](1.0,PointsL[iQ])  
        uFRef[2,iDoF,iQ,2] = uFe.phi[iDoF,2](1.0,PointsL[iQ])  
        uFRef[1,iDoF,iQ,3] = uFe.phi[iDoF,1](PointsL[iQ],1.0)  
        uFRef[2,iDoF,iQ,3] = uFe.phi[iDoF,2](PointsL[iQ],1.0)  
        uFRef[1,iDoF,iQ,4] = uFe.phi[iDoF,1](-1.0,PointsL[iQ])
        uFRef[2,iDoF,iQ,4] = uFe.phi[iDoF,2](-1.0,PointsL[iQ])
      end  
      for iDoF = 1 : FeT.DoF
        qRef[iDoF,iQ,1] = FeT.phi[iDoF,1](PointsL[iQ],-1.0)  
        qRef[iDoF,iQ,2] = FeT.phi[iDoF,1](1.0,PointsL[iQ])  
        qRef[iDoF,iQ,3] = FeT.phi[iDoF,1](PointsL[iQ],1.0)  
        qRef[iDoF,iQ,4] = FeT.phi[iDoF,1](-1.0,PointsL[iQ])  
      end  
    end  
  end  
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  uLoc = zeros(uFe.DoF)
  uFLoc = zeros(2)
  qLoc = zeros(FeT.DoF)

  @. q = 0
  for iF = 1 : Grid.NumFaces
    for iDoF = 1 : uFe.DoF
      ind = uFe.Glob[iDoF,iF]  
      uLoc[iDoF] = u[ind] 
    end  
    @. qLoc = 0
    for iE = 1 : length(Grid.Faces[iF].E)
      for iQ = 1 : NumQuadL  
        @. uFLoc = 0
        for iDoF = 1 : uFe.DoF  
          uFLoc[1] += uFRef[1,iDoF,iQ,iE] * uLoc[iDoF]  
          uFLoc[2] += uFRef[2,iDoF,iQ,iE] * uLoc[iDoF]  
        end  
        Jacobi!(DF,detDF,pinvDF,X,ElemType,PointsE[1,iQ,iE],PointsE[2,iQ,iE],Grid.Faces[iF], Grid)
        t = ((DF[1,1] * uFLoc[1] + DF[1,2] * uFLoc[2]) *
          (DF[1,1] * tBar[1,iE] + DF[1,2] * tBar[2,iE]) +
          (DF[2,1] * uFLoc[1] + DF[2,2] * uFLoc[2]) *
          (DF[2,1] * tBar[1,iE] + DF[2,2] * tBar[2,iE]) +
          (DF[3,1] * uFLoc[1] + DF[3,2] * uFLoc[2]) *
          (DF[3,1] * tBar[1,iE] + DF[3,2] * tBar[2,iE])) / detDF[1] #* Grid.Faces[iF].Orientation
        for iDoF = 1 : FeT.DoF  
          qLoc[iDoF] += t * qRef[iDoF,iQ,iE]  
        end  
      end
    end  
    for iDoF = 1 : FeT.DoF
      ind = FeT.Glob[iDoF,iF]
      q[ind] += qLoc[iDoF]
    end   
  end    
  ldiv!(FeT.LUM,q)
end
function CurlVel(q,FeT,u,uFe::HDivElement,QuadOrd,ElemType,Grid)
#
#
# int q*v dx = int Curl u * v dx = - int u * rot v dx + int_E (u*t)*q ds
#
#
  tBar = [1 0 1 0
          0 1 0 1]
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  if ElemType == Grids.Quad()
    uFRef  = zeros(uFe.Comp,uFe.DoF,NumQuadL,4)
    qRef  = zeros(FeT.DoF,NumQuadL,4)
    PointsE = zeros(2,NumQuadL,4)
    for iQ = 1 : NumQuadL
      PointsE[1,iQ,1] = PointsL[iQ]
      PointsE[2,iQ,1] = -1.0
      PointsE[1,iQ,2] = 1.0
      PointsE[2,iQ,2] = PointsL[iQ]
      PointsE[1,iQ,3] = PointsL[iQ]
      PointsE[2,iQ,3] = 1.0
      PointsE[1,iQ,4] = -1.0
      PointsE[2,iQ,4] = PointsL[iQ]
      for iDoF = 1 : uFe.DoF
        uFRef[1,iDoF,iQ,1] = uFe.phi[iDoF,1](PointsL[iQ],-1.0)  
        uFRef[2,iDoF,iQ,1] = uFe.phi[iDoF,2](PointsL[iQ],-1.0)  
        uFRef[1,iDoF,iQ,2] = uFe.phi[iDoF,1](1.0,PointsL[iQ])  
        uFRef[2,iDoF,iQ,2] = uFe.phi[iDoF,2](1.0,PointsL[iQ])  
        uFRef[1,iDoF,iQ,3] = uFe.phi[iDoF,1](PointsL[iQ],1.0)  
        uFRef[2,iDoF,iQ,3] = uFe.phi[iDoF,2](PointsL[iQ],1.0)  
        uFRef[1,iDoF,iQ,4] = uFe.phi[iDoF,1](-1.0,PointsL[iQ])
        uFRef[2,iDoF,iQ,4] = uFe.phi[iDoF,2](-1.0,PointsL[iQ])
      end  
      for iDoF = 1 : FeT.DoF
        qRef[iDoF,iQ,1] = FeT.phi[iDoF,1](PointsL[iQ],-1.0)  
        qRef[iDoF,iQ,2] = FeT.phi[iDoF,1](1.0,PointsL[iQ])  
        qRef[iDoF,iQ,3] = FeT.phi[iDoF,1](PointsL[iQ],1.0)  
        qRef[iDoF,iQ,4] = FeT.phi[iDoF,1](-1.0,PointsL[iQ])  
      end  
    end  
  end  
  DFL = zeros(3,2)
  detDFL = zeros(1)
  pinvDFL = zeros(3,2)
  XL = zeros(3)
  DFR = zeros(3,2)
  detDFR = zeros(1)
  pinvDFR = zeros(3,2)
  XR = zeros(3)
  uLocL = zeros(uFe.DoF)
  uLocR = zeros(uFe.DoF)
  uFLocL = zeros(2)
  uFLocR = zeros(2)
  uuFLocL = zeros(3)
  uuFLocR = zeros(3)
  ttL = zeros(3)
  ttR = zeros(3)
  qLocL = zeros(FeT.DoF)
  qLocR = zeros(FeT.DoF)

  @. q = 0
  for iE = 1 : Grid.NumEdges
    Edge = Grid.Edges[iE]
    if length(Edge.F) > 1
      iFL = Edge.F[1]
      EdgeTypeL = Edge.FE[1]
      iFR = Edge.F[2]
      EdgeTypeR = Edge.FE[2]
      #computation normales of edges
      @views tBarLocL = tBar[:, EdgeTypeL] #* Grid.Faces[iFL].Orientation
      @views tBarLocR = tBar[:, EdgeTypeR] #* Grid.Faces[iFR].Orientation
      for iDoF = 1 : uFe.DoF
        ind = uFe.Glob[iDoF,iFL]  
        uLocL[iDoF] = u[ind] 
        ind = uFe.Glob[iDoF,iFR]  
        uLocR[iDoF] = u[ind] 
      end  
      @. qLocL = 0
      @. qLocR = 0
      for iQ = 1 : NumQuadL  
        @. uFLocL = 0
        @. uFLocR = 0
        for iDoF = 1 : uFe.DoF  
          uFLocL[1] += uFRef[1,iDoF,iQ,EdgeTypeL] * uLocL[iDoF]  
          uFLocL[2] += uFRef[2,iDoF,iQ,EdgeTypeL] * uLocL[iDoF]  
          uFLocR[1] += uFRef[1,iDoF,iQ,EdgeTypeR] * uLocR[iDoF]  
          uFLocR[2] += uFRef[2,iDoF,iQ,EdgeTypeR] * uLocR[iDoF]  
        end  
        Jacobi!(DFL,detDFL,pinvDFL,XL,ElemType,PointsE[1,iQ,EdgeTypeL],
          PointsE[2,iQ,EdgeTypeL],Grid.Faces[iFL], Grid)
        uuFLocL[1] = (DFL[1,1] * uFLocL[1] + DFL[1,2] * uFLocL[2]) / detDFL[1]
        uuFLocL[2] = (DFL[2,1] * uFLocL[1] + DFL[2,2] * uFLocL[2]) / detDFL[1]
        uuFLocL[3] = (DFL[3,1] * uFLocL[1] + DFL[3,2] * uFLocL[2]) / detDFL[1]
        ttL[1] = DFL[1,1] * tBarLocL[1]  + DFL[1,2] * tBarLocL[2] 
        ttL[2] = DFL[2,1] * tBarLocL[1]  + DFL[2,2] * tBarLocL[2] 
        ttL[3] = DFL[3,1] * tBarLocL[1]  + DFL[3,2] * tBarLocL[2] 
        Jacobi!(DFR,detDFR,pinvDFR,XR,ElemType,PointsE[1,iQ,EdgeTypeR],
          PointsE[2,iQ,EdgeTypeR],Grid.Faces[iFR], Grid)
        uuFLocR[1] = (DFR[1,1] * uFLocR[1] + DFR[1,2] * uFLocR[2]) / detDFR[1]
        uuFLocR[2] = (DFR[2,1] * uFLocR[1] + DFR[2,2] * uFLocR[2]) / detDFR[1]
        uuFLocR[3] = (DFR[3,1] * uFLocR[1] + DFR[3,2] * uFLocR[2]) / detDFR[1]
        ttR[1] = DFR[1,1] * tBarLocR[1]  + DFR[1,2] * tBarLocR[2] 
        ttR[2] = DFR[2,1] * tBarLocR[1]  + DFR[2,2] * tBarLocR[2] 
        ttR[3] = DFR[3,1] * tBarLocR[1]  + DFR[3,2] * tBarLocR[2] 
        #=
        tL = ((DFL[1,1] * uFLocL[1] + DFL[1,2] * uFLocL[2]) *
          (DFL[1,1] * tBarLocL[1] + DFL[1,2] * tBarLocL[2]) +
          (DFL[2,1] * uFLocL[1] + DFL[2,2] * uFLocL[2]) *
          (DFL[2,1] * tBarLocL[1] + DFL[2,2] * tBarLocL[2]) +
          (DFL[3,1] * uFLocL[1] + DFL[3,2] * uFLocL[2]) *
          (DFL[3,1] * tBarLocL[1] + DFL[3,2] * tBarLocL[2])) / detDFL[1] #* Grid.Faces[iF].Orientation
        tR = ((DFR[1,1] * uFLocR[1] + DFR[1,2] * uFLocR[2]) *
          (DFR[1,1] * tBarLocR[1] + DFR[1,2] * tBarLocR[2]) +
          (DFR[2,1] * uFLocR[1] + DFR[2,2] * uFLocR[2]) *
          (DFR[2,1] * tBarLocR[1] + DFR[2,2] * tBarLocR[2]) +
          (DFR[3,1] * uFLocR[1] + DFR[3,2] * uFLocR[2]) *
          (DFR[3,1] * tBarLocR[1] + DFR[3,2] * tBarLocR[2])) / detDFR[1] #* Grid.Faces[iF].Orientation
        =#  
        t1 = uuFLocL' * ttL 
        t2 = uuFLocR' * ttR 
        for iDoF = 1 : FeT.DoF  
#         qLocL[iDoF] += 0.5 * (-tL + tR) * qRef[iDoF,iQ,EdgeTypeL]  
#         qLocR[iDoF] += 0.5 * (-tL + tR) * qRef[iDoF,iQ,EdgeTypeR]  
          qLocL[iDoF] += 0.5 * WeightsL[iQ] * (t1 + t2) * qRef[iDoF,iQ,EdgeTypeL]  
          qLocR[iDoF] -= 0.5 * WeightsL[iQ] * (t1 + t2) * qRef[iDoF,iQ,EdgeTypeR]  
        end  
      end
      for iDoF = 1 : FeT.DoF
        ind = FeT.Glob[iDoF,iFL]
        q[ind] += qLocL[iDoF]
        if ind == 33836
          @show qLocL[iDoF]
        end  
        ind = FeT.Glob[iDoF,iFR]
        q[ind] += qLocR[iDoF]
        if ind == 33836
          @show qLocR[iDoF]
        end  
      end   
    end    
  end    
  ldiv!(FeT.LUM,q)
end
