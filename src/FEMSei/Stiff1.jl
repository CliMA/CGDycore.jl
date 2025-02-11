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

function CurlVel(q,FeT,u,uFe::HDivElement,QuadOrd,ElemType,Grid,Jacobi)
#
#
# int q*v dx = int Curl u * v dx = - int u * rot v dx + int_E (u*t)*q ds
#
#
  @. q = 0

  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)  
  RotqRef  = zeros(2,FeT.DoF,NumQuad)
  uFRef  = zeros(uFe.Comp,uFe.DoF,NumQuad)
  for iQ = 1 : NumQuad
    for iDoF = 1 : uFe.DoF
      uFRef[1,iDoF,iQ] = uFe.phi[iDoF,1](Points[iQ,1],Points[iQ,2])
      uFRef[2,iDoF,iQ] = uFe.phi[iDoF,2](Points[iQ,1],Points[iQ,2])
    end  
    for iDoF = 1 : FeT.DoF
      RotqRef[1,iDoF,iQ] = -FeT.Gradphi[iDoF,1,2](Points[iQ,1],Points[iQ,2])
      RotqRef[2,iDoF,iQ] = FeT.Gradphi[iDoF,1,1](Points[iQ,1],Points[iQ,2])
    end  
  end  
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  uLoc = zeros(uFe.DoF)
  uFLoc = zeros(2)
  uuFLoc = zeros(3)
  qLoc = zeros(FeT.DoF)

  for iF = 1 : Grid.NumFaces
    for iDoF = 1 : uFe.DoF
      ind = uFe.Glob[iDoF,iF]  
      uLoc[iDoF] = u[ind] 
    end  
    @. qLoc = 0
    for iQ = 1 : NumQuad
      @. uFLoc = 0  
      for iDoF = 1 : uFe.DoF  
        uFLoc[1] += uFRef[1,iDoF,iQ] * uLoc[iDoF]  
        uFLoc[2] += uFRef[2,iDoF,iQ] * uLoc[iDoF]  
      end   
      Jacobi!(DF,detDF,pinvDF,X,ElemType,Points[iQ,1],Points[iQ,2],Grid.Faces[iF], Grid)
#     detDF[1] *= Grid.Faces[iF].Orientation
      uuFLoc[1] = (DF[1,1] * uFLoc[1] + DF[1,2] * uFLoc[2]) / detDF[1]
      uuFLoc[2] = (DF[2,1] * uFLoc[1] + DF[2,2] * uFLoc[2]) / detDF[1]
      uuFLoc[3] = (DF[3,1] * uFLoc[1] + DF[3,2] * uFLoc[2]) / detDF[1]
      uFLoc[1] = DF[1,1] * uuFLoc[1] + DF[2,1] * uuFLoc[2] + DF[3,1] * uuFLoc[3]
      uFLoc[2] = DF[1,2] * uuFLoc[1] + DF[2,2] * uuFLoc[2] + DF[3,2] * uuFLoc[3]
      for iDoF = 1 : FeT.DoF  
        qLoc[iDoF] +=  -Weights[iQ] * (uFLoc[1] * RotqRef[1,iDoF,iQ] +
          uFLoc[2] * RotqRef[2,iDoF,iQ])
      end  
    end  
    for iDoF = 1 : FeT.DoF
      ind = FeT.Glob[iDoF,iF]
      q[ind] += qLoc[iDoF]
    end
  end    
  @show sum(abs.(q))

  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  PointsE = zeros(2,NumQuadL,3)
  if ElemType == Grids.Tri()
    tBar = [1.0 -1.0  0.0
            0.0  1.0  1.0]
    NumE = 3
    PointsE = zeros(2,NumQuadL,NumE)
    for iQ = 1 : NumQuadL
      PointsE[1,iQ,1] = PointsL[iQ]
      PointsE[2,iQ,1] = -1.0
      PointsE[1,iQ,2] = -PointsL[iQ]
      PointsE[2,iQ,2] = PointsL[iQ]
      PointsE[1,iQ,3] = -1.0
      PointsE[2,iQ,3] = PointsL[iQ]
    end
  elseif ElemType == Grids.Quad()
    tBar = [1.0 0.0 1.0 0.0
            0.0 1.0 0.0 1.0]
    NumE = 4
    PointsE = zeros(2,NumQuadL,NumE)
    for iQ = 1 : NumQuadL
      PointsE[1,iQ,1] = PointsL[iQ]
      PointsE[2,iQ,1] = -1.0
      PointsE[1,iQ,2] = 1.0
      PointsE[2,iQ,2] = PointsL[iQ]
      PointsE[1,iQ,3] = PointsL[iQ]
      PointsE[2,iQ,3] = 1.0
      PointsE[1,iQ,4] = -1.0
      PointsE[2,iQ,4] = PointsL[iQ]
    end
  end  
  uFRef  = zeros(uFe.Comp,uFe.DoF,NumQuadL,NumE)
  qRef  = zeros(FeT.DoF,NumQuadL,NumE)
  for iQ = 1 : NumQuadL
    for iE = 1 : NumE  
      for iDoF = 1 : uFe.DoF
        uFRef[1,iDoF,iQ,iE] = uFe.phi[iDoF,1](PointsE[1,iQ,iE],PointsE[2,iQ,iE])
        uFRef[2,iDoF,iQ,iE] = uFe.phi[iDoF,2](PointsE[1,iQ,iE],PointsE[2,iQ,iE])
      end  
      for iDoF = 1 : FeT.DoF
        qRef[iDoF,iQ,iE] = FeT.phi[iDoF,1](PointsE[1,iQ,1],PointsE[2,iQ,1])
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
        detDFL[1] *= Grid.Faces[iFL].Orientation
        uuFLocL[1] = (DFL[1,1] * uFLocL[1] + DFL[1,2] * uFLocL[2]) / detDFL[1]
        uuFLocL[2] = (DFL[2,1] * uFLocL[1] + DFL[2,2] * uFLocL[2]) / detDFL[1]
        uuFLocL[3] = (DFL[3,1] * uFLocL[1] + DFL[3,2] * uFLocL[2]) / detDFL[1]
        ttL[1] = DFL[1,1] * tBarLocL[1]  + DFL[1,2] * tBarLocL[2] 
        ttL[2] = DFL[2,1] * tBarLocL[1]  + DFL[2,2] * tBarLocL[2] 
        ttL[3] = DFL[3,1] * tBarLocL[1]  + DFL[3,2] * tBarLocL[2] 
        Jacobi!(DFR,detDFR,pinvDFR,XR,ElemType,PointsE[1,iQ,EdgeTypeR],
          PointsE[2,iQ,EdgeTypeR],Grid.Faces[iFR], Grid)
        detDFR[1] *= Grid.Faces[iFR].Orientation
        uuFLocR[1] = (DFR[1,1] * uFLocR[1] + DFR[1,2] * uFLocR[2]) / detDFR[1]
        uuFLocR[2] = (DFR[2,1] * uFLocR[1] + DFR[2,2] * uFLocR[2]) / detDFR[1]
        uuFLocR[3] = (DFR[3,1] * uFLocR[1] + DFR[3,2] * uFLocR[2]) / detDFR[1]
        ttR[1] = DFR[1,1] * tBarLocR[1]  + DFR[1,2] * tBarLocR[2] 
        ttR[2] = DFR[2,1] * tBarLocR[1]  + DFR[2,2] * tBarLocR[2] 
        ttR[3] = DFR[3,1] * tBarLocR[1]  + DFR[3,2] * tBarLocR[2] 
        tL = 0.5 * (uuFLocL + uuFLocR)' * ttL * Grid.Faces[iFL].Orientation
        tR = 0.5 * (uuFLocL + uuFLocR)' * ttR * Grid.Faces[iFR].Orientation
        for iDoF = 1 : FeT.DoF  
          qLocL[iDoF] += Grid.Faces[iFL].OrientE[EdgeTypeL] * WeightsL[iQ] * 
            tL * qRef[iDoF,iQ,EdgeTypeL]  
          qLocR[iDoF] += Grid.Faces[iFR].OrientE[EdgeTypeR] * WeightsL[iQ] * 
            tR * qRef[iDoF,iQ,EdgeTypeR]  
        end  
      end
      for iDoF = 1 : FeT.DoF
        ind = FeT.Glob[iDoF,iFL]
        q[ind] += qLocL[iDoF]
        ind = FeT.Glob[iDoF,iFR]
        q[ind] += qLocR[iDoF]
      end   
    end    
  end    
  @show size(q)
  @show size(FeT.M)
  ldiv!(FeT.LUM,q)
end
function CurlVelSimple(q,FeT,u,uFe::HDivElement,QuadOrd,ElemType,Grid,F)
#
#
# int q*v dx = int Curl u * v dx = - int u * rot v dx + int_E (u*t)*q ds
#
#
  tBar = [1 0 1 0
          0 1 0 1]
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  @show NumQuadL
  @show WeightsL
  @show PointsL
  if ElemType == Grids.Quad()
    uFRef  = zeros(uFe.Comp,uFe.DoF,4)
    PointsE = zeros(2,4)
    PointsE[1,1] = 0.0
    PointsE[2,1] = -1.0
    PointsE[1,2] = 1.0
    PointsE[2,2] = 0.0
    PointsE[1,3] = 0.0
    PointsE[2,3] = 1.0
    PointsE[1,4] = -1.0
    PointsE[2,4] = 0.0
    for iDoF = 1 : uFe.DoF
      uFRef[1,iDoF,1] = uFe.phi[iDoF,1](0.0,-1.0)  
      uFRef[2,iDoF,1] = uFe.phi[iDoF,2](0.0,-1.0)  
      uFRef[1,iDoF,2] = uFe.phi[iDoF,1](1.0,0.0)
      uFRef[2,iDoF,2] = uFe.phi[iDoF,2](1.0,0.0)
      uFRef[1,iDoF,3] = uFe.phi[iDoF,1](0.0,1.0)  
      uFRef[2,iDoF,3] = uFe.phi[iDoF,2](0.0,1.0)  
      uFRef[1,iDoF,4] = uFe.phi[iDoF,1](-1.0,0.0)
      uFRef[2,iDoF,4] = uFe.phi[iDoF,2](-1.0,0.0)
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
  tLoc = zeros(3)
  VelSp = zeros(3)
  VelCa = zeros(3)

  @. q = 0
  q1 = similar(q)
  @. q1 = 0.0
  for iE = 1 : Grid.NumEdges
    Edge = Grid.Edges[iE]
    if length(Edge.F) > 1
      iFL = Edge.F[1]
      EdgeTypeL = Edge.FE[1]
      iFR = Edge.F[2]
      EdgeTypeR = Edge.FE[2]
      for iDoF = 1 : uFe.DoF
        ind = uFe.Glob[iDoF,iFL]  
        uLocL[iDoF] = u[ind] 
        ind = uFe.Glob[iDoF,iFR]  
        uLocR[iDoF] = u[ind] 
      end  
      @. uFLocL = 0
      @. uFLocR = 0
      for iDoF = 1 : uFe.DoF  
        uFLocL[1] += uFRef[1,iDoF,EdgeTypeL] * uLocL[iDoF]  
        uFLocL[2] += uFRef[2,iDoF,EdgeTypeL] * uLocL[iDoF]  
        uFLocR[1] += uFRef[1,iDoF,EdgeTypeR] * uLocR[iDoF]  
        uFLocR[2] += uFRef[2,iDoF,EdgeTypeR] * uLocR[iDoF]  
      end  
      Jacobi!(DFL,detDFL,pinvDFL,XL,ElemType,PointsE[1,EdgeTypeL],
        PointsE[2,EdgeTypeL],Grid.Faces[iFL], Grid)
      uuFLocL[1] = (DFL[1,1] * uFLocL[1] + DFL[1,2] * uFLocL[2]) 
      uuFLocL[2] = (DFL[2,1] * uFLocL[1] + DFL[2,2] * uFLocL[2]) 
      uuFLocL[3] = (DFL[3,1] * uFLocL[1] + DFL[3,2] * uFLocL[2])
      Jacobi!(DFR,detDFR,pinvDFR,XR,ElemType,PointsE[1,EdgeTypeR],
        PointsE[2,EdgeTypeR],Grid.Faces[iFR], Grid)
      uuFLocR[1] = (DFR[1,1] * uFLocR[1] + DFR[1,2] * uFLocR[2])
      uuFLocR[2] = (DFR[2,1] * uFLocR[1] + DFR[2,2] * uFLocR[2])
      uuFLocR[3] = (DFR[3,1] * uFLocR[1] + DFR[3,2] * uFLocR[2])
      tLoc[1] = Grid.Edges[iE].t.x * Grid.Edges[iE].a
      tLoc[2] = Grid.Edges[iE].t.y * Grid.Edges[iE].a
      tLoc[3] = Grid.Edges[iE].t.z * Grid.Edges[iE].a
      t = (uuFLocL + uuFLocR)' * tLoc / (detDFL[1] + detDFR[1])
      _,VelSp[1],VelSp[2],VelSp[3], = F(XL,0.0)
      lon,lat,r = Grids.cart2sphere(XL[1],XL[2],XL[3])
      VelCa = VelSphere2Cart(VelSp,lon,lat)
      t1 = VelCa' * tLoc
      q[iFL] += t * Grid.Faces[iFL].OrientE[EdgeTypeL] / Grid.Faces[iFL].Area
      q[iFR] += t * Grid.Faces[iFR].OrientE[EdgeTypeR] / Grid.Faces[iFR].Area
      q1[iFL] += t1 * Grid.Faces[iFL].OrientE[EdgeTypeL] / Grid.Faces[iFL].Area
      q1[iFR] += t1 * Grid.Faces[iFR].OrientE[EdgeTypeR] / Grid.Faces[iFR].Area
      if iFL == 38391
        @show iE  
        @show lon,lat
        @show VelCa
        @show uuFLocL
        @show uuFLocR
        @show t * Grid.Faces[iFL].OrientE[EdgeTypeL] / Grid.Faces[iFL].Area  
        @show t1 * Grid.Faces[iFL].OrientE[EdgeTypeL] / Grid.Faces[iFL].Area
        @show q[iFL],q1[iFL]
      end  
      if iFR == 38391
        @show iE  
        @show lon,lat
        @show VelCa
        @show uuFLocL
        @show uuFLocR
        @show t * Grid.Faces[iFR].OrientE[EdgeTypeR] / Grid.Faces[iFR].Area  
        @show t1 * Grid.Faces[iFR].OrientE[EdgeTypeR] / Grid.Faces[iFR].Area
        @show q[iFR],q1[iFR]
      end  
    end    
  end    
end

function CurlVelFun(q,FeT,u,uFe::HDivElement,QuadOrd,ElemType,Grid,F)
#
#
# int q*v dx = int Curl u * v dx = - int u * rot v dx + int_E (u*t)*q ds
#
#
  tBar = [1 0 1 0
          0 1 0 1]
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  @show NumQuadL
  @show WeightsL
  @show PointsL
  if ElemType == Grids.Quad()
    uFRef  = zeros(uFe.Comp,uFe.DoF,4)
    PointsE = zeros(2,4)
    PointsE[1,1] = 0.0
    PointsE[2,1] = -1.0
    PointsE[1,2] = 1.0
    PointsE[2,2] = 0.0
    PointsE[1,3] = 0.0
    PointsE[2,3] = 1.0
    PointsE[1,4] = -1.0
    PointsE[2,4] = 0.0
    for iDoF = 1 : uFe.DoF
      uFRef[1,iDoF,1] = uFe.phi[iDoF,1](0.0,-1.0)  
      uFRef[2,iDoF,1] = uFe.phi[iDoF,2](0.0,-1.0)  
      uFRef[1,iDoF,2] = uFe.phi[iDoF,1](1.0,0.0)
      uFRef[2,iDoF,2] = uFe.phi[iDoF,2](1.0,0.0)
      uFRef[1,iDoF,3] = uFe.phi[iDoF,1](0.0,1.0)  
      uFRef[2,iDoF,3] = uFe.phi[iDoF,2](0.0,1.0)  
      uFRef[1,iDoF,4] = uFe.phi[iDoF,1](-1.0,0.0)
      uFRef[2,iDoF,4] = uFe.phi[iDoF,2](-1.0,0.0)
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
  X = zeros(3)
  uLocL = zeros(uFe.DoF)
  uLocR = zeros(uFe.DoF)
  uFLocL = zeros(2)
  uFLocR = zeros(2)
  uuFLocL = zeros(3)
  uuFLocR = zeros(3)
  VelSp = zeros(3)
  VelCa = zeros(3)
  tLoc = zeros(3)

  @. q = 0
  for iE = 1 : Grid.NumEdges
    Edge = Grid.Edges[iE]
    if length(Edge.F) > 1
      iFL = Edge.F[1]
      EdgeTypeL = Edge.FE[1]
      iFR = Edge.F[2]
      EdgeTypeR = Edge.FE[2]
      for iDoF = 1 : uFe.DoF
        ind = uFe.Glob[iDoF,iFL]  
        uLocL[iDoF] = u[ind] 
        ind = uFe.Glob[iDoF,iFR]  
        uLocR[iDoF] = u[ind] 
      end  
      @. uFLocL = 0
      @. uFLocR = 0
      for iDoF = 1 : uFe.DoF  
        uFLocL[1] += uFRef[1,iDoF,EdgeTypeL] * uLocL[iDoF]  
        uFLocL[2] += uFRef[2,iDoF,EdgeTypeL] * uLocL[iDoF]  
        uFLocR[1] += uFRef[1,iDoF,EdgeTypeR] * uLocR[iDoF]  
        uFLocR[2] += uFRef[2,iDoF,EdgeTypeR] * uLocR[iDoF]  
      end  
      Jacobi!(DFL,detDFL,pinvDFL,XL,ElemType,PointsE[1,EdgeTypeL],
        PointsE[2,EdgeTypeL],Grid.Faces[iFL], Grid)
      X[1] = Grid.Edges[iE].Mid.x
      X[2] = Grid.Edges[iE].Mid.y
      X[3] = Grid.Edges[iE].Mid.z
      _,VelSp[1],VelSp[2],VelSp[3], = F(X,0.0)
      lon,lat,r = Grids.cart2sphere(X[1],X[2],X[3])
      VelCa = VelSphere2Cart(VelSp,lon,lat)
      tLoc[1] = Grid.Edges[iE].t.x * Grid.Edges[iE].a
      tLoc[2] = Grid.Edges[iE].t.y * Grid.Edges[iE].a
      tLoc[3] = Grid.Edges[iE].t.z * Grid.Edges[iE].a
      t = VelCa' * tLoc
      q[iFL] += t * Grid.Faces[iFL].OrientE[EdgeTypeL] / Grid.Faces[iFL].Area
      q[iFR] += t * Grid.Faces[iFR].OrientE[EdgeTypeR] / Grid.Faces[iFR].Area
    end    
  end    
end
