function Interpolate1!(backend,FTB,uN,Fe::HDivConfElement,ElemType,Grid,QuadOrd,Jacobi,F)
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  X = zeros(3)
  VelSp = zeros(3)
  VelCa = zeros(3)
  nR = zeros(3)
  nL = zeros(3)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)

  if ElemType == Grids.Quad()
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
    end
  elseif ElemType == Grids.Tri()
    PointsE = zeros(2,NumQuadL,3)
    for iQ = 1 : NumQuadL
      PointsE[1,iQ,1] = PointsL[iQ]
      PointsE[2,iQ,1] = -1.0
      PointsE[1,iQ,2] = -PointsL[iQ]
      PointsE[2,iQ,2] = PointsL[iQ]
      PointsE[1,iQ,3] = -1.0
      PointsE[2,iQ,3] = PointsL[iQ]
    end
  end  
  @. uN = 0
  for iE = 1 : Grid.NumEdges
    Edge = Grid.Edges[iE]
    if length(Edge.F) > 1
      iFL = Edge.F[1]
      EdgeTypeL = Edge.FE[1]
      for iQ = 1 : NumQuadL
        Jacobi!(DF,detDF,pinvDF,X,ElemType,PointsE[1,iQ,EdgeTypeL],
          PointsE[2,iQ,EdgeTypeL],Grid.Faces[iFL], Grid)
        nBarL = Grid.nBar[:, EdgeTypeL]
        detDFL = detDF[1] * Grid.Faces[iFL].Orientation
        nL[1] = (pinvDF[1,1] * nBarL[1] + pinvDF[1,2] * nBarL[2]) * detDFL
        nL[2] = (pinvDF[2,1] * nBarL[1] + pinvDF[2,2] * nBarL[2]) * detDFL
        nL[3] = (pinvDF[3,1] * nBarL[1] + pinvDF[3,2] * nBarL[2]) * detDFL
        _,VelSp[1],VelSp[2],VelSp[3], = F(X,0.0)
        lon,lat,r = Grids.cart2sphere(X[1],X[2],X[3])
        VelCa = VelSphere2Cart(VelSp,lon,lat)
        ind = Fe.Glob[EdgeTypeL,iFL]
        uN[iE] -= 0.5 * nL' * VelCa * WeightsL[iQ] 
      end 
    end  
  end  
end

function Interpolate!(backend,FTB,uN,Fe::HDivConfElement,ElemType,Grid,QuadOrd,Jacobi,F)
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  X = zeros(3)
  VelSp = zeros(3)
  VelCa = zeros(3)
  n = zeros(3)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)

  if ElemType == Grids.Quad()
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
    end
  end  
  @. uN = 0
  for iE = 1 : Grid.NumEdges
    Edge = Grid.Edges[iE]
    if length(Edge.F) > 1
      iF = Edge.F[1]
      EdgeType = Edge.FE[1]
      for iQ = 1 : NumQuadL
        Jacobi!(DF,detDF,pinvDF,X,ElemType,PointsE[1,iQ,EdgeType],
          PointsE[2,iQ,EdgeType],Grid.Faces[iF], Grid)
        nBar = Grid.nBar[:, EdgeType]
        n[1] = (pinvDF[1,1] * nBar[1] + pinvDF[1,2] * nBar[2]) * detDF[1]
        n[2] = (pinvDF[2,1] * nBar[1] + pinvDF[2,2] * nBar[2]) * detDF[1]
        n[3] = (pinvDF[3,1] * nBar[1] + pinvDF[3,2] * nBar[2]) * detDF[1]
        _,VelSp[1],VelSp[2],VelSp[3], = F(X,0.0)
        lon,lat,r = Grids.cart2sphere(X[1],X[2],X[3])
        VelCa = VelSphere2Cart(VelSp,lon,lat)
        ind = Fe.Glob[EdgeType,iF]
        uN[iE] += - 0.5 * n' * VelCa * WeightsL[iQ]
      end  
    end  
  end  
end

function InterpolateCons!(backend,FTB,uN,Fe::HDivConfElement,Grid,QuadOrd,Jacobi,F)
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  X = zeros(3)
  VelSp = zeros(3)
  for iE = 1 : Grid.NumEdges
    X[1] = Grid.Edges[iE].Mid.x  
    X[2] = Grid.Edges[iE].Mid.y  
    X[3] = Grid.Edges[iE].Mid.z  
    h,VelSp[1],VelSp[2],VelSp[3], = F(X,0.0)
    lon,lat,r = Grids.cart2sphere(X[1],X[2],X[3])
    VelCa = VelSphere2Cart(VelSp,lon,lat) * h
    uN[iE] = -0.5 * Grid.Edges[iE].a * (Grid.Edges[iE].n.x * VelCa[1] +   
      Grid.Edges[iE].n.y * VelCa[2] + Grid.Edges[iE].n.z * VelCa[3])  
  end  
end

function Project(backend,FTB,Fe::ScalarElement,Grid,QuadOrd,Jacobi,F)
  NumQuad,Weights,Points = QuadRule(Fe.Type,QuadOrd)
  fRef  = zeros(Fe.Comp,Fe.DoF,length(Weights))

  p=zeros(Fe.NumG)
  for i = 1 : length(Weights)
    for iComp = 1 : Fe.Comp
      for iD = 1 : Fe.DoF
        fRef[iComp,iD,i] = Fe.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  for iF = 1 : Grid.NumFaces
    pLoc = zeros(Fe.Comp,Fe.DoF)
    for i = 1 : length(Weights)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      f, = F(X,0.0)
      pLoc += abs(detDFLoc)*Weights[i]*(fRef[:,:,i]*f)
    end
    @. p[Fe.Glob[:,iF]] += pLoc[Fe.Comp,:]
  end
  ldiv!(Fe.LUM,p)
  return p
end

