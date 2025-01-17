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
      iF = Edge.F[1]
      EdgeType = Edge.FE[1]
      for iQ = 1 : NumQuadL
        Jacobi!(DF,detDF,pinvDF,X,ElemType,PointsE[1,iQ,EdgeType],
          PointsE[2,iQ,EdgeType],Grid.Faces[iF], Grid)
        nBar = Grid.nBar[:, EdgeType]
        n[1] = (pinvDF[1,1] * nBar[1] + pinvDF[1,2] * nBar[2]) * detDF[1]
        n[2] = (pinvDF[2,1] * nBar[1] + pinvDF[2,2] * nBar[2]) * detDF[1]
        n[3] = (pinvDF[3,1] * nBar[1] + pinvDF[3,2] * nBar[2]) * detDF[1]
        n *= Grid.Faces[iF].Orientation
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

function InterpolateRT!(u,FE,Jacobi,Grid,ElemType::Grids.Tri,QuadOrd,F)
  NumQuadT, WeightsT, PointsT = QuadRule(ElemType,QuadOrd)
  NumQuadL, WeightsL, PointsL = FEMSei.QuadRule(Grids.Line(),QuadOrd)
  k = FE.Order
  DoF = FE.DoF
  s = @polyvar x[1:2]

  if k > 0
    P_km1 = Polynomial_k(k-1,s)
    lP_km1 = length(P_km1)
  else
    lP_km1 = 0
  end  
  ValP_km1=zeros(NumQuadT,lP_km1)
  for iQ = 1 : NumQuadT
    for i = 1 : lP_km1
      ValP_km1[iQ,i] = P_km1[i](PointsT[iQ,1],PointsT[iQ,2])  
    end
  end  
  @polyvar t
  phiL = CGLine(k,t)
  l_phiL = length(phiL)
  ValphiL=zeros(NumQuadL,l_phiL)
  for iQ = 1 : NumQuadL
    for i = 1 : l_phiL
      ValphiL[iQ,i] = phiL[i](PointsL[iQ])  
    end
  end  
  uLoc = zeros(DoF)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  uP = zeros(2)
  VelSp = zeros(3)
  for iF = 1 : Grid.NumFaces
    iDoF = 1
    # Compute functional over edges
    # Edge 1 (-1,-1) -> (1,-1)
    @. uLoc = 0
    for iQ = 1 : NumQuadL
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,PointsL[iQ,1],-1,Grid.Faces[iF], Grid)
      h,VelSp[1],VelSp[2],VelSp[3], = F(X,0.0)
      lon,lat,r = Grids.cart2sphere(X[1],X[2],X[3])
      VelCa = VelSphere2Cart(VelSp,lon,lat)
      uP .= detDF[1] * pinvDF' * VelCa
      for i = 0 : k
        uLoc[iDoF+i] += 0.5 * Grid.Faces[iF].Orientation * uP[2] * ValphiL[iQ,i+1] * WeightsL[iQ]
      end
    end
    iDoF += k + 1
    # Edge 2 (1,-1) -> (-1,1)
    for iQ = 1 : NumQuadL
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,-PointsL[iQ,1],PointsL[iQ,1],Grid.Faces[iF], Grid)
      h,VelSp[1],VelSp[2],VelSp[3], = F(X,0.0)
      lon,lat,r = Grids.cart2sphere(X[1],X[2],X[3])
      VelCa = VelSphere2Cart(VelSp,lon,lat)
      uP .= detDF[1] * pinvDF' * VelCa
      for i = 0 : k
        uLoc[iDoF+i] += - 0.5 * Grid.Faces[iF].Orientation * (uP[1] + uP[2]) * ValphiL[iQ,i+1] * WeightsL[iQ] 
      end
    end
    iDoF += k + 1
    # Edge 3 (-1,1) -> (-1,-1)
    for iQ = 1 : NumQuadL
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,-1,-PointsL[iQ,1],Grid.Faces[iF], Grid)
      h,VelSp[1],VelSp[2],VelSp[3], = F(X,0.0)
      lon,lat,r = Grids.cart2sphere(X[1],X[2],X[3])
      VelCa = VelSphere2Cart(VelSp,lon,lat)
      uP .= detDF[1] * pinvDF' * VelCa
      for i = 0 : k
        uLoc[iDoF+i] += - 0.5 * Grid.Faces[iF].Orientation * uP[1] * ValphiL[iQ,i+1] * WeightsL[iQ]
      end
    end
    iDoF += k 
# Interior  
    for iQ = 1 : NumQuadT
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,PointsT[iQ,1],PointsT[iQ,2],Grid.Faces[iF], Grid)
       h,VelSp[1],VelSp[2],VelSp[3], = F(X,0.0)
      lon,lat,r = Grids.cart2sphere(X[1],X[2],X[3])
      VelCa = VelSphere2Cart(VelSp,lon,lat)
      uP .= detDF[1] * pinvDF' * VelCa
      for i = 1 : lP_km1
        phiLoc =  ValP_km1[iQ,i]
        uLoc[iDoF+2*i-1] += + 0.25 * Grid.Faces[iF].Orientation * uP[1] * phiLoc * WeightsT[iQ]
        uLoc[iDoF+2*i] += + 0.25 * Grid.Faces[iF].Orientation * uP[2] * phiLoc * WeightsT[iQ] 
      end
    end
    @. u[FE.Glob[:,iF]] = uLoc
  end  
end

function InterpolateRT!(u,FE,Jacobi,Grid,ElemType::Grids.Quad,QuadOrd,F)
  NumQuadT, WeightsT, PointsT = QuadRule(ElemType,QuadOrd)
  NumQuadL, WeightsL, PointsL = FEMSei.QuadRule(Grids.Line(),QuadOrd)
  k = FE.Order
  DoF = FE.DoF
  s = @polyvar x[1:2]
  P_kp1x1 = Polynomial_1D(k+1,s,1)
  P_kx1 = Polynomial_1D(k,s,1)
  P_kp1x2 = Polynomial_1D(k+1,s,2)
  P_kx2 = Polynomial_1D(k,s,2)
  if k > 0
    P_km1x1 = Polynomial_1D(k-1,s,1)
    P_km1x2 = Polynomial_1D(k-1,s,2)
  end
  
  DoF = 2 * (k+2) * (k+1)

  phi = Array{Polynomial,2}(undef,DoF,2)
  phiB = Array{Polynomial,2}(undef,DoF,2)
  rounded_poly = Array{Polynomial,2}(undef,DoF,2)
  Divphi = Array{Polynomial,2}(undef,DoF,1)
  rounded_Divphi = Array{Polynomial,2}(undef,DoF,1)
  iDoF = 1 
  for i = 1 : k+2
    for j = 1 : k+1
      phi[iDoF,1] = P_kp1x1[i] * P_kx2[j] 
      phi[iDoF,2] = 0.0 * x[1] + 0.0 * x[2]
      iDoF += 1
      phi[iDoF,1] = 0.0 * x[1] + 0.0 * x[2]
      phi[iDoF,2] = P_kp1x2[i] * P_kx1[j] 
      iDoF += 1
    end
  end
  @polyvar t
  phiL = CGLine(k,t)
  QuadOrd = 3
  NumQuadL, WeightsL, PointsL = FEMSei.QuadRule(Grids.Line(),QuadOrd)
  I = zeros(DoF,DoF)
  rDoF = 1
  ########
  phiL = CGLine(k,t)
  uLoc = zeros(DoF)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  uP = zeros(2)
  VelSp = zeros(3)
  for iF = 1 : Grid.NumFaces
    iDoF = 1
    rDoF = 1
    # Compute functional over edges
    # Edge 1 (-1,-1) -> (1,-1)
    @. uLoc = 0
    for iQ = 1 : NumQuadL
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,PointsL[iQ,1],-1,Grid.Faces[iF], Grid)
      h,VelSp[1],VelSp[2],VelSp[3], = F(X,0.0)
      lon,lat,r = Grids.cart2sphere(X[1],X[2],X[3])
      VelCa = VelSphere2Cart(VelSp,lon,lat) 
      uP .= detDF[1] * pinvDF' * VelCa
      for i = 0 : k
#       uLoc[iDoF+i] += 0.5 * Grid.Faces[iF].Orientation * uP[2] * phiL[i+1](PointsL[iQ]) * WeightsL[iQ]
        uLoc[iDoF+i] += 0.5 * uP[2] * phiL[i+1](PointsL[iQ]) * WeightsL[iQ]
      end
    end
    iDoF += k + 1
    # Edge 2 (1,-1) -> (1,1)
    for iQ = 1 : NumQuadL
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,1,PointsL[iQ,1],Grid.Faces[iF], Grid)
      h,VelSp[1],VelSp[2],VelSp[3], = F(X,0.0)
      lon,lat,r = Grids.cart2sphere(X[1],X[2],X[3])
      VelCa = VelSphere2Cart(VelSp,lon,lat) 
      uP .= detDF[1] * pinvDF' * VelCa
      for i = 0 : k
#       uLoc[iDoF+i] += 0.5 * Grid.Faces[iF].Orientation * uP[1] * phiL[i+1](PointsL[iQ]) * WeightsL[iQ] 
        uLoc[iDoF+i] += -0.5 * uP[1] * phiL[i+1](PointsL[iQ]) * WeightsL[iQ] 
      end
    end
    iDoF += k + 1
    # Edge 3 (-1,1) -> (1,1)
    for iQ = 1 : NumQuadL
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,PointsL[iQ,1],1,Grid.Faces[iF], Grid)
      h,VelSp[1],VelSp[2],VelSp[3], = F(X,0.0)
      lon,lat,r = Grids.cart2sphere(X[1],X[2],X[3])
      VelCa = VelSphere2Cart(VelSp,lon,lat)
      uP .= detDF[1] * pinvDF' * VelCa
      for i = 0 : k
#       uLoc[iDoF+i] += -0.5 * Grid.Faces[iF].Orientation * uP[2] * phiL[i+1](PointsL[iQ]) * WeightsL[iQ]
        uLoc[iDoF+i] += 0.5 * uP[2] * phiL[i+1](PointsL[iQ]) * WeightsL[iQ]
      end
    end
    iDoF += k + 1
    # Edge 4 (-1,-1) -> (-1,1)
    for iQ = 1 : NumQuadL
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,-1,PointsL[iQ,1],Grid.Faces[iF], Grid)
      h,VelSp[1],VelSp[2],VelSp[3], = F(X,0.0)
      lon,lat,r = Grids.cart2sphere(X[1],X[2],X[3])
      VelCa = VelSphere2Cart(VelSp,lon,lat)
      uP .= detDF[1] * pinvDF' * VelCa
      for i = 0 : k
#       uLoc[iDoF+i] += -0.5 * Grid.Faces[iF].Orientation * uP[1] * phiL[i+1](PointsL[iQ]) * WeightsL[iQ]
        uLoc[iDoF+i] += -0.5 * uP[1] * phiL[i+1](PointsL[iQ]) * WeightsL[iQ]
      end
    end
    iDoF += k + 1
# Interior  
    for iQ = 1 : NumQuadT
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,PointsT[iQ,1],PointsT[iQ,2],Grid.Faces[iF], Grid)
       h,VelSp[1],VelSp[2],VelSp[3], = F(X,0.0)
      lon,lat,r = Grids.cart2sphere(X[1],X[2],X[3])
      VelCa = VelSphere2Cart(VelSp,lon,lat)
      uP .= detDF[1] * pinvDF' * VelCa
      for i = 0 : 2 * (k - 1)
        phiLoc = phi[i+1](PointsT[iQ,1],PointsT[iQ,2])  
        uLoc[iDoF+i] += + 0.25 * Grid.Faces[iF].Orientation * uP[1] * phiLoc * WeightsT[iQ]
        uLoc[iDoF+i+1] += + 0.25 * Grid.Faces[iF].Orientation * uP[2] * phiLoc * WeightsT[iQ] 
      end
    end
    @. u[FE.Glob[:,iF]] = uLoc
  end  
end


