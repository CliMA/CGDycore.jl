function InterpolateDG!(u,FE::DGStruct,Jacobi,Grid,ElemType,F)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  for iF = 1 : Grid.NumFaces
    for iDoF = 1 : FE.DoF
      Jacobi!(DF,detDF,pinvDF,X,ElemType,FE.points[iDoF,1],FE.points[iDoF,2],Grid.Faces[iF],Grid)  
      h, = F(X,0.0)
      ind = FE.Glob[iDoF,iF]
      u[ind] = h
    end    
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
        uLoc[iDoF+i] += -0.5 * Grid.Faces[iF].Orientation * (uP[1] + uP[2]) * ValphiL[iQ,i+1] * WeightsL[iQ] 
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
        uLoc[iDoF+i] += -0.5 * Grid.Faces[iF].Orientation * uP[1] * ValphiL[iQ,i+1] * WeightsL[iQ]
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
        uLoc[iDoF+2*i-1] += + 0.5 * Grid.Faces[iF].Orientation * uP[1] * phiLoc * WeightsT[iQ]
        uLoc[iDoF+2*i] += + 0.5 * Grid.Faces[iF].Orientation * uP[2] * phiLoc * WeightsT[iQ] 
      end
    end
    @. u[FE.Glob[:,iF]] = uLoc
  end  
end

function InterpolatehRT!(u,FE,Jacobi,Grid,ElemType::Grids.Tri,QuadOrd,F)
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
      uP .= detDF[1] * pinvDF' * VelCa * h
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
      uP .= detDF[1] * pinvDF' * VelCa * h
      for i = 0 : k
        uLoc[iDoF+i] += -0.5 * Grid.Faces[iF].Orientation * (uP[1] + uP[2]) * ValphiL[iQ,i+1] * WeightsL[iQ] 
      end
    end
    iDoF += k + 1
    # Edge 3 (-1,1) -> (-1,-1)
    for iQ = 1 : NumQuadL
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,-1,-PointsL[iQ,1],Grid.Faces[iF], Grid)
      h,VelSp[1],VelSp[2],VelSp[3], = F(X,0.0)
      lon,lat,r = Grids.cart2sphere(X[1],X[2],X[3])
      VelCa = VelSphere2Cart(VelSp,lon,lat)
      uP .= detDF[1] * pinvDF' * VelCa * h
      for i = 0 : k
        uLoc[iDoF+i] += -0.5 * Grid.Faces[iF].Orientation * uP[1] * ValphiL[iQ,i+1] * WeightsL[iQ]
      end
    end
    iDoF += k 
# Interior  
    for iQ = 1 : NumQuadT
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,PointsT[iQ,1],PointsT[iQ,2],Grid.Faces[iF], Grid)
      h,VelSp[1],VelSp[2],VelSp[3], = F(X,0.0)
      lon,lat,r = Grids.cart2sphere(X[1],X[2],X[3])
      VelCa = VelSphere2Cart(VelSp,lon,lat)
      uP .= detDF[1] * pinvDF' * VelCa * h
      for i = 1 : lP_km1
        phiLoc =  ValP_km1[iQ,i]
        uLoc[iDoF+2*i-1] += + 0.5 * Grid.Faces[iF].Orientation * uP[1] * phiLoc * WeightsT[iQ]
        uLoc[iDoF+2*i] += + 0.5 * Grid.Faces[iF].Orientation * uP[2] * phiLoc * WeightsT[iQ] 
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

function InterpolatehRT!(u,FE,Jacobi,Grid,ElemType::Grids.Quad,QuadOrd,F)
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
      uP .= detDF[1] * pinvDF' * VelCa * h
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
      uP .= detDF[1] * pinvDF' * VelCa * h
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
      uP .= detDF[1] * pinvDF' * VelCa * h
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
      uP .= detDF[1] * pinvDF' * VelCa * h
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
      uP .= detDF[1] * pinvDF' * VelCa * h
      for i = 0 : 2 * (k - 1)
        phiLoc = phi[i+1](PointsT[iQ,1],PointsT[iQ,2])  
        uLoc[iDoF+i] += + 0.25 * Grid.Faces[iF].Orientation * uP[1] * phiLoc * WeightsT[iQ]
        uLoc[iDoF+i+1] += + 0.25 * Grid.Faces[iF].Orientation * uP[2] * phiLoc * WeightsT[iQ] 
      end
    end
    @. u[FE.Glob[:,iF]] = uLoc
  end  
end

function InterpolateScalarHDivVecDG!(backend,FTB,uP,uFeP::VectorElement,h,hFe::ScalarElement,u,uFe::HDivElement,Grid,
  ElemType::Grids.ElementType,QuadOrd,Jacobi)
  @. uP = 0
  NumP = size(uFeP.points,1)
  ufRef  = zeros(uFe.Comp,uFe.DoF,NumP)
  hfRef  = zeros(hFe.Comp,hFe.DoF,NumP)
  for iP = 1 : NumP
    @inbounds for iComp = 1 : uFe.Comp
      @inbounds for iD = 1 : uFe.DoF
        ufRef[iComp,iD,iP] = uFe.phi[iD,iComp](uFeP.points[iP,1],uFeP.points[iP,2])
      end
    end
  end
  for iP = 1 : NumP
    @inbounds for iD = 1 : hFe.DoF
      hfRef[1,iD,iP] = hFe.phi[iD,1](uFeP.points[iP,1],uFeP.points[iP,2])
    end
  end
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  cPLoc = zeros(uFeP.DoF)
  ufRefLoc = zeros(2)
  @inbounds for iF = 1 : Grid.NumFaces
    iDoFVecDG = 1
    for iP = 1 : NumP
      @. ufRefLoc = 0
      @inbounds for iDoF = 1 : uFe.DoF
        ind = uFe.Glob[iDoF,iF]  
        @views ufRefLoc += ufRef[:,iDoF,iP] * u[ind]
      end  
      hfRefLoc = 0
      @inbounds for iDoF = 1 : hFe.DoF
        ind = hFe.Glob[iDoF,iF]  
        hfRefLoc += hfRef[1,iDoF,iP] * h[ind]
      end  
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,uFeP.points[iP,1],uFeP.points[iP,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1] * Grid.Faces[iF].Orientation
      uPLoc = 1 / detDFLoc * (DF * ufRefLoc) / hfRefLoc
      
      ind = uFeP.Glob[iDoFVecDG,iF]  
      uP[ind] = uPLoc[1]
      iDoFVecDG += 1

      ind = uFeP.Glob[iDoFVecDG,iF]  
      uP[ind] = uPLoc[2]
      iDoFVecDG += 1

      ind = uFeP.Glob[iDoFVecDG,iF] 
      uP[ind] = uPLoc[3]
      iDoFVecDG += 1
    end  
  end
end



