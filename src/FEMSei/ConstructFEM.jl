#=
function ConstructFEM(k)


# P0 1
# P1 1 + 2
# P2 1 + 2 + 3
# Pk 1 + 2 + ... + (k+1) = (k+1)*(k+2)/2
@polyvar x1 x2 
#RT0 
#P0^2 + (x,y)*P0
Dim = 2
k = 0
DoF = 2 * (k + 1) * (k + 2) / 2 + (k + 1)
phi = Array{Polynomial,2}(undef,DoF,Dim)
phi[1,1] = 1.0 
phi[1,2] = 0.0
k = 1
DoF = 2 * (k + 1) * (k + 2) / 2 + (k + 1)

end
=#

function ConstructRT_k(k)

  s = @polyvar x[1:2]

  P_km1 = Polynomial_k(k,s)
  H_km1 = HomegenuousPolynomial(k,s)

  lP_km1 = length(P_km1)
  lH_km1 = length(H_km1)
  DoF = 2 * lP_km1 + lH_km1

  phi = Array{Polynomial,2}(undef,DoF,2)
  phiB = Array{Polynomial,2}(undef,DoF,2)
  iDoF = 1 
  for i = 1 : lP_km1
    phi[iDoF,1] = P_km1[i]  
    phi[iDoF,2] = 0.0 * x[1] + 0.0 * x[2]
    iDoF += 1
    phi[iDoF,2] = P_km1[i]  
    phi[iDoF,1] = 0.0 * x[1] + 0.0 * x[2]
    iDoF += 1
  end  
  for i = 1 : lH_km1
    phi[iDoF,1] = H_km1[i] * x[1]
    phi[iDoF,2] = H_km1[i] * x[2]
    iDoF += 1
  end  
  @polyvar t
  phiL = CGLine(k,t)
  QuadOrd = 3
  NumQuadL, WeightsL, PointsL = FEMSei.QuadRule(Grids.Line(),QuadOrd)
  I = zeros(DoF,DoF)
  rDoF = 1
# Compute functional over edges
  # Edge 1 (-1,-1) -> (1,-1)
  for iDoF = 1 : DoF
    phiE2 = subs(phi[iDoF,2], x[1] => t, x[2] => -1.0)
    for i = 0 : k
      for iQ = 1 : NumQuadL
        I[rDoF+i,iDoF] += 0.5 * phiE2(PointsL[iQ]) * phiL[i+1](PointsL[iQ]) * WeightsL[iQ]
      end
    end
  end
  rDoF += k + 1
  # Edge 2 (1,-1) -> (-1,1)
  for iDoF = 1 : DoF
    phiE1 = subs(phi[iDoF,1], x[1] => -t, x[2] => t)
    phiE2 = subs(phi[iDoF,2], x[1] => -t, x[2] => t)
    for i = 0 : k
      for iQ = 1 : NumQuadL
        I[rDoF+i,iDoF] += -0.5 * (phiE1(PointsL[iQ]) + phiE2(PointsL[iQ])) * phiL[i+1](PointsL[iQ]) * WeightsL[iQ] 
      end
    end
  end
  rDoF += k + 1
# Edge 3 (-1,1) -> (-1,-1)
  for iDoF = 1 : DoF
    phiE1 = subs(phi[iDoF,1], x[1] => -1, x[2] => -t)
    for i = 0 : k
      for iQ = 1 : NumQuadL
        I[rDoF+i,iDoF] += -0.5 * phiE1(PointsL[iQ]) * phiL[i+1](PointsL[iQ]) * WeightsL[iQ]  
      end  
    end  
  end  
  rDoF += k + 1
  NumQuadT, WeightsT, PointsT = FEMSei.QuadRule(Grids.Tri(),QuadOrd)
# Interior  
  for iDoF = 1 : DoF
    phiI1 = phi[iDoF,1]  
    phiI2 = phi[iDoF,2]  
    for i = 0 : 2 * (k - 1)
      for iQ = 1 : NumQuadT
        I[rDoF+i,iDoF] += 0.25 * phiI1(PointsT[iQ,1],PointsT[iQ,2]) * WeightsT[iQ]
        I[rDoF+i+1,iDoF] += 0.25 * phiI2(PointsT[iQ,1],PointsT[iQ,2]) * WeightsT[iQ]
      end
      phiI1 = phiI1 * x[1]
      phiI2 = phiI2 * x[2]
    end
  end
  for iDoF = 1 : DoF  
    for jDoF = 1 : DoF  
      if abs(I[iDoF,jDoF]) < 1.e-12
        I[iDoF,jDoF] = 0
      end
    end
  end  
  r = zeros(DoF)
  for iDoF = 1 : DoF  
    r[iDoF] = 1
    c = I \ r
    phiB[iDoF,1] = 0.0 * x[1] + 0.0 * x[2]
    phiB[iDoF,2] = 0.0 * x[1] + 0.0 * x[2]
    for jDoF = 1 : DoF  
      phiB[iDoF,:] += c[jDoF] * phi[jDoF,:]  
    end  
    r[iDoF] = 0
  end  
  return phi, phiB
end

function CGLine(k,x)
  phi = Array{Polynomial,1}(undef,k+1)
  for i = 0 : k
    phi[i+1] = 0.5^k * (1.0-x)^(k-i)*(1.0+x)^i
  end  
  return phi
end

function HomegenuousPolynomial(k,x)
  phi = Array{Polynomial,1}(undef,k+1)
  for i = 0 : k
    phi[i+1] = x[1][1]^(k-i) * x[1][2]^i + 0.0 
  end  
  return phi
end  

function Polynomial_k(k,x)
  DoF::Int = (k + 1) *(k + 2) / 2
  phi = Array{Polynomial,1}(undef,DoF)
  iDoF = 1
  for i = 0 : k
    for j = 0 : i  
      phi[iDoF] = x[1][1]^(i-j) * x[1][2]^j + 0.0 
      iDoF += 1
    end
  end
  return phi
end

function InterpolateRT!(u,FE,Jacobi,Grid,F)

  k = FE.Order
  DoF = FE.DoF
  s = @polyvar x[1:2]

  P_km1 = Polynomial_k(k,s)
  H_km1 = HomegenuousPolynomial(k,s)

  lP_km1 = length(P_km1)
  lH_km1 = length(H_km1)
  DoF = 2 * lP_km1 + lH_km1

  phi = Array{Polynomial,2}(undef,DoF,2)
  phiB = Array{Polynomial,2}(undef,DoF,2)
  iDoF = 1
  for i = 1 : lP_km1
    phi[iDoF,1] = P_km1[i]
    phi[iDoF,2] = 0.0 * x[1] + 0.0 * x[2]
    iDoF += 1
    phi[iDoF,2] = P_km1[i]
    phi[iDoF,1] = 0.0 * x[1] + 0.0 * x[2]
    iDoF += 1
  end
  for i = 1 : lH_km1
    phi[iDoF,1] = H_km1[i] * x[1]
    phi[iDoF,2] = H_km1[i] * x[2]
    iDoF += 1
  end
  @polyvar t
  phiL = CGLine(k,t)
  QuadOrd = 3
  NumQuadL, WeightsL, PointsL = FEMSei.QuadRule(Grids.Line(),QuadOrd)
  NumQuadT, WeightsT, PointsT = FEMSei.QuadRule(Grids.Tri(),QuadOrd)
  uLoc = zeros(DoF)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  uP = zeros(2)
  VelSp = zeros(3)
  @show phi
  for iF = 1 : Grid.NumFaces
    iDoF = 1
    rDoF = 1
    # Compute functional over edges
    # Edge 1  
    @. uLoc = 0
    for iQ = 1 : NumQuadL
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,PointsL[iQ,1],-1,Grid.Faces[iF], Grid)
      _,VelSp[1],VelSp[2],VelSp[3], = F(X,0.0)
      lon,lat,r = Grids.cart2sphere(X[1],X[2],X[3])
      VelCa = VelSphere2Cart(VelSp,lon,lat)
      uP .= detDF[1] * pinvDF' * VelCa
      for i = 0 : k
        uLoc[iDoF+i] += Grid.Faces[iF].Orientation * uP[2] * phiL[i+1](PointsL[iQ]) * WeightsL[iQ]
      end
    end
    iDoF += k + 1
    # Edge 2 (1,-1) -> (-1,1)
    for iQ = 1 : NumQuadL
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,-PointsL[iQ,1],PointsL[iQ,1],Grid.Faces[iF], Grid)
      _,VelSp[1],VelSp[2],VelSp[3], = F(X,0.0)
      lon,lat,r = Grids.cart2sphere(X[1],X[2],X[3])
      VelCa = VelSphere2Cart(VelSp,lon,lat)
      uP .= detDF[1] * pinvDF' * VelCa
      for i = 0 : k
        uLoc[iDoF+i] += -Grid.Faces[iF].Orientation * (uP[1] + uP[2]) * phiL[i+1](PointsL[iQ]) * WeightsL[iQ] 
      end
    end
    iDoF += k + 1
    # Edge 3 (-1,1) -> (-1,-1)
    for iQ = 1 : NumQuadL
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,-1,-PointsL[iQ,1],Grid.Faces[iF], Grid)
      _,VelSp[1],VelSp[2],VelSp[3], = F(X,0.0)
      lon,lat,r = Grids.cart2sphere(X[1],X[2],X[3])
      VelCa = VelSphere2Cart(VelSp,lon,lat)
      uP .= detDF[1] * pinvDF' * VelCa
      for i = 0 : k
        uLoc[iDoF+i] += -Grid.Faces[iF].Orientation * uP[1] * phiL[i+1](PointsL[iQ]) * WeightsL[iQ]
      end
    end
    iDoF += k + 1
#   if Grid.Faces[iF].E[1] == 1
#     @show uLoc[1],Grid.Faces[iF].OrientE[1],Grid.Faces[iF].Orientation   
#   elseif Grid.Faces[iF].E[2] == 1
#     @show uLoc[2],Grid.Faces[iF].OrientE[2],Grid.Faces[iF].Orientation   
#   elseif Grid.Faces[iF].E[3] == 1
#     @show uLoc[3],Grid.Faces[iF].OrientE[3],Grid.Faces[iF].Orientation   
#   end  
# Interior  
    rDoF += k + 1
    for iQ = 1 : NumQuadT
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,PointsT[iQ,1],PointsT[iQ,2],Grid.Faces[iF], Grid)
       _,VelSp[1],VelSp[2],VelSp[3], = F(X,0.0)
      lon,lat,r = Grids.cart2sphere(X[1],X[2],X[3])
      VelCa = VelSphere2Cart(VelSp,lon,lat)
      uP .= detDF[1] * pinvDF' * VelCa
      for i = 0 : 2 * (k - 1)
        phiLoc = phi[i+1](PointsT[iQ,1],PointsT[iQ,2])  
        uLoc[rDoF+i] += -Grid.Faces[iF].Orientation * uP[1] * phiLoc * WeightsT[iQ]
        uLoc[rDoF+i+1] += -Grid.Faces[iF].Orientation * uP[2] * phiLoc * WeightsT[iQ]
      end
    end
    @. u[FE.Glob[:,iF]] = uLoc
  end  
end
