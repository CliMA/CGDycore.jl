function ConvertVelocityCart!(backend,FTB,VelCa,Vel,Fe::HDivElement,Grid,Jacobi)
  if Grid.Type == Grids.Tri()
    ksi1 = -1/3 
    ksi2 = -1/3
  elseif Grid.Type == Grids.Quad()
    ksi1 = 0
    ksi2 = 0
  end  
      
  fRef = zeros(Fe.Comp,Fe.DoF)
  for iComp = 1 : Fe.Comp
    for iD = 1 : Fe.DoF
      fRef[iComp,iD] = Fe.phi[iD,iComp](ksi1,ksi2)
    end
  end

  VelLoc = zeros(Fe.DoF)
  for iF = 1 : Grid.NumFaces
    DF, detJ = Jacobi(Grid.Type,ksi1,ksi2,Grid.Faces[iF],Grid)
    @views VelCa[iF,:]  .= (1 / detJ) * DF * (fRef[:, :] * Vel[Fe.Glob[:,iF]])
  end  
end

function ConvertVelocitySp!(backend,FTB,VelSp,Vel,Fe::HDivElement,Grid,Jacobi)
  if Grid.Type == Grids.Tri()
    ksi1 = -1/3 
    ksi2 = -1/3
  elseif Grid.Type == Grids.Quad()
    ksi1 = 0
    ksi2 = 0
  end  
      
  fRef = zeros(Fe.Comp,Fe.DoF)
  for iComp = 1 : Fe.Comp
    for iD = 1 : Fe.DoF
      fRef[iComp,iD] = Fe.phi[iD,iComp](ksi1,ksi2)
    end
  end

  VelLoc = zeros(Fe.DoF)
  VelCa = zeros(3)
  for iF = 1 : Grid.NumFaces
    VelLoc = Vel[Fe.Glob[:,iF]]  
    DF, detJ, X = Jacobi(Grid.Type,ksi1,ksi2,Grid.Faces[iF],Grid)
    VelCa .= (1 / detJ) * DF * (fRef[:, :] * VelLoc) 
    lon,lat,_ = Grids.cart2sphere(X[1],X[2],X[3])
    VelSpLoc = VelCart2Sphere(VelCa,lon,lat)
    VelSp[iF,1] = VelSpLoc[1]
    VelSp[iF,2] = VelSpLoc[2]
  end  
end

function VelCart2Sphere(VelCa,lon,lat)

  Rotate = zeros(3,3)
  Rotate[1,1] =           -sin(lon)
  Rotate[2,1] =-sin(lat) * cos(lon)
  Rotate[3,1] = cos(lat) * cos(lon)
  Rotate[1,2] =           cos(lon)
  Rotate[2,2] =-sin(lat) * sin(lon)
  Rotate[3,2] = cos(lat) * sin(lon)
  Rotate[1,3] =          0.0
  Rotate[2,3] = cos(lat)
  Rotate[3,3] = sin(lat)
  VelSp = zeros(3,1)
  VelSp[1] = Rotate[1,1] * VelCa[1] +
    Rotate[1,2] * VelCa[2] +
    Rotate[1,3] * VelCa[3]
  VelSp[2] = Rotate[2,1] * VelCa[1] +
    Rotate[2,2] * VelCa[2] +
    Rotate[2,3] * VelCa[3]
  VelSp[3] = Rotate[3,1]*VelCa[1] +
    Rotate[3,2] * VelCa[2] +
    Rotate[3,3] * VelCa[3]
  return VelSp  
end

function VelSphere2Cart(VelSp,lon,lat)
  rot=zeros(3,3);
  rot[1,1]=-sin(lon);
  rot[2,1]=-sin(lat)*cos(lon);
  rot[3,1]= cos(lat)*cos(lon);
  rot[1,2]=cos(lon);
  rot[2,2]=-sin(lat)*sin(lon);
  rot[3,2]=cos(lat)*sin(lon);
  rot[1,3]=0.0;
  rot[2,3]=cos(lat);
  rot[3,3]=sin(lat);
  VelCa = rot'*VelSp
end
