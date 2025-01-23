function ConvertScalar!(backend,FTB,pC,p,Fe::ScalarElement,Grid,Jacobi)
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

  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  for iF = 1 : Grid.NumFaces
    @views pC[iF]  = fRef[1, :]' * p[Fe.Glob[:,iF]]
  end
end

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
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  for iF = 1 : Grid.NumFaces
    Jacobi!(DF,detDF,pinvDF,X,Grid.Type,ksi1,ksi2,Grid.Faces[iF],Grid)
    detDFLoc = detDF[1]
    @views VelCa[iF,:]  .= (Grid.Faces[iF].Orientation / detDFLoc) * DF * (fRef[:, :] * Vel[Fe.Glob[:,iF]])
  VelCa[iF,3] = norm(VelCa[iF,:])
  end  
end

function ConvertVelocityCart!(backend,FTB,VelCa,Vel,Fe::VectorElement,Grid,Jacobi)
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
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  for iF = 1 : Grid.NumFaces
    Jacobi!(DF,detDF,pinvDF,X,Grid.Type,ksi1,ksi2,Grid.Faces[iF],Grid)
    detDFLoc = detDF[1]
    @views VelCa[iF,1] = (fRef[:, 1]' * Vel[Fe.Glob[:,iF]])
    @views VelCa[iF,2] = (fRef[:, 2]' * Vel[Fe.Glob[:,iF]])
    @views VelCa[iF,3] = (fRef[:, 3]' * Vel[Fe.Glob[:,iF]])
  end  
end

function ConvertVelocityCart!(backend,FTB,VelCa,Vel,Fe::HCurlElement,Grid,Jacobi)
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
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  for iF = 1 : Grid.NumFaces
    Jacobi!(DF,detDF,pinvDF,X,Grid.Type,ksi1,ksi2,Grid.Faces[iF],Grid)
    detDFLoc = detDF[1]
    @views VelCa[iF,:]  .= pinvDF * (fRef[:, :] * Vel[Fe.Glob[:,iF]])
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
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  for iF = 1 : Grid.NumFaces
    @. VelLoc = Vel[Fe.Glob[:,iF]] 
    Jacobi!(DF,detDF,pinvDF,X,Grid.Type,ksi1,ksi2,Grid.Faces[iF],Grid)
    detDFLoc = detDF[1] * Grid.Faces[iF].Orientation
    VelCa .= 1.0 / detDFLoc * DF * (fRef[:, :] * VelLoc) 
    lon,lat,_ = Grids.cart2sphere(X[1],X[2],X[3])
    VelSpLoc = VelCart2Sphere(VelCa,lon,lat)
    VelSp[iF,1] = VelSpLoc[1]
    VelSp[iF,2] = VelSpLoc[2]
  end  
end

function ConvertScalarVelocitySp!(backend,FTB,VelSp,Vel,Fe::HDivElement,h,hFe::ScalarElement,Grid,Jacobi)
  if Grid.Type == Grids.Tri()
    ksi1 = -1/3 
    ksi2 = -1/3
  elseif Grid.Type == Grids.Quad()
    ksi1 = 0
    ksi2 = 0
  end  
      
  fRef = zeros(Fe.Comp,Fe.DoF)
  hfRef = zeros(hFe.Comp,hFe.DoF)
  for iComp = 1 : Fe.Comp
    for iD = 1 : Fe.DoF
      fRef[iComp,iD] = Fe.phi[iD,iComp](ksi1,ksi2)
    end
  end
  for iComp = 1 : hFe.Comp
    for iD = 1 : hFe.DoF
      hfRef[iComp,iD] = hFe.phi[iD,iComp](ksi1,ksi2)
    end
  end

  VelLoc = zeros(Fe.DoF)
  VelCa = zeros(3)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  for iF = 1 : Grid.NumFaces
    VelLoc = Vel[Fe.Glob[:,iF]]  
    hLoc = h[hFe.Glob[:,iF]]
    Jacobi!(DF,detDF,pinvDF,X,Grid.Type,ksi1,ksi2,Grid.Faces[iF],Grid)
    detDFLoc = detDF[1]
    VelCa .= (Grid.Faces[iF].Orientation / detDFLoc) * DF * (fRef[:, :] * VelLoc) / (hfRef[:, :] * hLoc)
    lon,lat,_ = Grids.cart2sphere(X[1],X[2],X[3])
    VelSpLoc = VelCart2Sphere(VelCa,lon,lat)
    VelSp[iF,1] = VelSpLoc[1]
    VelSp[iF,2] = VelSpLoc[2]
  end  
end

function ConvertVelocitySp!(backend,FTB,VelSp,Vel,Fe::VectorElement,Grid,Jacobi)
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

  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  VelCa = zeros(3)
  X = zeros(3)
  for iF = 1 : Grid.NumFaces
    VelCa .= (fRef[:, :] * Vel[Fe.Glob[:,iF]])
    Jacobi!(DF,detDF,pinvDF,X,Grid.Type,ksi1,ksi2,Grid.Faces[iF],Grid)
    lon,lat,_ = Grids.cart2sphere(X[1],X[2],X[3])
    VelSpLoc = VelCart2Sphere(VelCa,lon,lat)
    VelSp[iF,1] = VelSpLoc[1]
    VelSp[iF,2] = VelSpLoc[2]
  end  
end

function ConvertVelocitySp!(backend,FTB,VelSp,Vel,Fe::HCurlElement,Grid,Jacobi)
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
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  for iF = 1 : Grid.NumFaces
    VelLoc = Vel[Fe.Glob[:,iF]]  
    Jacobi!(DF,detDF,pinvDF,X,Grid.Type,ksi1,ksi2,Grid.Faces[iF],Grid)
    detDFLoc = detDF[1]
    VelCa .= pinvDF * (fRef[:, :] * VelLoc) 
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
