function Vorticity!(backend,FTB,Vort,VortFE::ScalarElement,hu,uFE::HDivElement,h,
  hFE::ScalarElement,ND::HCurlElement,Curl,Grid,ElemType,QuadOrd,Jacobi)
  UCacheu = zeros(size(hu))
  UCachep = zeros(VortFE.NumG)
  ProjectHDivScalarHCurl!(backend,FTB,UCacheu,ND,hu,uFE,h,hFE,
    Grid,ElemType,QuadOrd,Jacobi)
  mul!(UCachep,Curl,UCacheu)
  ldiv!(VortFE.LUM,UCachep)
  ConvertScalar!(backend,FTB,Vort,UCachep,VortFE,Grid,Jacobi)
end

function Vorticity!(backend,FTB,Vort,VortFE::ScalarElement,hu,uFE::HDivElement,h,
  hFE::ScalarElement,ND::HCurlElement,Curl,Grid,ElemType,QuadOrd,Jacobi,ksi)
  UCacheu = zeros(size(hu))
  UCachep = zeros(VortFE.NumG)
  ProjectHDivScalarHCurl!(backend,FTB,UCacheu,ND,hu,uFE,h,hFE,
    Grid,ElemType,QuadOrd,Jacobi)
  mul!(UCachep,Curl,UCacheu)
  ldiv!(VortFE.LUM,UCachep)
  ConvertScalar!(backend,FTB,Vort,UCachep,VortFE,Grid,Jacobi,ksi)
end

function Vorticity!(backend,FTB,Vort,VortFE::ScalarElement,u,uFE::HDivElement,
  Grid,ElemType,QuadOrd,Jacobi,ksi)
  UCachep = zeros(VortFE.NumG)
  CurlVel!(UCachep,VortFE,u,uFE,QuadOrd,ElemType,Grid,Jacobi)
  ConvertScalar!(backend,FTB,Vort,UCachep,VortFE,Grid,Jacobi,ksi)
end

function Vorticity!(backend,FTB,Vort,VortFE::ScalarElement,u,uFE::HDivElement,
  h,hFE::ScalarElement,Grid,ElemType,QuadOrd,Jacobi,ksi)
  UCachep = zeros(VortFE.NumG)
  CurlVel!(UCachep,VortFE,u,uFE,h,hFE,QuadOrd,ElemType,Grid,Jacobi)
  ConvertScalar!(backend,FTB,Vort,UCachep,VortFE,Grid,Jacobi,ksi)
end

function ConvertScalar!(backend,FTB,pC,p,Fe::ScalarElement,Grid,Jacobi)
  if Grid.Type == Grids.Tri()
    ksi1 = -1/3
    ksi2 = -1/3
  elseif Grid.Type == Grids.Quad()
    ksi1 = 0
    ksi2 = 0
  end

  fRef = zeros(Fe.Comp,Fe.DoF)
  @inbounds for iComp = 1 : Fe.Comp
    @inbounds for iD = 1 : Fe.DoF
      fRef[iComp,iD] = Fe.phi[iD,iComp](ksi1,ksi2)
    end
  end

  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  @inbounds for iF = 1 : Grid.NumFaces
    @views pC[iF]  = fRef[1, :]' * p[Fe.Glob[:,iF]]
  end
end

function ConvertScalar!(backend,FTB,pC,p,Fe::ScalarElement,Grid,Jacobi,ksi)
  Numksi = size(ksi,1)
  fRef = zeros(Numksi,Fe.Comp,Fe.DoF)
  @inbounds for iksi = 1 : Numksi
    @inbounds for iComp = 1 : Fe.Comp
      @inbounds for iD = 1 : Fe.DoF
        fRef[iksi,iComp,iD] = Fe.phi[iD,iComp](ksi[iksi,1],ksi[iksi,2])
      end
    end
  end
  iS = 0
  @inbounds for iF = 1 : Grid.NumFaces
    @inbounds for iksi = 1 : Numksi
      iS += 1  
      @views pC[iS]  = fRef[iksi,1, :]' * p[Fe.Glob[:,iF]]
    end  
  end
end

function ConvertVelocity!(backend,FTB,VelCa,Vel,Fe::HDivElement,Grid,Jacobi,::Grids.CartesianGrid)
  if Grid.Type == Grids.Tri()
    ksi1 = -1/3 
    ksi2 = -1/3
  elseif Grid.Type == Grids.Quad()
    ksi1 = 0
    ksi2 = 0
  end  
      
  fRef = zeros(Fe.Comp,Fe.DoF)
  @inbounds for iComp = 1 : Fe.Comp
    @inbounds for iD = 1 : Fe.DoF
      fRef[iComp,iD] = Fe.phi[iD,iComp](ksi1,ksi2)
    end
  end

  VelLoc = zeros(Fe.DoF)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  @inbounds for iF = 1 : Grid.NumFaces
    Jacobi(DF,detDF,pinvDF,X,Grid.Type,ksi1,ksi2,Grid.Faces[iF],Grid)
    detDFLoc = detDF[1]
    @views VelCa[iF,:]  .= (Grid.Faces[iF].Orientation / detDFLoc) * DF * (fRef[:, :] * Vel[Fe.Glob[:,iF]])
  end  
end

function ConvertVelocity!(backend,FTB,VelSp,Vel,Fe::HDivElement,Grid,Jacobi,::Grids.SphericalGrid)
  if Grid.Type == Grids.Tri()
    ksi1 = -1/3
    ksi2 = -1/3
  elseif Grid.Type == Grids.Quad()
    ksi1 = 0
    ksi2 = 0
  end

  fRef = zeros(Fe.Comp,Fe.DoF)
  @inbounds for iComp = 1 : Fe.Comp
    @inbounds for iD = 1 : Fe.DoF
      fRef[iComp,iD] = Fe.phi[iD,iComp](ksi1,ksi2)
    end
  end

  VelLoc = zeros(Fe.DoF)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  VelCa = zeros(3)
  @inbounds for iF = 1 : Grid.NumFaces
    Jacobi(DF,detDF,pinvDF,X,Grid.Type,ksi1,ksi2,Grid.Faces[iF],Grid)
    detDFLoc = detDF[1]
    @views VelCa  .= (Grid.Faces[iF].Orientation / detDFLoc) * DF * (fRef[:, :] * Vel[Fe.Glob[:,iF]])
    lon,lat,_ = Grids.cart2sphere(X[1],X[2],X[3])
    VelSpLoc = Grids.VelCart2Sphere(VelCa,lon,lat)
    VelSp[iF,1] = VelSpLoc[1]
    VelSp[iF,2] = VelSpLoc[2]
  end
end

function ConvertVelocity!(backend,FTB,VelCa,Vel,Fe::VectorElement,Grid,Jacobi,::Grids.CartesianGrid)
  if Grid.Type == Grids.Tri()
    ksi1 = -1/3 
    ksi2 = -1/3
  elseif Grid.Type == Grids.Quad()
    ksi1 = 0
    ksi2 = 0
  end  
      
  fRef = zeros(Fe.Comp,Fe.DoF)
  @inbounds for iComp = 1 : Fe.Comp
    @inbounds for iD = 1 : Fe.DoF
      fRef[iComp,iD] = Fe.phi[iD,iComp](ksi1,ksi2)
    end
  end

  @inbounds for iF = 1 : Grid.NumFaces
    @views VelCa[iF,1] = (fRef[:, 1]' * Vel[Fe.Glob[:,iF]])
    @views VelCa[iF,2] = (fRef[:, 2]' * Vel[Fe.Glob[:,iF]])
    @views VelCa[iF,3] = (fRef[:, 3]' * Vel[Fe.Glob[:,iF]])
  end  
end

function ConvertVelocityCart!(backend,FTB,VelCa,Vel,Fe::VectorElement,Grid,Jacobi,ksi)

  Numksi = size(ksi,1)    
  fRef = zeros(Numksi,Fe.Comp,Fe.DoF)
  @inbounds for iksi = 1 : Numksi
    @inbounds for iComp = 1 : Fe.Comp
      @inbounds for iD = 1 : Fe.DoF
        fRef[iksi,iComp,iD] = Fe.phi[iD,iComp](ksi[iksi,1],ksi[iksi,2])
      end
    end
  end

  iS = 0
  @inbounds for iF = 1 : Grid.NumFaces
    @inbounds for iksi = 1 : Numksi
      iS += 1
      @views VelCa[iS,1] = (fRef[iksi,:, 1]' * Vel[Fe.Glob[:,iF]])
      @views VelCa[iS,2] = (fRef[iksi,:, 2]' * Vel[Fe.Glob[:,iF]])
      @views VelCa[iS,3] = (fRef[iksi,:, 3]' * Vel[Fe.Glob[:,iF]])
    end  
  end
end

function ConvertVelocity!(backend,FTB,VelSp,Vel,Fe::VectorElement,Grid,Jacobi,ksi,::Grids.SphericalGrid)

  Numksi = size(ksi,1)
  fRef = zeros(Numksi,Fe.Comp,Fe.DoF)
  @inbounds for iksi = 1 : Numksi
    @inbounds for iComp = 1 : Fe.Comp
      @inbounds for iD = 1 : Fe.DoF
        fRef[iksi,iComp,iD] = Fe.phi[iD,iComp](ksi[iksi,1],ksi[iksi,2])
      end
    end
  end

  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  iS = 0
  VelCa = zeros(3)
  VelSpLoc = zeros(3)
  @inbounds for iF = 1 : Grid.NumFaces
    @inbounds for iksi = 1 : Numksi
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,ksi[iksi,1],ksi[iksi,2],Grid.Faces[iF],Grid)
      @views VelCa[1] = (fRef[iksi, 1,:]' * Vel[Fe.Glob[:,iF]])
      @views VelCa[2] = (fRef[iksi, 2,:]' * Vel[Fe.Glob[:,iF]])
      @views VelCa[3] = (fRef[iksi, 3,:]' * Vel[Fe.Glob[:,iF]])
      lon,lat,_ = Grids.cart2sphere(X[1],X[2],X[3])
      VelSpLoc = Grids.VelCart2Sphere(VelCa,lon,lat)
      iS += 1
      VelSp[iS,1] = VelSpLoc[1]
      VelSp[iS,2] = VelSpLoc[2]
    end
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
  @inbounds for iComp = 1 : Fe.Comp
    @inbounds for iD = 1 : Fe.DoF
      fRef[iComp,iD] = Fe.phi[iD,iComp](ksi1,ksi2)
    end
  end

  VelLoc = zeros(Fe.DoF)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  @inbounds for iF = 1 : Grid.NumFaces
    Jacobi(DF,detDF,pinvDF,X,Grid.Type,ksi1,ksi2,Grid.Faces[iF],Grid)
    detDFLoc = detDF[1]
    @views VelCa[iF,:]  .= pinvDF * (fRef[:, :] * Vel[Fe.Glob[:,iF]])
  end
end

function ConvertVelocitySp!(backend,FTB,VelSp,Vel,Fe::HDivElement,Grid,Jacobi,ksi)
  
  Numksi = size(ksi,1)    
  fRef = zeros(Numksi,Fe.Comp,Fe.DoF)
  @inbounds for iksi = 1 : Numksi
    @inbounds for iComp = 1 : Fe.Comp
      @inbounds for iD = 1 : Fe.DoF
        fRef[iksi,iComp,iD] = Fe.phi[iD,iComp](ksi[iksi,1],ksi[iksi,2])
      end
    end
  end

  VelLoc = zeros(Fe.DoF)
  VelCa = zeros(3)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  iS = 0
  @inbounds for iF = 1 : Grid.NumFaces
    @. VelLoc = Vel[Fe.Glob[:,iF]]
    @inbounds for iksi = 1 : Numksi
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,ksi[iksi,1],ksi[iksi,2],Grid.Faces[iF],Grid)
      detDFLoc = detDF[1]
      VelCa .= (1.0 / detDFLoc) * DF * (fRef[iksi,:, :] * VelLoc) * Grid.Faces[iF].Orientation
      lon,lat,_ = Grids.cart2sphere(X[1],X[2],X[3])
      VelSpLoc = Grids.VelCart2Sphere(VelCa,lon,lat)
      iS += 1
      VelSp[iS,1] = VelSpLoc[1]
      VelSp[iS,2] = VelSpLoc[2]
    end  
  end  
end
  
   # VelCa .= 1.0 / detDFLoc * DF * (fRef[:, :] * VelLoc) 
   

function ConvertScalarVelocityCart!(backend,FTB,VelCart,Vel,Fe::HDivElement,h,hFe::ScalarElement,
  Grid,Jacobi,ksi)
      
  Numksi = size(ksi,1)    
  fRef = zeros(Numksi,Fe.Comp,Fe.DoF)
  hfRef = zeros(Numksi,hFe.Comp,hFe.DoF)
  @inbounds for iksi = 1 : Numksi
    @inbounds for iComp = 1 : Fe.Comp
      @inbounds for iD = 1 : Fe.DoF
        fRef[iksi,iComp,iD] = Fe.phi[iD,iComp](ksi[iksi,1],ksi[iksi,2])
      end
    end
    @inbounds for iComp = 1 : hFe.Comp
      @inbounds for iD = 1 : hFe.DoF
        hfRef[iksi,iComp,iD] = hFe.phi[iD,iComp](ksi[iksi,1],ksi[iksi,2])
      end
    end
  end

  VelLoc = zeros(Fe.DoF)
  VelCa = zeros(3)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  iS = 0
  @inbounds for iF = 1 : Grid.NumFaces
    VelLoc = Vel[Fe.Glob[:,iF]]  
    hLoc = h[hFe.Glob[:,iF]]
    @inbounds for iksi = 1 : Numksi
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,ksi[iksi,1],ksi[iksi,2],Grid.Faces[iF],Grid)
      detDFLoc = detDF[1]
      VelCa .= (Grid.Faces[iF].Orientation / detDFLoc) * DF * (fRef[iksi,:, :] * VelLoc) / (hfRef[iksi,:, :] * hLoc)
      iS += 1
      VelCart[iS,1] = VelCa[1]
      VelCart[iS,2] = VelCa[2]
      VelCart[iS,3] = VelCa[3]
    end  
  end  
end

function ConvertVelocityCart!(backend,FTB,VelCart,Vel,Fe::HDivElement,
  Grid,Jacobi,ksi)
      
  Numksi = size(ksi,1)    
  fRef = zeros(Numksi,Fe.Comp,Fe.DoF)
  @inbounds for iksi = 1 : Numksi
    @inbounds for iComp = 1 : Fe.Comp
      @inbounds for iD = 1 : Fe.DoF
        fRef[iksi,iComp,iD] = Fe.phi[iD,iComp](ksi[iksi,1],ksi[iksi,2])
      end
    end
  end

  VelLoc = zeros(Fe.DoF)
  VelCa = zeros(3)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  iS = 0
  @inbounds for iF = 1 : Grid.NumFaces
    VelLoc = Vel[Fe.Glob[:,iF]]  
    @inbounds for iksi = 1 : Numksi
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,ksi[iksi,1],ksi[iksi,2],Grid.Faces[iF],Grid)
      detDFLoc = detDF[1]
      VelCa .= (Grid.Faces[iF].Orientation / detDFLoc) * DF * (fRef[iksi,:, :] * VelLoc) 
      iS += 1
      VelCart[iS,1] = VelCa[1]
      VelCart[iS,2] = VelCa[2]
      VelCart[iS,3] = VelCa[3]
    end  
  end  
end

function ConvertScalarVelocity!(backend,FTB,VelCart,Vel,Fe::HDivElement,h,hFe::ScalarElement,
  Grid,Jacobi,::Grids.CartesianGrid)
  if Grid.Type == Grids.Tri()
    ksi1 = -1/3 
    ksi2 = -1/3
  elseif Grid.Type == Grids.Quad()
    ksi1 = 0
    ksi2 = 0
  end  
      
  fRef = zeros(Fe.Comp,Fe.DoF)
  hfRef = zeros(hFe.Comp,hFe.DoF)
  @inbounds for iComp = 1 : Fe.Comp
    @inbounds for iD = 1 : Fe.DoF
      fRef[iComp,iD] = Fe.phi[iD,iComp](ksi1,ksi2)
    end
  end
  @inbounds for iComp = 1 : hFe.Comp
    @inbounds for iD = 1 : hFe.DoF
      hfRef[iComp,iD] = hFe.phi[iD,iComp](ksi1,ksi2)
    end
  end

  VelLoc = zeros(Fe.DoF)
  VelCa = zeros(3)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  @inbounds for iF = 1 : Grid.NumFaces
    VelLoc = Vel[Fe.Glob[:,iF]]  
    hLoc = h[hFe.Glob[:,iF]]
    Jacobi(DF,detDF,pinvDF,X,Grid.Type,ksi1,ksi2,Grid.Faces[iF],Grid)
    detDFLoc = detDF[1]
    VelCa .= (Grid.Faces[iF].Orientation / detDFLoc) * DF * (fRef[:, :] * VelLoc) / (hfRef[:, :] * hLoc)
    VelCart[iF,1] = VelCa[1]
    VelCart[iF,2] = VelCa[2]
    VelCart[iF,3] = VelCa[3]
  end  
end

function ConvertScalarVelocity!(backend,FTB,VelSp,Vel,Fe::HDivElement,h,hFe::ScalarElement,
  Grid,Jacobi,::Grids.SphericalGrid)
  if Grid.Type == Grids.Tri()
    ksi1 = -1/3 
    ksi2 = -1/3
  elseif Grid.Type == Grids.Quad()
    ksi1 = 0
    ksi2 = 0
  end  
      
  fRef = zeros(Fe.Comp,Fe.DoF)
  hfRef = zeros(hFe.Comp,hFe.DoF)
  @inbounds for iComp = 1 : Fe.Comp
    @inbounds for iD = 1 : Fe.DoF
      fRef[iComp,iD] = Fe.phi[iD,iComp](ksi1,ksi2)
    end
  end
  @inbounds for iComp = 1 : hFe.Comp
    @inbounds for iD = 1 : hFe.DoF
      hfRef[iComp,iD] = hFe.phi[iD,iComp](ksi1,ksi2)
    end
  end

  VelLoc = zeros(Fe.DoF)
  VelCa = zeros(3)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  @inbounds for iF = 1 : Grid.NumFaces
    VelLoc = Vel[Fe.Glob[:,iF]]  
    hLoc = h[hFe.Glob[:,iF]]
    Jacobi(DF,detDF,pinvDF,X,Grid.Type,ksi1,ksi2,Grid.Faces[iF],Grid)
    detDFLoc = detDF[1]
    VelCa .= (Grid.Faces[iF].Orientation / detDFLoc) * DF * (fRef[:, :] * VelLoc) / (hfRef[:, :] * hLoc)
    lon,lat,_ = Grids.cart2sphere(X[1],X[2],X[3])
    VelSpLoc = Grids.VelCart2Sphere(VelCa,lon,lat)
    VelSp[iF,1] = VelSpLoc[1]
    VelSp[iF,2] = VelSpLoc[2]
  end  
end

function ConvertScalarVelocitySp!(backend,FTB,VelSp,Vel,Fe::HDivElement,h,hFe::ScalarElement,Grid,Jacobi,ksi)

  Numksi = size(ksi,1)    
  fRef = zeros(Numksi,Fe.Comp,Fe.DoF)
  hfRef = zeros(Numksi,hFe.Comp,hFe.DoF)
  @inbounds for iksi = 1 : Numksi
    @inbounds for iComp = 1 : Fe.Comp
      @inbounds for iD = 1 : Fe.DoF
        fRef[iksi,iComp,iD] = Fe.phi[iD,iComp](ksi[iksi,1],ksi[iksi,2])
      end 
    end 
    @inbounds for iComp = 1 : hFe.Comp
      @inbounds for iD = 1 : hFe.DoF
        hfRef[iksi,iComp,iD] = hFe.phi[iD,iComp](ksi[iksi,1],ksi[iksi,2])
      end 
    end 
  end 
      
  VelLoc = zeros(Fe.DoF)
  VelCa = zeros(3)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  iS = 0
  @inbounds for iF = 1 : Grid.NumFaces
    VelLoc = Vel[Fe.Glob[:,iF]]  
    hLoc = h[hFe.Glob[:,iF]]
    @inbounds for iksi = 1 : Numksi
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,ksi[iksi,1],ksi[iksi,2],Grid.Faces[iF],Grid)
      detDFLoc = detDF[1]
      VelCa .= (Grid.Faces[iF].Orientation / detDFLoc) * DF * (fRef[iksi,:, :] * VelLoc) / (hfRef[iksi,:, :] * hLoc)
      lon,lat,_ = Grids.cart2sphere(X[1],X[2],X[3])
      VelSpLoc = Grids.VelCart2Sphere(VelCa,lon,lat)
      iS += 1
      VelSp[iS,1] = VelSpLoc[1]
      VelSp[iS,2] = VelSpLoc[2]
    end  
  end  
end

function ConvertVelocity!(backend,FTB,VelSp,Vel,Fe::VectorElement,Grid,Jacobi,::Grids.SphericalGrid)
  if Grid.Type == Grids.Tri()
    ksi1 = -1/3 
    ksi2 = -1/3
  elseif Grid.Type == Grids.Quad()
    ksi1 = 0
    ksi2 = 0
  end  
      
  fRef = zeros(Fe.Comp,Fe.DoF)
  @inbounds for iComp = 1 : Fe.Comp
    @inbounds for iD = 1 : Fe.DoF
      fRef[iComp,iD] = Fe.phi[iD,iComp](ksi1,ksi2)
    end
  end

  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  VelCa = zeros(3)
  X = zeros(3)
  @inbounds for iF = 1 : Grid.NumFaces
    VelCa .= (fRef[:, :] * Vel[Fe.Glob[:,iF]])
    Jacobi(DF,detDF,pinvDF,X,Grid.Type,ksi1,ksi2,Grid.Faces[iF],Grid)
    lon,lat,_ = Grids.cart2sphere(X[1],X[2],X[3])
    VelSpLoc = Grids.VelCart2Sphere(VelCa,lon,lat)
    VelSp[iF,1] = VelSpLoc[1]
    VelSp[iF,2] = VelSpLoc[2]
  end  
end

function ConvertVelocitySp!(backend,FTB,VelSp,Vel,Fe::HCurlElement,Grid,Jacobi,ksi)

  Numksi = size(ksi,1)    
  fRef = zeros(Numksi,Fe.Comp,Fe.DoF)
  @inbounds for iksi = 1 : Numksi
    @inbounds for iComp = 1 : Fe.Comp
      @inbounds for iD = 1 : Fe.DoF
        fRef[iksi,iComp,iD] = Fe.phi[iD,iComp](ksi[iksi,1],ksi[iksi,2])
      end
    end
  end

  VelLoc = zeros(Fe.DoF)
  VelCa = zeros(3)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  iS = 0
  @inbounds for iF = 1 : Grid.NumFaces
    VelLoc = Vel[Fe.Glob[:,iF]]
    @inbounds for iksi = 1 : Numksi
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,ksi[iksi,1],ksi[iksi,2],Grid.Faces[iF],Grid)
      detDFLoc = detDF[1]
      VelCa .= pinvDF * (fRef[iksi,:, :] * VelLoc)
      lon,lat,_ = Grids.cart2sphere(X[1],X[2],X[3])
      VelSpLoc = Grids.VelCart2Sphere(VelCa,lon,lat)
      iS += 1
      VelSp[iS,1] = VelSpLoc[1]
      VelSp[iS,2] = VelSpLoc[2]
    end  
  end  
end

function ConvertVelocityCart!(backend,FTB,VelCart,Vel,Fe::HCurlElement,Grid,Jacobi,ksi)

  Numksi = size(ksi,1)    
  fRef = zeros(Numksi,Fe.Comp,Fe.DoF)
  @inbounds for iksi = 1 : Numksi
    @inbounds for iComp = 1 : Fe.Comp
      @inbounds for iD = 1 : Fe.DoF
        fRef[iksi,iComp,iD] = Fe.phi[iD,iComp](ksi[iksi,1],ksi[iksi,2])
      end
    end
  end

  VelLoc = zeros(Fe.DoF)
  VelCa = zeros(3)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  iS = 0
  @inbounds for iF = 1 : Grid.NumFaces
    VelLoc = Vel[Fe.Glob[:,iF]]
    @inbounds for iksi = 1 : Numksi
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,ksi[iksi,1],ksi[iksi,2],Grid.Faces[iF],Grid)
      detDFLoc = detDF[1]
      VelCa .= pinvDF * (fRef[iksi,:, :] * VelLoc)
      iS += 1
      VelCart[iS,1] = VelCa[1]
      VelCart[iS,2] = VelCa[2]
      VelCart[iS,3] = VelCa[3]
    end  
  end  
end

