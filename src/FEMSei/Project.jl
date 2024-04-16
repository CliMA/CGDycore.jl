function Project(backend,FTB,Fe::ScalarElement,Grid,QuadOrd,Jacobi,F)
  QQ = FEMSei.QuadRule{FTB}(Grid.Type,backend,QuadOrd)
  Weights = QQ.Weights
  Points = QQ.Points
  fRef  = zeros(Fe.Comp,Fe.DoF,length(Weights))

  p=zeros(Fe.NumG)
  for i = 1 : length(Weights)
    for iComp = 1 : Fe.Comp
      for iD = 1 : Fe.DoF
        fRef[iComp,iD,i] = Fe.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  for iF = 1 : Grid.NumFaces
    pLoc = zeros(Fe.Comp,Fe.DoF)
    for i = 1 : length(Weights)
      _, detJ, _,X = Jacobi(Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      f, = F(X,0.0)
      pLoc += abs(detJ)*Weights[i]*(fRef[:,:,i]*f)
    end
    @. p[Fe.Glob[:,iF]] += pLoc[Fe.Comp,:]
  end
  p = Fe.M\p
  return p
end

function Project!(backend,FTB,p,Fe::ScalarElement,Grid,QuadOrd,Jacobi,F)
  QQ = FEMSei.QuadRule{FTB}(Grid.Type,backend,QuadOrd)
  Weights = QQ.Weights
  Points = QQ.Points
  fRef  = zeros(Fe.Comp,Fe.DoF,length(Weights))

  for i = 1 : length(Weights)
    for iComp = 1 : Fe.Comp
      for iD = 1 : Fe.DoF
        fRef[iComp,iD,i] = Fe.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  @. p = 0
  for iF = 1 : Grid.NumFaces
    pLoc = zeros(Fe.Comp,Fe.DoF)
    for i = 1 : length(Weights)
      _, detJ, _,X = Jacobi(Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      f, = F(X,0.0)
      pLoc += abs(detJ)*Weights[i]*(fRef[:,:,i]*f)
    end
    @. p[Fe.Glob[:,iF]] += pLoc[Fe.Comp,:]
  end
  p .= Fe.M \ p
end

function Project(backend,FTB,Fe::HDivKiteDElement,Grid,QuadOrd,Jacobi,F)
  QQ = FEMSei.QuadRule{FTB}(Grid.Type,backend,QuadOrd)
  Weights = QQ.Weights
  Points = QQ.Points
  fRef  = zeros(Fe.Comp,Fe.DoF,length(Weights))

  p=zeros(Fe.NumG)
  f = zeros(3)
  for i = 1 : length(Weights)
    for iComp = 1 : Fe.Comp
      for iD = 1 : Fe.DoF
        fRef[iComp,iD,i] = Fe.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  for iF = 1 : Grid.NumFaces
    pLoc = zeros(Fe.DoF)
    for i = 1 : length(Weights)
      DF, detJ,_,X = Jacobi(Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      _,f[1],f[2],f[3], = F(X,0.0)
      pLoc += sign(detJ)*Weights[i]*(fRef[:,:,i]' * (DF' * f))
    end
    @. p[Fe.Glob[:,iF]] += pLoc[:]
  end
  p = Fe.M\p
  return p
end

function Project!(backend,FTB,p,Fe::HDivKiteDElement,Grid,QuadOrd,Jacobi,F)
  QQ = FEMSei.QuadRule{FTB}(Grid.Type,backend,QuadOrd)
  Weights = QQ.Weights
  Points = QQ.Points
  fRef  = zeros(Fe.Comp,Fe.DoF,length(Weights))

  pp=zeros(Fe.NumG)
  VelSp = zeros(3)
  for i = 1 : length(Weights)
    for iComp = 1 : Fe.Comp
      for iD = 1 : Fe.DoF
        fRef[iComp,iD,i] = Fe.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  for iF = 1 : Grid.NumFaces
    pLoc = zeros(Fe.DoF)
    for i = 1 : length(Weights)
      DF, detJ,_,X = Jacobi(Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      _,VelSp[1],VelSp[2],VelSp[3], = F(X,0.0)
      lon,lat,r = Grids.cart2sphere(X[1],X[2],X[3])
      VelCa = VelSphere2Cart(VelSp,lon,lat)
      pLoc += sign(detJ)*Weights[i]*(fRef[:,:,i]' * (DF' * VelCa))
    end
    @. pp[Fe.Glob[:,iF]] += pLoc[:]
  end
  @show "HDiv",maximum(pp)
  pp = Fe.M \ pp
  @. p = pp
end

function Project!(backend,FTB,p,Fe::HCurlKiteDElement,Grid,QuadOrd,Jacobi,F)
  QQ = FEMSei.QuadRule{FTB}(Grid.Type,backend,QuadOrd)
  Weights = QQ.Weights
  Points = QQ.Points
  fRef  = zeros(Fe.Comp,Fe.DoF,length(Weights))

  pp=zeros(Fe.NumG)
  VelSp = zeros(3)
  for i = 1 : length(Weights)
    for iComp = 1 : Fe.Comp
      for iD = 1 : Fe.DoF
        fRef[iComp,iD,i] = Fe.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  for iF = 1 : Grid.NumFaces
    pLoc = zeros(Fe.DoF)
    for i = 1 : length(Weights)
      _, detJ,invPDF,X = Jacobi(Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      _,VelSp[1],VelSp[2],VelSp[3], = F(X,0.0)
      lon,lat,r = Grids.cart2sphere(X[1],X[2],X[3])
      VelCa = VelSphere2Cart(VelSp,lon,lat)
      @views pLoc .+= detJ * Weights[i] * (fRef[:,:,i]' * (invPDF' * VelCa))
    end
    @views @. pp[Fe.Glob[:,iF]] += pLoc[:]
  end
  @show "HCurl",maximum(pp)
  pp = Fe.M \ pp
  @. p = pp
end

function ProjectHDivHCurl!(backend,FTB,uCurl,Fe::HCurlKiteDElement,Grid,QuadOrd,Jacobi,
  FeF::HDivKiteDElement,uDiv)
  QQ = FEMSei.QuadRule{FTB}(Grid.Type,backend,QuadOrd)
  Weights = QQ.Weights
  Points = QQ.Points
  fRef  = zeros(Fe.Comp,Fe.DoF,length(Weights))
  fFRef  = zeros(FeF.Comp,FeF.DoF,length(Weights))

  pp=zeros(Fe.NumG)
  @inbounds for i = 1 : length(Weights)
    @inbounds for iComp = 1 : Fe.Comp
      @inbounds for iD = 1 : Fe.DoF
        fRef[iComp,iD,i] = Fe.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  @inbounds for i = 1 : length(Weights)
    @inbounds for iComp = 1 : Fe.Comp
      @inbounds for iD = 1 : Fe.DoF
        fFRef[iComp,iD,i] = FeF.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  pLoc = zeros(Fe.DoF)
  ppF = zeros(FeF.DoF)
  @inbounds for iF = 1 : Grid.NumFaces
    @. pLoc = 0
    @views ppF .= uDiv[FeF.Glob[:,iF]]
    @inbounds for i = 1 : length(Weights)
      @views pLoc += Weights[i] * ((fRef[:,:,i])' * (fFRef[:,:,i] * ppF))
    end
    @. pp[Fe.Glob[:,iF]] += pLoc[:]
  end
  pp = Fe.M \ pp
  @. uCurl = pp
end


function ComputeScalar(backend,FTB,Fe::ScalarElement,Grid,p)
  fRef  = zeros(Fe.Comp,Fe.DoF)

  for iComp = 1 : Fe.Comp
    for iD = 1 : Fe.DoF
      fRef[iComp,iD] = Fe.phi[iD,iComp](0.0,0.0)
    end
  end
  pM = zeros(Grid.NumFaces,1)
  for iF = 1 : Grid.NumFaces
    pLoc = p[Fe.Glob[:,iF]]
    pM[iF,:] = fRef[:,:]*pLoc
  end
  return pM
end

function ComputeVector(backend,FTB,Fe::HDivKiteDElement,Grid,Jacobi,p)
  fRef  = zeros(Fe.Comp,Fe.DoF)

  for iComp = 1 : Fe.Comp
    for iD = 1 : Fe.DoF
      fRef[iComp,iD] = Fe.phi[iD,iComp](0.0,0.0)
    end
  end

  pM = zeros(Grid.NumFaces,3)
  for iF = 1 : Grid.NumFaces
    DF, detJ,_,X = Jacobi(Grid.Type,0.0,0.0,Grid.Faces[iF], Grid)
    pLoc = p[Fe.Glob[:,iF]]
    pM[iF,:] = 1/detJ * DF * (fRef[:,:]*pLoc)
  end
  return pM
end
