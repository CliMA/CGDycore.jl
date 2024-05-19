function Project(backend,FTB,Fe::ScalarElement,Grid,QuadOrd,Jacobi,F)
  NumQuad,Weights,Points = FEMSei.QuadRule(Grid.Type,QuadOrd)
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
  p = Fe.M\p
  return p
end

function Project!(backend,FTB,p,Fe::ScalarElement,Grid,QuadOrd,Jacobi,F)
  NumQuad,Weights,Points = FEMSei.QuadRule(Grid.Type,QuadOrd)
  fRef  = zeros(Fe.Comp,Fe.DoF,length(Weights))

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
  @. p = 0
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
  p .= Fe.M \ p
end

function Project(backend,FTB,Fe::HDivKiteDElement,Grid,QuadOrd,Jacobi,F)
  NumQuad,Weights,Points = FEMSei.QuadRule(Grid.Type,QuadOrd)
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
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  for iF = 1 : Grid.NumFaces
    pLoc = zeros(Fe.DoF)
    for i = 1 : length(Weights)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      _,f[1],f[2],f[3], = F(X,0.0)
      pLoc += sign(detDFLoc)*Weights[i]*(fRef[:,:,i]' * (DF' * f))
    end
    @. p[Fe.Glob[:,iF]] += pLoc[:]
  end
  p = Fe.M\p
  return p
end

function Project!(backend,FTB,p,Fe::HDivElement,Grid,QuadOrd,Jacobi,F)
  NumQuad,Weights,Points = FEMSei.QuadRule(Grid.Type,QuadOrd)
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
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  for iF = 1 : Grid.NumFaces
    pLoc = zeros(Fe.DoF)
    for i = 1 : length(Weights)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      _,VelSp[1],VelSp[2],VelSp[3], = F(X,0.0)
      lon,lat,r = Grids.cart2sphere(X[1],X[2],X[3])
      VelCa = VelSphere2Cart(VelSp,lon,lat)
      pLoc += sign(detDFLoc)*Weights[i]*(fRef[:,:,i]' * (DF' * VelCa))
    end
    @. pp[Fe.Glob[:,iF]] += pLoc[:]
  end
  pp = Fe.M \ pp
  @. p = pp
end

function Project!(backend,FTB,p,Fe::HDivKiteDElement,Grid,QuadOrd,Jacobi,F)
  NumQuad,Weights,Points = FEMSei.QuadRule(Grid.Type,QuadOrd)
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
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  for iF = 1 : Grid.NumFaces
    pLoc = zeros(Fe.DoF)
    for i = 1 : length(Weights)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      _,VelSp[1],VelSp[2],VelSp[3], = F(X,0.0)
      lon,lat,r = Grids.cart2sphere(X[1],X[2],X[3])
      VelCa = VelSphere2Cart(VelSp,lon,lat)
      pLoc += sign(detDFLoc)*Weights[i]*(fRef[:,:,i]' * (DF' * VelCa))
    end
    @. pp[Fe.Glob[:,iF]] += pLoc[:]
  end
  pp = Fe.M \ pp
  @. p = pp
end


function Project!(backend,FTB,p,Fe::HCurlKiteDElement,Grid,QuadOrd,Jacobi,F)
  NumQuad,Weights,Points = FEMSei.QuadRule(Grid.Type,QuadOrd)
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
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  for iF = 1 : Grid.NumFaces
    pLoc = zeros(Fe.DoF)
    for i = 1 : length(Weights)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      _,VelSp[1],VelSp[2],VelSp[3], = F(X,0.0)
      lon,lat,r = Grids.cart2sphere(X[1],X[2],X[3])
      VelCa = VelSphere2Cart(VelSp,lon,lat)
      @views pLoc .+= detDFLoc * Weights[i] * (fRef[:,:,i]' * (invPDF' * VelCa))
    end
    @views @. pp[Fe.Glob[:,iF]] += pLoc[:]
  end
  pp = Fe.M \ pp
  @. p = pp
end

function ProjectHDivHCurl!(backend,FTB,uCurl,Fe::HCurlKiteDElement,
  uDiv,FeF::HDivKiteDElement,Grid,ElemType::Grids.ElementType,QuadOrd,Jacobi)
  NumQuad,Weights,Points = FEMSei.QuadRule(ElemType,QuadOrd)
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

#=
function ProjectHDivHVort!(backend,FTB,qVort, u,FeF::HDivKiteDElement,
  FeT::ScalarElement,Grid,ElemType::Grids.ElementType,QuadOrd,Jacobi)
  NumQuad,Weights,Points = FEMSei.QuadRule(Grid.Type,QuadOrd)
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
=#

function ProjectHDivHDivScalar!(backend,FTB,uh,Fe::HDivKiteDElement,
  u,uFeF::HDivKiteDElement,h,hFeF::ScalarElement,Grid,ElemType::Grids.ElementType,QuadOrd,Jacobi)
  NumQuad,Weights,Points = FEMSei.QuadRule(ElemType,QuadOrd)
  fRef  = zeros(Fe.Comp,Fe.DoF,NumQuad)
  uFFRef  = zeros(uFeF.Comp,uFeF.DoF,NumQuad)
  hFFRef  = zeros(hFeF.Comp,hFeF.DoF,NumQuad)

  pp=zeros(Fe.NumG)
  @inbounds for iQ = 1 : NumQuad
    @inbounds for iComp = 1 : Fe.Comp
      @inbounds for iD = 1 : Fe.DoF
        fRef[iComp,iD,iQ] = Fe.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    @inbounds for iComp = 1 : uFeF.Comp
      @inbounds for iD = 1 : uFeF.DoF
        uFFRef[iComp,iD,iQ] = uFeF.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    @inbounds for iComp = 1 : hFeF.Comp
      @inbounds for iD = 1 : hFeF.DoF
        hFFRef[iComp,iD,iQ] = hFeF.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
  end
  pLoc = zeros(Fe.DoF)
  uFLoc = zeros(2,NumQuad)
  uLoc = zeros(2)
  uLoc3 = zeros(3)
  hFLoc = zeros(NumQuad)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  @inbounds for iF = 1 : Grid.NumFaces
    @. uFLoc = 0
    @inbounds for iDoFuFeF = 1 : uFeF.DoF
      ind = uFeF.Glob[iDoFuFeF,iF]  
      @inbounds for iQ = 1 : NumQuad
        @views @. uFLoc[:,iQ] += uFFRef[:,iDoFuFeF,iQ] * u[ind]
      end  
    end  
    @. hFLoc = 0
    @inbounds for iDoFhFeF = 1 : hFeF.DoF
      ind = hFeF.Glob[iDoFhFeF,iF]  
      @inbounds for iQ = 1 : NumQuad
        hFLoc[iQ] += hFFRef[1,iDoFhFeF,iQ] * h[ind]
      end  
    end  
    @. pLoc = 0
    @inbounds for iQ = 1 : length(Weights)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[iQ,1],Points[iQ,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      uLoc31 = DF[1,1] * uFLoc[1,iQ] + DF[1,2] * uFLoc[2,iQ]
      uLoc32 = DF[2,1] * uFLoc[1,iQ] + DF[2,2] * uFLoc[2,iQ]
      uLoc33 = DF[3,1] * uFLoc[1,iQ] + DF[3,2] * uFLoc[2,iQ]
      uLoc1 = DF[1,1] * uLoc31 + DF[2,1] * uLoc32 + DF[3,1] * uLoc33
      uLoc2 = DF[1,2] * uLoc31 + DF[2,2] * uLoc32 + DF[3,2] * uLoc33
      uLoc1  = hFLoc[iQ] / detDFLoc * uLoc1
      uLoc2  = hFLoc[iQ] / detDFLoc * uLoc2
      for iDoFFe = 1 : Fe.DoF
        pLoc[iDoFFe] +=  uLoc1 * fRef[1,iDoFFe,iQ] + uLoc2 * fRef[2,iDoFFe,iQ]  
      end  
    end  
    for iDoFFe = 1 : Fe.DoF
      ind = Fe.Glob[iDoFFe,iF] 
      pp[ind] += pLoc[iDoFFe]
    end  
  end
  pp = Fe.M \ pp
  @. uh = pp
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

  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  pM = zeros(Grid.NumFaces,3)
  for iF = 1 : Grid.NumFaces
    Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
    detDFLoc = detDF[1]
    pLoc = p[Fe.Glob[:,iF]]
    pM[iF,:] = 1/detDFLoc * DF * (fRef[:,:]*pLoc)
  end
  return pM
end
