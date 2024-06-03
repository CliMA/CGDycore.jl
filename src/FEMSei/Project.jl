function Project(backend,FTB,Fe::ScalarElement,Grid,QuadOrd,Jacobi,F)
  NumQuad,Weights,Points = QuadRule(Grid.Type,QuadOrd)
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

function Project!(backend,FTB,p,Fe::ScalarElement,Grid,QuadOrd,Jacobi,F)
  NumQuad,Weights,Points = QuadRule(Grid.Type,QuadOrd)
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
  ldiv!(Fe.LUM,p)
end


function Project!(backend,FTB,p,Fe::HDivElement,Grid,QuadOrd,Jacobi,F)
  NumQuad,Weights,Points = QuadRule(Grid.Type,QuadOrd)
  fRef  = zeros(Fe.Comp,Fe.DoF,length(Weights))

  @. p = 0
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
      pLoc += Grid.Faces[iF].Orientation * Weights[i] * (fRef[:,:,i]' * (DF' * VelCa))
    end
    @. p[Fe.Glob[:,iF]] += pLoc[:]
  end
  ldiv!(Fe.LUM,p)
end

function Project!(backend,FTB,p,Fe::HCurlElement,Grid,QuadOrd,Jacobi,F)
  NumQuad,Weights,Points = QuadRule(Grid.Type,QuadOrd)
  fRef  = zeros(Fe.Comp,Fe.DoF,length(Weights))

  @. p = 0
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
      @views pLoc .+= detDFLoc * Weights[i] * (fRef[:,:,i]' * (pinvDF' * VelCa))
    end
    @views @. p[Fe.Glob[:,iF]] += pLoc[:]
  end
  ldiv!(Fe.LUM,p)
end

function ProjectHDivHCurl!(backend,FTB,uCurl,Fe::HCurlElement,
  uDiv,FeF::HDivElement,Grid,ElemType::Grids.ElementType,QuadOrd,Jacobi)
  NumQuad,Weights,Points = QuadRule(ElemType,QuadOrd)
  fRef  = zeros(Fe.Comp,Fe.DoF,length(Weights))
  fFRef  = zeros(FeF.Comp,FeF.DoF,length(Weights))

  @. uCurl = 0
  @inbounds for i = 1 : length(Weights)
    @inbounds for iComp = 1 : Fe.Comp
      @inbounds for iD = 1 : Fe.DoF
        fRef[iComp,iD,i] = Fe.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  @inbounds for i = 1 : length(Weights)
    @inbounds for iComp = 1 : FeF.Comp
      @inbounds for iD = 1 : FeF.DoF
        fFRef[iComp,iD,i] = FeF.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  uCurlLoc = zeros(Fe.DoF)
  uuF = zeros(FeF.DoF)

  @inbounds for iF = 1 : Grid.NumFaces
    @. uCurlLoc = 0
    for iDoF = 1 : FeF.DoF
      ind = FeF.Glob[iDoF,iF]  
      uuF[iDoF] = uDiv[ind]
    end  
    for i = 1 : length(Weights)
      fFRefLoc1 = 0.0
      fFRefLoc2 = 0.0
      for iDoF = 1 : FeF.DoF
        fFRefLoc1 += fFRef[1,iDoF,i] * uuF[iDoF]  
        fFRefLoc2 += fFRef[2,iDoF,i] * uuF[iDoF]  
      end  
      for iDoF = 1 : Fe.DoF
        uCurlLoc[iDoF] += Grid.Faces[iF].Orientation * Weights[i] * (fRef[1,iDoF,i] * fFRefLoc1 +
          fRef[2,iDoF,i] * fFRefLoc2)
      end  
    end
    for iDoF = 1 : Fe.DoF
      ind = Fe.Glob[iDoF,iF]  
      uCurl[ind] += uCurlLoc[iDoF]
    end  
  end
  ldiv!(Fe.LUM,uCurl)
end

function ProjecthScalaruHDivHDiv!(backend,FTB,huDiv,Fe::HDivElement,
  h,hFeF::ScalarElement,uDiv,uFeF::HDivElement,Grid,ElemType::Grids.ElementType,QuadOrd,Jacobi)
  NumQuad,Weights,Points = QuadRule(ElemType,QuadOrd)
  fRef  = zeros(Fe.Comp,Fe.DoF,length(Weights))
  ufFRef  = zeros(uFeF.Comp,uFeF.DoF,length(Weights))
  hfFRef  = zeros(hFeF.Comp,hFeF.DoF,length(Weights))

  @. huDiv = 0
  @inbounds for i = 1 : length(Weights)
    @inbounds for iComp = 1 : Fe.Comp
      @inbounds for iD = 1 : Fe.DoF
        fRef[iComp,iD,i] = Fe.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  @inbounds for i = 1 : length(Weights)
    @inbounds for iComp = 1 : hFeF.Comp
      @inbounds for iD = 1 : hFeF.DoF
        hfFRef[iComp,iD,i] = hFeF.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  @inbounds for i = 1 : length(Weights)
    for iComp = 1 : uFeF.Comp
      for iD = 1 : uFeF.DoF
        ufFRef[iComp,iD,i] = uFeF.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  huDivLoc = zeros(Fe.DoF)
  uuF = zeros(uFeF.DoF)
  hhF = zeros(hFeF.DoF)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  @inbounds for iF = 1 : Grid.NumFaces
    @. huDivLoc = 0
    for iDoF = 1 : uFeF.DoF
      ind = uFeF.Glob[iDoF,iF]  
      uuF[iDoF] = uDiv[ind]
    end  
    for iDoF = 1 : hFeF.DoF
      ind = hFeF.Glob[iDoF,iF]  
      hhF[iDoF] = h[ind]
    end  
    for iQ = 1 : length(Weights)
      ufFRefLoc1 = 0.0
      ufFRefLoc2 = 0.0
      for iDoF = 1 : uFeF.DoF
        ufFRefLoc1 += ufFRef[1,iDoF,iQ] * uuF[iDoF]  
        ufFRefLoc2 += ufFRef[2,iDoF,iQ] * uuF[iDoF]  
      end  
      hfFRefLoc = 0.0
      for iDoF = 1 : hFeF.DoF
        hfFRefLoc += hfFRef[1,iDoF,iQ] * hhF[iDoF]  
      end  
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[iQ,1],Points[iQ,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      uLoc31 = DF[1,1] * ufFRefLoc1 + DF[1,2] * ufFRefLoc2
      uLoc32 = DF[2,1] * ufFRefLoc1 + DF[2,2] * ufFRefLoc2
      uLoc33 = DF[3,1] * ufFRefLoc1 + DF[3,2] * ufFRefLoc2
      uLoc1 = DF[1,1] * uLoc31 + DF[2,1] * uLoc32 + DF[3,1] * uLoc33
      uLoc2 = DF[1,2] * uLoc31 + DF[2,2] * uLoc32 + DF[3,2] * uLoc33
      for iDoF = 1 : Fe.DoF
        huDivLoc[iDoF] += Weights[iQ] * hfFRefLoc * (fRef[1,iDoF,iQ] * uLoc1 +
          fRef[2,iDoF,iQ] * uLoc2) / abs(detDFLoc)
      end  
    end
    for iDoF = 1 : Fe.DoF
      ind = Fe.Glob[iDoF,iF]  
      huDiv[ind] += huDivLoc[iDoF]
    end  
  end
  ldiv!(Fe.LUM,huDiv)
end

#=
function ProjectHDivHVort!(backend,FTB,qVort, u,FeF::HDivKiteDElement,
  FeT::ScalarElement,Grid,ElemType::Grids.ElementType,QuadOrd,Jacobi)
  NumQuad,Weights,Points = QuadRule(Grid.Type,QuadOrd)
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
  NumQuad,Weights,Points = QuadRule(ElemType,QuadOrd)
  fRef  = zeros(Fe.Comp,Fe.DoF,NumQuad)
  uFFRef  = zeros(uFeF.Comp,uFeF.DoF,NumQuad)
  hFFRef  = zeros(hFeF.Comp,hFeF.DoF,NumQuad)

  @. uh = 0
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
      uh[ind] += pLoc[iDoFFe]
    end  
  end
  ldiv!(Fe.LUM,uh)
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
