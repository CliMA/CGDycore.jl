function Interpolate!(backend,FTB,uN,Fe::HDivConfElement,Grid,QuadOrd,Jacobi,F)
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  X = zeros(3)
  VelSp = zeros(3)
  for iE = 1 : Grid.NumEdges
    X[1] = Grid.Edges[iE].Mid.x  
    X[2] = Grid.Edges[iE].Mid.y  
    X[3] = Grid.Edges[iE].Mid.z  
    _,VelSp[1],VelSp[2],VelSp[3], = F(X,0.0)
    lon,lat,r = Grids.cart2sphere(X[1],X[2],X[3])
    VelCa = VelSphere2Cart(VelSp,lon,lat)
    uN[iE] = Grid.Edges[iE].a * (Grid.Edges[iE].n.x * VelCa[1] +   
      Grid.Edges[iE].n.y * VelCa[2] + Grid.Edges[iE].n.z * VelCa[3])  
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
    uN[iE] = Grid.Edges[iE].a * (Grid.Edges[iE].n.x * VelCa[1] +   
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

function Project!(backend,FTB,p,Fe::ScalarElement,Grid,QuadOrd,Jacobi,F)
  NumQuad,Weights,Points = QuadRule(Fe.Type,QuadOrd)
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

function ProjectTr!(backend,FTB,p,Fe::ScalarElement,Grid,QuadOrd,Jacobi,F)
  NumQuad,Weights,Points = QuadRule(Fe.Type,QuadOrd)
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
      _,_,_,_,f = F(X,0.0)
      pLoc += abs(detDFLoc)*Weights[i]*(fRef[:,:,i]*f)
    end
    @. p[Fe.Glob[:,iF]] += pLoc[Fe.Comp,:]
  end
  ldiv!(Fe.LUM,p)
  @show size(p), size(Fe.LUM)
end

function Project!(backend,FTB,p,Fe::HDivElement,Grid,QuadOrd,Jacobi,F)
  NumQuad,Weights,Points = QuadRule(Fe.Type,QuadOrd)
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
      for iD = 1 : Fe.DoF
        pLoc[iD] += Grid.Faces[iF].Orientation * Weights[i] * (fRef[:,iD,i]' * (DF' * VelCa))
      end  
    end
    @views @. p[Fe.Glob[:,iF]] += pLoc
  end
  ldiv!(Fe.LUM,p)
  @show size(Fe.LUM), size(p), "Project!"
end

function Project!(backend,FTB,p,Fe::HCurlElement,Grid,QuadOrd,Jacobi,F)
  NumQuad,Weights,Points = QuadRule(Fe.Type,QuadOrd)
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

#projection of scalar to scalar, i.e. CG1->DG1
function ProjectScalarScalar!(backend,FTB,cP,FeP::ScalarElement,c,Fe::ScalarElement,Grid,ElemType::Grids.ElementType,QuadOrd,Jacobi)
  NumQuad,Weights,Points = QuadRule(ElemType,QuadOrd)
  fRef  = zeros(Fe.Comp,Fe.DoF,length(Weights))
  fPRef  = zeros(FeP.Comp,FeP.DoF,length(Weights))

  @. cP = 0
  @inbounds for i = 1 : length(Weights)
    @inbounds for iComp = 1 : Fe.Comp
      @inbounds for iD = 1 : Fe.DoF
        fRef[iComp,iD,i] = Fe.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  @inbounds for i = 1 : length(Weights)
    @inbounds for iComp = 1 : FeP.Comp
      @inbounds for iD = 1 : FeP.DoF
        fPRef[iComp,iD,i] = FeP.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  cPLoc = zeros(FeP.DoF)
  cc = zeros(Fe.DoF)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)


  @inbounds for iF = 1 : Grid.NumFaces
    @. cPLoc = 0
    for iDoF = 1 : Fe.DoF
      ind = Fe.Glob[iDoF,iF]  
      cc[iDoF] = c[ind]
    end  
    for iQ = 1 : length(Weights)
      fRefLoc = 0.0
      for iDoF = 1 : Fe.DoF
        fRefLoc += fRef[1,iDoF,iQ] * cc[iDoF]  
      end  
      #determinant
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[iQ,1],Points[iQ,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      for iDoF = 1 : FeP.DoF
        cPLoc[iDoF] +=  detDFLoc * Weights[iQ] * (fPRef[1,iDoF,iQ] * fRefLoc)
      end  
    end
    for iDoF = 1 : FeP.DoF
      ind = FeP.Glob[iDoF,iF]  
      cP[ind] += cPLoc[iDoF]
    end  
  end
  ldiv!(FeP.LUM,cP)
end

function ProjectHDivVecDG1!(backend,FTB,cP,FeP::VectorElement,c,Fe::HDivElement,Grid,
  ElemType::Grids.ElementType,QuadOrd,Jacobi)
  @. cP = 0
  fRef  = zeros(Fe.Comp,Fe.DoF)
  @inbounds for iComp = 1 : Fe.Comp
    @inbounds for iD = 1 : Fe.DoF
      fRef[iComp,iD] = Fe.phi[iD,iComp](0.0,0.0)
    end
  end
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  cPLoc = zeros(FeP.DoF)
  fRefLoc = zeros(2)
  @inbounds for iF = 1 : Grid.NumFaces
    @. fRefLoc = 0
    @inbounds for iDoF = 1 : Fe.DoF
      ind = Fe.Glob[iDoF,iF]  
      @views fRefLoc += fRef[:,iDoF] * c[ind]
    end  
    Jacobi!(DF,detDF,pinvDF,X,Grid.Type,0.0,0.0,Grid.Faces[iF], Grid)
    cPLoc = 1 / detDF[1] * DF * fRefLoc
    @inbounds for iDoF = 1 : FeP.DoF
      ind = FeP.Glob[iDoF,iF]  
      cP[ind] += cPLoc[iDoF]
    end  
  end
end

function ProjectScalarHDivVecDG1!(backend,FTB,uP,uFeP::VectorElement,h,hFe::ScalarElement,u,uFe::HDivElement,Grid,
  ElemType::Grids.ElementType,QuadOrd,Jacobi)
  @. uP = 0
  ufRef  = zeros(uFe.Comp,uFe.DoF)
  hfRef  = zeros(hFe.Comp,hFe.DoF)
  @inbounds for iComp = 1 : uFe.Comp
    @inbounds for iD = 1 : uFe.DoF
      ufRef[iComp,iD] = uFe.phi[iD,iComp](0.0,0.0)
    end
  end
  @inbounds for iD = 1 : hFe.DoF
    hfRef[1,iD] = hFe.phi[iD,1](0.0,0.0)
  end

  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  cPLoc = zeros(uFeP.DoF)
  ufRefLoc = zeros(2)
  @inbounds for iF = 1 : Grid.NumFaces
    @. ufRefLoc = 0
    @inbounds for iDoF = 1 : uFe.DoF
      ind = uFe.Glob[iDoF,iF]  
      @views ufRefLoc += ufRef[:,iDoF] * u[ind]
    end  
    hfRefLoc = 0
    @inbounds for iDoF = 1 : hFe.DoF
      ind = hFe.Glob[iDoF,iF]  
      hfRefLoc += hfRef[1,iDoF] * h[ind]
    end  
    Jacobi!(DF,detDF,pinvDF,X,Grid.Type,0.0,0.0,Grid.Faces[iF], Grid)
    uPLoc = 1 / detDF[1] * (DF * ufRefLoc)
    uPLoc = 1 / hfRefLoc
    @inbounds for iDoF = 1 : uFeP.DoF
      ind = uFeP.Glob[iDoF,iF]  
      uP[ind] += uPLoc[iDoF]
    end  
  end
end

function ProjectHDivVecDG!(backend,FTB,cP,FeP::VectorElement,c,Fe::HDivElement,Grid,ElemType::Grids.ElementType,QuadOrd,Jacobi)
  NumQuad,Weights,Points = QuadRule(ElemType,QuadOrd)
  fRef  = zeros(Fe.Comp,Fe.DoF,NumQuad)
  fPRef  = zeros(FeP.Comp,FeP.DoF,NumQuad)

  @. cP = 0
  @inbounds for iQ = 1 : NumQuad
    @inbounds for iComp = 1 : Fe.Comp
      @inbounds for iD = 1 : Fe.DoF
        fRef[iComp,iD,iQ] = Fe.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
  end
  @inbounds for iQ = 1 : NumQuad
    @inbounds for iComp = 1 : FeP.Comp
      @inbounds for iD = 1 : FeP.DoF
        fPRef[iComp,iD,iQ] = FeP.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
  end
  cPLoc = zeros(FeP.DoF)
  cc = zeros(Fe.DoF)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  fRefLoc = zeros(2)
  DFfRefLoc = zeros(3)

  @inbounds for iF = 1 : Grid.NumFaces
    @. cPLoc = 0
    @inbounds for iDoF = 1 : Fe.DoF
      ind = Fe.Glob[iDoF,iF]  
      cc[iDoF] = c[ind]
    end  
    @inbounds for iQ = 1 : NumQuad
      @. fRefLoc = 0.0
      @inbounds for iDoF = 1 : Fe.DoF
        @inbounds for iComp = 1 : Fe.Comp
          fRefLoc[iComp] += fRef[iComp,iDoF,iQ] * cc[iDoF]  
        end  
      end  
      #determinant
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[iQ,1],Points[iQ,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      DFfRefLoc[1] = DF[1,1] * fRefLoc[1] + DF[1,2] * fRefLoc[2]
      DFfRefLoc[2] = DF[2,1] * fRefLoc[1] + DF[2,2] * fRefLoc[2]
      DFfRefLoc[3] = DF[3,1] * fRefLoc[1] + DF[3,2] * fRefLoc[2]
      #DF=3x2 
      @inbounds for iDoF = 1 : FeP.DoF
        cPLoc[iDoF] += Grid.Faces[iF].Orientation * Weights[iQ] * (DFfRefLoc[1] * fPRef[1,iDoF,iQ] +
          DFfRefLoc[2] * fPRef[2,iDoF,iQ] + DFfRefLoc[3] * fPRef[3,iDoF,iQ])
      end  
    end
    @inbounds for iDoF = 1 : FeP.DoF
      ind = FeP.Glob[iDoF,iF]  
      cP[ind] += cPLoc[iDoF]
    end  
  end
  ldiv!(FeP.LUM,cP)
end

function ProjectHDivHCurl!(backend,FTB,uCurl,Fe::HCurlElement,
  uDiv,FeF::HDivElement,Grid,ElemType::Grids.ElementType,QuadOrd,Jacobi)
  NumQuad,Weights,Points = QuadRule(Fe.Type,QuadOrd)
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
  @time @inbounds for iF = 1 : Grid.NumFaces
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
          fRef[2,iDoF,iQ] * uLoc2) / detDFLoc
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
  NumQuad,Weights,Points = QuadRule(Fe.Type,QuadOrd)
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

function ProjectVectorScalarVectorHDiv(backend,FTB,u,Fe::VectorElement,
  h,hFeF::ScalarElement,huDiv,huFeF::HDivElement,Grid,ElemType::Grids.ElementType,QuadOrd,Jacobi)
  NumQuad,Weights,Points = QuadRule(ElemType,QuadOrd)
  fRef  = zeros(Fe.Comp,Fe.DoF,length(Weights))
  hufFRef  = zeros(huFeF.Comp,huFeF.DoF,length(Weights))
  hfFRef  = zeros(hFeF.Comp,hFeF.DoF,length(Weights))

  @. u = 0
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
    for iComp = 1 : huFeF.Comp
      for iD = 1 : huFeF.DoF
        hufFRef[iComp,iD,i] = huFeF.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  uLoc = zeros(Fe.DoF)
  huuF = zeros(huFeF.DoF)
  hhF = zeros(hFeF.DoF)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  @time @inbounds for iF = 1 : Grid.NumFaces
    @. uLoc = 0
    for iDoF = 1 : huFeF.DoF
      ind = huFeF.Glob[iDoF,iF]  
      huuF[iDoF] = huDiv[ind]
    end  
    for iDoF = 1 : hFeF.DoF
      ind = hFeF.Glob[iDoF,iF]  
      hhF[iDoF] = h[ind]
    end  
    for iQ = 1 : length(Weights)
      hufFRefLoc1 = 0.0
      hufFRefLoc2 = 0.0
      for iDoF = 1 : huFeF.DoF
        hufFRefLoc1 += hufFRef[1,iDoF,iQ] * huuF[iDoF]  
        hufFRefLoc2 += hufFRef[2,iDoF,iQ] * huuF[iDoF]   
      end  
      hfFRefLoc = 0.0
      for iDoF = 1 : hFeF.DoF
        hfFRefLoc += hfFRef[1,iDoF,iQ] * hhF[iDoF]  
      end  
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[iQ,1],Points[iQ,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      huLoc1 = DF[1,1] * hufFRefLoc1 + DF[1,2] * hufFRefLoc2
      huLoc2 = DF[2,1] * hufFRefLoc1 + DF[2,2] * hufFRefLoc2
      huLoc3 = DF[3,1] * hufFRefLoc1 + DF[3,2] * hufFRefLoc2
      for iDoF = 1 : Fe.DoF
        uLoc[iDoF] += Weights[iQ] / hfFRefLoc * (fRef[1,iDoF,iQ] * huLoc1 +
          fRef[2,iDoF,iQ] * huLoc2 + fRef[3,iDoF,iQ] * huLoc3)
      end  
    end
    for iDoF = 1 : Fe.DoF
      ind = Fe.Glob[iDoF,iF]  
      u[ind] += uLoc[iDoF]
    end  
  end
  ldiv!(Fe.LUM,u)
end

function ProjectScalar!(backend,FTB,p,Fe::HDivElement,Grid,QuadOrd,Jacobi,F)
  NumQuad,Weights,Points = QuadRule(Fe.Type,QuadOrd)
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
      h,VelSp[1],VelSp[2],VelSp[3], = F(X,0.0)
      lon,lat,r = Grids.cart2sphere(X[1],X[2],X[3])
      VelCa = VelSphere2Cart(VelSp,lon,lat)
      pLoc += Grid.Faces[iF].Orientation * Weights[i] * h * (fRef[:,:,i]' * (DF' * VelCa))
    end
    @. p[Fe.Glob[:,iF]] += pLoc[:]
  end
  ldiv!(Fe.LUM,p)
end
