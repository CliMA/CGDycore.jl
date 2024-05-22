function VortCrossVel!(backend,FTB,Rhs,u,uFe::HDivElement,q,qFe::ScalarElement,FeT::HDivElement,
  Grid,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(Grid.Type,QuadOrd)
  fTRef  = zeros(FeT.Comp,FeT.DoF,length(Weights))
  qfFRef  = zeros(qFe.Comp,qFe.DoF,length(Weights))

  for i = 1 : length(Weights)
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeT.DoF
        fTRef[iComp,iD,i] = FeT.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
    for iComp = 1 : qFe.Comp
      for iD = 1 : qFe.DoF
        qfFRef[iComp,iD,i] = qFe.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  ufFRef = fTRef
  RhsLoc = zeros(FeT.DoF)
  uLoc = zeros(uFe.Comp)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  for iF = 1 : Grid.NumFaces
    RhsLoc .= 0
    for i = 1 : length(Weights)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      qLoc = qfFRef[1,1,i] * q[qFe.Glob[1,iF]]
      for j = 2 : qFe.DoF 
        qLoc += qfFRef[1,j,i] * q[qFe.Glob[j,iF]]
      end  
      uLoc .= [-ufFRef[2,1,i],ufFRef[1,1,i]] * u[uFe.Glob[1,iF]]
      for j = 2 : FeT.DoF 
        uLoc .+= [-ufFRef[2,j,i],ufFRef[1,j,i]] * u[uFe.Glob[j,iF]]
      end  
      RhsLoc += 1 / detDFLoc * Weights[i] * qLoc * ((DF * uLoc)' * (DF * fTRef[:,:,i]))'
    end
    @views @. Rhs[FeT.Glob[:,iF]] += RhsLoc
  end
end  

function DivMatrix(backend,FTB,FeF::HDivElement,FeT::ScalarElement,Grid,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(Grid.Type,QuadOrd)
  fFRef  = zeros(FeT.Comp,FeF.DoF,length(Weights))
  fTRef  = zeros(FeT.Comp,FeT.DoF,length(Weights))

  for i = 1 : length(Weights)
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeT.DoF
        fTRef[iComp,iD,i] = FeT.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeF.DoF
        fFRef[iComp,iD,i] = FeF.Divphi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  RowInd = Int64[]
  ColInd = Int64[]
  Val = Float64[]
  DivLoc = zeros(FeT.DoF,FeF.DoF)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)

  for iF = 1 : Grid.NumFaces
    DivLoc .= 0
    for i = 1 : length(Weights)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      DivLoc += sign(detDFLoc) * Weights[i] * (fTRef[:,:,i]' * fFRef[:,:,i])
    end
    for j = 1 : size(DivLoc,2)
      for i = 1 : size(DivLoc,1)
        push!(RowInd,FeT.Glob[i,iF])
        push!(ColInd,FeF.Glob[j,iF])
        push!(Val,DivLoc[i,j])
      end
    end
  end
  Div = sparse(RowInd, ColInd, Val)
  return Div
end

function CurlMatrix(backend,FTB,FeF::HCurlElement,FeT::ScalarElement,Grid,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(Grid.Type,QuadOrd)
  fFRef  = zeros(FeT.Comp,FeF.DoF,length(Weights))
  fTRef  = zeros(FeT.Comp,FeT.DoF,length(Weights))

  for i = 1 : length(Weights)
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeT.DoF
        fTRef[iComp,iD,i] = FeT.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeF.DoF
        fFRef[iComp,iD,i] = FeF.Curlphi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  RowInd = Int64[]
  ColInd = Int64[]
  Val = Float64[]
  CurlLoc = zeros(FeT.DoF,FeF.DoF)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)

  for iF = 1 : Grid.NumFaces
    CurlLoc .= 0
    for i = 1 : length(Weights)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      CurlLoc += sign(detDFLoc) * Weights[i] * (fTRef[:,:,i]' * fFRef[:,:,i])
    end
    for j = 1 : size(CurlLoc,2)
      for i = 1 : size(CurlLoc,1)
        push!(RowInd,FeT.Glob[i,iF])
        push!(ColInd,FeF.Glob[j,iF])
        push!(Val,CurlLoc[i,j])
      end
    end
  end
  Curl = sparse(RowInd, ColInd, Val)
  return Curl
end

function DivMatrix(backend,FTB,FeF::HDivKiteDElement,FeT::ScalarElement,Grid,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(Grid.Type,QuadOrd)
  fFRef  = zeros(FeT.Comp,FeF.DoF,length(Weights))
  fTRef  = zeros(FeT.Comp,FeT.DoF,length(Weights))

  for i = 1 : length(Weights)
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeT.DoF
        fTRef[iComp,iD,i] = FeT.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeF.DoF
        fFRef[iComp,iD,i] = FeF.Divphi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  fFRefX  = zeros(FeT.Comp,FeF.DoF,length(WeightsL))
  fFRefY  = zeros(FeT.Comp,FeF.DoF,length(WeightsL))
  fTRefX  = zeros(FeT.Comp,FeT.DoF,length(WeightsL))
  fTRefY  = zeros(FeT.Comp,FeT.DoF,length(WeightsL))
  for i = 1 : length(WeightsL)
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeT.DoF
        fTRefX[iComp,iD,i] = FeT.phi[iD,iComp](-1.0,PointsL[i])
        fTRefY[iComp,iD,i] = FeT.phi[iD,iComp](PointsL[i],-1.0)
      end
    end
    for iD = 1 : FeF.DoF
      fFRefX[1,iD,i] = FeF.phi[iD,1](-1.0,PointsL[i])
      fFRefY[1,iD,i] = FeF.phi[iD,2](PointsL[i],-1.0)
    end
  end

  RowInd = Int64[]
  ColInd = Int64[]
  Val = Float64[]
  DivLoc = zeros(FeT.DoF,FeF.DoF)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)

  for iF = 1 : Grid.NumFaces
    DivLoc .= 0
    for i = 1 : length(Weights)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      DivLoc += sign(detDFLoc) * Weights[i] * (fTRef[:,:,i]' * fFRef[:,:,i])
    end
    for i = 1 : length(WeightsL)
       DivLoc += WeightsL[i] * (fTRefX[:,:,i]' * fFRefX[:,:,i] +  
         fTRefY[:,:,i]' * fFRefY[:,:,i])  
    end   
    for j = 1 : size(DivLoc,2)
      for i = 1 : size(DivLoc,1)
        push!(RowInd,FeT.Glob[i,iF])
        push!(ColInd,FeF.Glob[j,iF])
        push!(Val,DivLoc[i,j])
      end
    end
  end
  Div = sparse(RowInd, ColInd, Val)
  return Div
end

function CurlMatrix(backend,FTB,FeF::HCurlKiteDElement,FeT::ScalarElement,Grid,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(Grid.Type,QuadOrd)
  fFRef  = zeros(FeF.Comp,FeF.DoF,length(Weights))
  fTRef  = zeros(FeF.Comp,FeT.DoF,length(Weights))

  for i = 1 : length(Weights)
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeT.DoF
        fTRef[iComp,iD,i] = FeT.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeF.DoF
        fFRef[iComp,iD,i] = FeF.Curlphi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  fFRefX  = zeros(FeT.Comp,FeF.DoF,length(WeightsL))
  fFRefY  = zeros(FeT.Comp,FeF.DoF,length(WeightsL))
  fTRefX  = zeros(FeT.Comp,FeT.DoF,length(WeightsL))
  fTRefY  = zeros(FeT.Comp,FeT.DoF,length(WeightsL))
  for i = 1 : length(WeightsL)
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeT.DoF
        fTRefX[iComp,iD,i] = FeT.phi[iD,iComp](-1.0,PointsL[i])
        fTRefY[iComp,iD,i] = FeT.phi[iD,iComp](PointsL[i],-1.0)
      end
    end
    for iD = 1 : FeF.DoF
      fFRefX[1,iD,i] = FeF.phi[iD,2](-1.0,PointsL[i])
      fFRefY[1,iD,i] = -FeF.phi[iD,1](PointsL[i],-1.0)
    end
  end

  RowInd = Int64[]
  ColInd = Int64[]
  Val = Float64[]
  CurlLoc = zeros(FeT.DoF,FeF.DoF)
  CurlLoc1 = zeros(FeT.DoF,FeF.DoF)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)

  for iF = 1 : Grid.NumFaces
    CurlLoc .= 0
    for i = 1 : length(Weights)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      CurlLoc += -sign(detDFLoc) * Weights[i] * (fTRef[:,:,i]' * fFRef[:,:,i])
    end
    for i = 1 : length(WeightsL)
      CurlLoc += WeightsL[i] * (fTRefX[:,:,i]' * fFRefX[:,:,i] +  
        fTRefY[:,:,i]' * fFRefY[:,:,i])  
    end   
    for j = 1 : size(CurlLoc,2)
      for i = 1 : size(CurlLoc,1)
        push!(RowInd,FeT.Glob[i,iF])
        push!(ColInd,FeF.Glob[j,iF])
        push!(Val,CurlLoc[i,j])
      end
    end
  end
  Curl = sparse(RowInd, ColInd, Val)
  return Curl
end

function GradMatrix(backend,FTB,FeF::ScalarElement,FeT::HDivElement,Grid,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(Grid.Type,QuadOrd)
  fFRef  = zeros(FeT.Comp,FeF.DoF,length(Weights))
  fTRef  = zeros(FeT.Comp,FeT.DoF,length(Weights))

  for i = 1 : length(Weights)
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeT.DoF
        fTRef[iComp,iD,i] = FeT.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeF.DoF
        fFRef[iComp,iD,i] = FeF.Gradphi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  RowInd = Int64[]
  ColInd = Int64[]
  Val = Float64[]
  GradLoc = zeros(FeT.DoF,FeF.DoF)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)

  for iF = 1 : Grid.NumFaces
    GradLoc .= 0
    for i = 1 : length(Weights)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      GradLoc += sign(detDFLoc) * Weights[i] * (fTRef[:,:,i]' * fFRef[:,:,i])
    end
    for j = 1 : size(GradLoc,2)
      for i = 1 : size(GradLoc,1)
        push!(RowInd,FeT.Glob[i,iF])
        push!(ColInd,FeF.Glob[j,iF])
        push!(Val,GradLoc[i,j])
      end
    end
  end
  Grad = sparse(RowInd, ColInd, Val)
  return Grad
end
    
function GradMatrix(backend,FTB,FeF::ScalarElement,FeT::HDivKiteDElement,Grid,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(Grid.Type,QuadOrd)
  fFRef  = zeros(FeT.Comp,FeF.DoF,length(Weights))
  fTRef  = zeros(FeT.Comp,FeT.DoF,length(Weights))

  for i = 1 : length(Weights)
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeT.DoF
        fTRef[iComp,iD,i] = FeT.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeF.DoF
        fFRef[iComp,iD,i] = FeF.Gradphi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  fTRefX  = zeros(FeF.Comp,FeT.DoF,length(WeightsL))
  fTRefY  = zeros(FeF.Comp,FeT.DoF,length(WeightsL))
  fFRefX  = zeros(FeF.Comp,FeF.DoF,length(WeightsL))
  fFRefY  = zeros(FeF.Comp,FeF.DoF,length(WeightsL))
  for i = 1 : length(WeightsL)
    for iComp = 1 : FeF.Comp
      for iD = 1 : FeF.DoF
        fFRefX[iComp,iD,i] = FeF.phi[iD,iComp](1.0,PointsL[i])
        fFRefY[iComp,iD,i] = FeF.phi[iD,iComp](PointsL[i],1.0)
      end
    end
    for iD = 1 : FeT.DoF
      fTRefX[1,iD,i] = FeT.phi[iD,1](1.0,PointsL[i])
      fTRefY[1,iD,i] = FeT.phi[iD,2](PointsL[i],1.0)
    end
  end
  RowInd = Int64[]
  ColInd = Int64[]
  Val = Float64[]
  GradLoc = zeros(FeT.DoF,FeF.DoF)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)

  for iF = 1 : Grid.NumFaces
    GradLoc .= 0
    for i = 1 : length(Weights)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      GradLoc += sign(detDFLoc) * Weights[i] * (fTRef[:,:,i]' * fFRef[:,:,i])
    end
    for i = 1 : length(WeightsL)
       GradLoc -= WeightsL[i] * (fTRefX[:,:,i]' * fFRefX[:,:,i] +  
         fTRefY[:,:,i]' * fFRefY[:,:,i])  
    end   
    for j = 1 : size(GradLoc,2)
      for i = 1 : size(GradLoc,1)
        push!(RowInd,FeT.Glob[i,iF])
        push!(ColInd,FeF.Glob[j,iF])
        push!(Val,GradLoc[i,j])
      end
    end
  end
  Grad = sparse(RowInd, ColInd, Val)
  return Grad
end

function GradKinHeight!(backend,FTB,Rhs,h,hFeF::ScalarElement,u,uFeF::HDivKiteDElement,
  FeT::HDivKiteDElement,Grid,ElemType::Grids.ElementType,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  hFRef  = zeros(FeT.Comp,hFeF.DoF,length(Weights))
  fTRef  = zeros(FeT.Comp,FeT.DoF,length(Weights))

  for iQ = 1 : NumQuad
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeT.DoF
        fTRef[iComp,iD,iQ] = FeT.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    for iComp = 1 : FeT.Comp
      for iD = 1 : hFeF.DoF
        hFRef[iComp,iD,iQ] = hFeF.Gradphi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
  end
  uFRef = fTRef
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  fTRefX  = zeros(FeT.Comp,FeT.DoF,NumQuadL)
  fTRefY  = zeros(FeT.Comp,FeT.DoF,NumQuadL)
  hFRefX  = zeros(hFeF.Comp,hFeF.DoF,NumQuadL)
  hFRefY  = zeros(hFeF.Comp,hFeF.DoF,NumQuadL)
  for iQ = 1 : NumQuadL
    for iComp = 1 : hFeF.Comp
      for iD = 1 : hFeF.DoF
        hFRefX[iComp,iD,iQ] = hFeF.phi[iD,iComp](1.0,PointsL[iQ])
        hFRefY[iComp,iD,iQ] = hFeF.phi[iD,iComp](PointsL[iQ],1.0)
      end
    end
    @show size(FeT.phi)
    for iD = 1 : FeT.DoF
      fTRefX[1,iD,iQ] = FeT.phi[iD,1](1.0,PointsL[iQ])
      fTRefX[2,iD,iQ] = FeT.phi[iD,2](1.0,PointsL[iQ])
      fTRefY[1,iD,iQ] = FeT.phi[iD,1](PointsL[iQ],1.0)
      fTRefY[2,iD,iQ] = FeT.phi[iD,2](PointsL[iQ],1.0)
    end
  end
  uFRefX = fTRefX
  uFRefY = fTRefY
  GradLoc = zeros(FeT.DoF)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)

  for iF = 1 : Grid.NumFaces
    GradLoc .= 0
    for iQ = 1 : length(Weights)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[iQ,1],Points[iQ,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      u1 = 0.0
      u2 = 0.0
      for iDoF = 1 : uFeF.DoF
        ind = uFeF.Glob[iDoF,iF]  
        @views u1 += uFRef[1,iDoF,iQ] * u[ind]
        @views u2 += uFRef[2,iDoF,iQ] * u[ind]
      end  
      uLoc1 = 1 / detDFLoc * (DF[1,1] * u1 + DF[1,2] * u2)
      uLoc2 = 1 / detDFLoc * (DF[2,1] * u1 + DF[2,2] * u2)
      uLoc3 = 1 / detDFLoc * (DF[3,1] * u1 + DF[3,2] * u2)
      hLoc = 0.5 * (uLoc1 * uLoc1 + uLoc2 * uLoc2 + uLoc3 * uLoc3)
      for iDoF = 1 : hFeF.DoF
        ind = hFeF.Glob[iDoF,iF]  
        hLoc += hFRef[1,iDoF,iQ] * h[ind]
      end  
      for iDoF = 1 : FeT.DoF
        GradLoc[iDoF] += sign(detDFLoc) * Weights[iQ] * fTRef[1,iDoF,iQ] * hLoc
      end  
    end
    for iQ = 1 : length(WeightsL)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,1.0,PointsL[iQ],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      u1X = 0.0
      u2X = 0.0
      for iDoF = 1 : uFeF.DoF
        ind = uFeF.Glob[iDoF,iF]  
        u1X += uFRefX[1,iDoF,iQ] * u[ind]
        u2X += uFRefX[2,iDoF,iQ] * u[ind]
      end  
      uLoc1X = 1 / detDFLoc * (DF[1,1] * u1X + DF[1,2] * u2X)
      uLoc2X = 1 / detDFLoc * (DF[2,1] * u1X + DF[2,2] * u2X)
      uLoc3X = 1 / detDFLoc * (DF[3,1] * u1X + DF[3,2] * u2X)
      hLocX = 0.5 * (uLoc1X * uLoc1X + uLoc2X * uLoc2X + uLoc3X * uLoc3X)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,PointsL[iQ],1.0,Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      u1Y = 0.0
      u2Y = 0.0
      for iDoF = 1 : uFeF.DoF
        ind = uFeF.Glob[iDoF,iF]  
        u1Y += uFRefY[1,iDoF,iQ] * u[ind]
        u2Y += uFRefY[2,iDoF,iQ] * u[ind]
      end  
      uLoc1Y = 1 / detDFLoc * (DF[1,1] * u1Y + DF[1,2] * u2Y)
      uLoc2Y = 1 / detDFLoc * (DF[2,1] * u1Y + DF[2,2] * u2Y)
      uLoc3Y = 1 / detDFLoc * (DF[3,1] * u1Y + DF[3,2] * u2Y)
      hLocY = 0.5 * (uLoc1Y * uLoc1Y + uLoc2Y * uLoc2Y + uLoc3Y * uLoc3Y)
      for iDoF = 1 : hFeF.DoF
        ind = hFeF.Glob[iDoF,iF]  
        hLocX += hFRefX[1,iDoF] * h[ind]
        hLocY += hFRefY[1,iDoF] * h[ind]
      end  
      for iDoF = 1 : FeT.DoF
        GradLoc[iDoF] += -WeightsL[iQ] * (fTRefX[1,iDoF,iQ] * hLocX + fTRefY[1,iDoF,iQ] * hLocY)
      end  
    end   
    for iDoF = 1 : FeT.DoF
      ind = FeT.Glob[iDoF,iF]  
      Rhs[ind] += GradLoc[iDoF]
    end  
  end
end

function DivRhs!(backend,FTB,Div,u,uFeF::HDivKiteDElement,h,hFeF::ScalarElement,FeT::ScalarElement,Grid,
  ElemType::Grids.ElementType,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  uFFRef  = zeros(FeT.Comp,uFeF.DoF,NumQuad)
  hFFRef  = zeros(FeT.Comp,hFeF.DoF,NumQuad)
  fTRef  = zeros(FeT.Comp,FeT.DoF,NumQuad)

  for iQ = 1 : NumQuad
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeT.DoF
        fTRef[iComp,iD,iQ] = FeT.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    for iComp = 1 : FeT.Comp
      for iD = 1 : uFeF.DoF
        uFFRef[iComp,iD,iQ] = uFeF.Divphi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    for iComp = 1 : FeT.Comp
      for iD = 1 : hFeF.DoF
        hFFRef[iComp,iD,iQ] = hFeF.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
  end
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  uFFRefX  = zeros(FeT.Comp,uFeF.DoF,NumQuadL)
  uFFRefY  = zeros(FeT.Comp,uFeF.DoF,NumQuadL)
  hFFRefX  = zeros(FeT.Comp,hFeF.DoF,NumQuadL)
  hFFRefY  = zeros(FeT.Comp,hFeF.DoF,NumQuadL)
  fTRefX  = zeros(FeT.Comp,FeT.DoF,NumQuadL)
  fTRefY  = zeros(FeT.Comp,FeT.DoF,NumQuadL)
  for iQ = 1 : NumQuadL
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeT.DoF
        fTRefX[iComp,iD,iQ] = FeT.phi[iD,iComp](-1.0,PointsL[iQ])
        fTRefY[iComp,iD,iQ] = FeT.phi[iD,iComp](PointsL[iQ],-1.0)
      end
    end
    for iD = 1 : uFeF.DoF
      uFFRefX[1,iD,iQ] = uFeF.phi[iD,1](-1.0,PointsL[iQ])
      uFFRefY[1,iD,iQ] = uFeF.phi[iD,2](PointsL[iQ],-1.0)
    end
    for iD = 1 : hFeF.DoF
      hFFRefX[1,iD,iQ] = hFeF.phi[iD,1](-1.0,PointsL[iQ])
      hFFRefY[1,iD,iQ] = hFeF.phi[iD,1](PointsL[iQ],-1.0)
    end
  end

  DivLoc = zeros(FeT.DoF)
  uFLoc = zeros(NumQuad)
  uFLocX = zeros(NumQuadL)
  uFLocY = zeros(NumQuadL)
  hFLoc = zeros(NumQuad)
  hFLocX = zeros(NumQuadL)
  hFLocY = zeros(NumQuadL)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)

  @time for iF = 1 : Grid.NumFaces
    @. uFLoc = 0 
    @. uFLocX = 0 
    @. uFLocY = 0 
    @. hFLoc = 0 
    @. hFLocX = 0 
    @. hFLocY = 0 
    @. DivLoc = 0
    @inbounds for iDoFuFeF = 1 : uFeF.DoF
      ind = uFeF.Glob[iDoFuFeF,iF]  
      @inbounds for iQ = 1 : NumQuad
        uFLoc[iQ] += uFFRef[1,iDoFuFeF,iQ] * u[ind]
      end  
      @inbounds for iQ = 1 : NumQuadL
        uFLocX[iQ] += uFFRefX[1,iDoFuFeF,iQ] * u[ind]
        uFLocY[iQ] += uFFRefY[1,iDoFuFeF,iQ] * u[ind]
      end  
    end   
    @inbounds for iDoFhFeF = 1 : hFeF.DoF
      ind = hFeF.Glob[iDoFhFeF,iF]  
      @inbounds for iQ = 1 : NumQuad
        hFLoc[iQ] += hFFRef[1,iDoFhFeF,iQ] * h[ind]
      end  
      @inbounds for iQ = 1 : NumQuadL
        hFLocX[iQ] += hFFRefX[1,iDoFhFeF,iQ] * h[ind]
        hFLocY[iQ] += hFFRefY[1,iDoFhFeF,iQ] * h[ind]
      end  
    end   
    for iQ = 1 : NumQuad
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[iQ,1],Points[iQ,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      for iDoFFeT = 1 : FeT.DoF
        DivLoc[iDoFFeT] +=  sign(detDFLoc) * Weights[iQ] * fTRef[1,iDoFFeT,iQ] * uFLoc[iQ] * hFLoc[iQ]
      end
    end
    for iQ = 1 : NumQuadL
      for iDoFFeT = 1 : FeT.DoF
        DivLoc[iDoFFeT] +=  WeightsL[iQ] * (fTRefX[1,iDoFFeT,iQ] * uFLocX[iQ] * hFLocX[iQ] +
          fTRefY[1,iDoFFeT,iQ] * uFLocY[iQ] * hFLocY[iQ])
      end
    end   
    for iDoFFeT = 1 : FeT.DoF
      ind = FeT.Glob[iDoFFeT,iF]
      Div[ind] += DivLoc[iDoFFeT]
    end
  end
end

function GradRhs!(backend,FTB,Grad,h,hFeF::ScalarElement,u,uFeF::HDivKiteDElement,FeT::HDivKiteDElement,Grid,
  ElemType::Grids.ElementType,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  uFFRef  = zeros(uFeF.Comp,uFeF.DoF,NumQuad)
  hFFRef  = zeros(hFeF.Comp,hFeF.DoF,NumQuad)
  fTRef  = zeros(hFeF.Comp,FeT.DoF,NumQuad)

  for iQ = 1 : NumQuad
    for iComp = 1 : hFeF.Comp
      for iD = 1 : FeT.DoF
        fTRef[iComp,iD,iQ] = FeT.Divphi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    for iComp = 1 : uFeF.Comp
      for iD = 1 : uFeF.DoF
        uFFRef[iComp,iD,iQ] = uFeF.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    for iComp = 1 : hFeF.Comp
      for iD = 1 : hFeF.DoF
        hFFRef[iComp,iD,iQ] = hFeF.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
  end
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  uFFRefX  = zeros(uFeF.Comp,uFeF.DoF,NumQuadL)
  uFFRefY  = zeros(uFeF.Comp,uFeF.DoF,NumQuadL)
  hFFRefX  = zeros(hFeF.Comp,hFeF.DoF,NumQuadL)
  hFFRefY  = zeros(hFeF.Comp,hFeF.DoF,NumQuadL)
  fTRefX  = zeros(hFeF.Comp,FeT.DoF,NumQuadL)
  fTRefY  = zeros(hFeF.Comp,FeT.DoF,NumQuadL)
  for iQ = 1 : NumQuadL
    for iD = 1 : FeT.DoF
      fTRefX[1,iD,iQ] = FeT.phi[iD,1](-1.0,PointsL[iQ])
      fTRefY[1,iD,iQ] = FeT.phi[iD,2](PointsL[iQ],-1.0)
    end
    for iComp = 1 : uFeF.Comp
      for iD = 1 : uFeF.DoF
        uFFRefX[iComp,iD,iQ] = uFeF.phi[iD,iComp](-1.0,PointsL[iQ])
        uFFRefY[iComp,iD,iQ] = uFeF.phi[iD,iComp](PointsL[iQ],-1.0)
      end  
    end
    for iD = 1 : hFeF.DoF
      hFFRefX[1,iD,iQ] = hFeF.phi[iD,1](-1.0,PointsL[iQ])
      hFFRefY[1,iD,iQ] = hFeF.phi[iD,1](PointsL[iQ],-1.0)
    end
  end

  GradLoc = zeros(FeT.DoF)
  uLoc = zeros(uFeF.DoF)
  uFLoc = zeros(2)
  KLoc = zeros(NumQuad)
  KLocX = zeros(NumQuadL)
  KLocY = zeros(NumQuadL)
  hFLoc = zeros(NumQuad)
  hFLocX = zeros(NumQuadL)
  hFLocY = zeros(NumQuadL)
  DF = zeros(3,2)
  detDF = zeros(1)
  detDFLoc = zeros(NumQuad)
  pinvDF = zeros(3,2)
  X = zeros(3)

  @time for iF = 1 : Grid.NumFaces
    @. GradLoc = 0
    @. hFLoc = 0
    @. hFLocX = 0
    @. hFLocY = 0
    @inbounds for iDoFuFeF = 1 : uFeF.DoF
      ind = uFeF.Glob[iDoFuFeF,iF]  
      uLoc[iDoFuFeF] = u[ind]
    end  
    @inbounds for iQ = 1 : NumQuad
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[iQ,1],Points[iQ,2],Grid.Faces[iF], Grid)
      detDFLoc[iQ] = detDF[1]
      @. uFLoc = 0
      @inbounds for iDoFuFeF = 1 : uFeF.DoF
        uFLoc[1] += uFFRef[1,iDoFuFeF,iQ] * uLoc[iDoFuFeF]
        uFLoc[2] += uFFRef[2,iDoFuFeF,iQ] * uLoc[iDoFuFeF]
      end
      u1 = (DF[1,1] * uFLoc[1] + DF[1,2] * uFLoc[2]) / detDFLoc[iQ]
      u2 = (DF[2,1] * uFLoc[1] + DF[2,2] * uFLoc[2]) / detDFLoc[iQ]
      KLoc[iQ] = 0.5 * (u1 * u1 + u2 * u2)
    end  
    @inbounds for iQ = 1 : NumQuadL
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,-1.0,Points[iQ,2],Grid.Faces[iF], Grid)
      @. uFLoc = 0
      @inbounds for iDoFuFeF = 1 : uFeF.DoF
        uFLoc[1] += uFFRefX[1,iDoFuFeF,iQ] * uLoc[iDoFuFeF]
        uFLoc[2] += uFFRefX[2,iDoFuFeF,iQ] * uLoc[iDoFuFeF]
      end
      u1 = (DF[1,1] * uFLoc[1] + DF[1,2] * uFLoc[2]) / detDF[1]
      u2 = (DF[2,1] * uFLoc[1] + DF[2,2] * uFLoc[2]) / detDF[1]
      KLocX[iQ] = 0.5 * (u1 * u1 + u2 * u2)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[iQ,1],-1.0,Grid.Faces[iF], Grid)
      @. uFLoc = 0
      @inbounds for iDoFuFeF = 1 : uFeF.DoF
        uFLoc[1] += uFFRefY[1,iDoFuFeF,iQ] * uLoc[iDoFuFeF]
        uFLoc[2] += uFFRefY[2,iDoFuFeF,iQ] * uLoc[iDoFuFeF]
      end
      u1 = (DF[1,1] * uFLoc[1] + DF[1,2] * uFLoc[2]) / detDF[1]
      u2 = (DF[2,1] * uFLoc[1] + DF[2,2] * uFLoc[2]) / detDF[1]
      KLocY[iQ] = 0.5 * (u1 * u1 + u2 * u2)
    end  
      
    @inbounds for iDoFhFeF = 1 : hFeF.DoF
      ind = hFeF.Glob[iDoFhFeF,iF]  
      @inbounds for iQ = 1 : NumQuad
        hFLoc[iQ] += hFFRef[1,iDoFhFeF,iQ] * h[ind]
      end  
      @inbounds for iQ = 1 : NumQuadL
        hFLocX[iQ] += hFFRefX[1,iDoFhFeF,iQ] * h[ind]
        hFLocY[iQ] += hFFRefY[1,iDoFhFeF,iQ] * h[ind]
      end  
    end   
    for iQ = 1 : NumQuad
      for iDoFFeT = 1 : FeT.DoF
        GradLoc[iDoFFeT] +=  sign(detDFLoc[iQ]) * Weights[iQ] * fTRef[1,iDoFFeT,iQ] * (hFLoc[iQ] + KLoc[iQ])
      end
    end
    for iQ = 1 : NumQuadL
      for iDoFFeT = 1 : FeT.DoF
        GradLoc[iDoFFeT] +=  WeightsL[iQ] * (fTRefX[1,iDoFFeT,iQ] * (hFLocX[iQ] + KLocX[iQ]) +
          fTRefY[1,iDoFFeT,iQ] * (hFLocY[iQ] + KLocY[iQ]))
      end
    end   
    for iDoFFeT = 1 : FeT.DoF
      ind = FeT.Glob[iDoFFeT,iF]
      Grad[ind] += GradLoc[iDoFFeT]
    end
  end
end

function CrossRhs!(backend,FTB,Cross,q,qFeF::ScalarElement,u,uFeF::HDivKiteDElement,FeT::HDivKiteDElement,Grid,
  ElemType::Grids.ElementType,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  qFFRef  = zeros(qFeF.Comp,qFeF.DoF,NumQuad)
  fTRef  = zeros(FeT.Comp,FeT.DoF,NumQuad)

  for iQ = 1 : NumQuad
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeT.DoF
        fTRef[iComp,iD,iQ] = FeT.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    for iComp = 1 : qFeF.Comp
      for iD = 1 : qFeF.DoF
        qFFRef[iComp,iD,iQ] = qFeF.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
  end
  uFFRef = fTRef

  CrossLoc = zeros(FeT.DoF)
  uLoc = zeros(uFeF.DoF)
  uFLoc = zeros(2)
  qFLoc = zeros(NumQuad)
  cu1 = zeros(NumQuad)
  cu2 = zeros(NumQuad)
  DF = zeros(3,2)
  detDF = zeros(1)
  detDFLoc = zeros(NumQuad)
  pinvDF = zeros(3,2)
  X = zeros(3)
  Omega = 2 * pi / 24.0 / 3600.0
  Omega = 0.0

  @time for iF = 1 : Grid.NumFaces
    @. CrossLoc = 0
    @inbounds for iDoFuFeF = 1 : uFeF.DoF
      ind = uFeF.Glob[iDoFuFeF,iF]  
      uLoc[iDoFuFeF] = u[ind]
    end  
    @inbounds for iQ = 1 : NumQuad
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[iQ,1],Points[iQ,2],Grid.Faces[iF], Grid)
      detDFLoc[iQ] = detDF[1]
      @. uFLoc = 0
      @inbounds for iDoFuFeF = 1 : uFeF.DoF
        uFLoc[1] += uFFRef[1,iDoFuFeF,iQ] * uLoc[iDoFuFeF]
        uFLoc[2] += uFFRef[2,iDoFuFeF,iQ] * uLoc[iDoFuFeF]
      end
      u1 = -(DF[2,1] * uFLoc[1] + DF[2,2] * uFLoc[2]) / detDFLoc[iQ]
      u2 = (DF[1,1] * uFLoc[1] + DF[1,2] * uFLoc[2]) / detDFLoc[iQ]
      cu1[iQ] = (DF[1,1] * u1 + DF[2,1] * u2)
      cu2[iQ] = (DF[1,2] * u1 + DF[2,2] * u2)
      sinlat = X[3] / sqrt(X[1]^2 + X[2]^2 + X[3]^2)
      qFLoc[iQ] = 2 * Omega * sinlat
    end  
      
    @inbounds for iDoFqFeF = 1 : qFeF.DoF
      ind = qFeF.Glob[iDoFqFeF,iF]  
      @inbounds for iQ = 1 : NumQuad
        qFLoc[iQ] += qFFRef[1,iDoFqFeF,iQ] * q[ind]
      end  
    end   
    @inbounds for iQ = 1 : NumQuad
      @inbounds for iDoFFeT = 1 : FeT.DoF
        CrossLoc[iDoFFeT] +=  Weights[iQ] * qFLoc[iQ] * (fTRef[1,iDoFFeT,iQ] * cu1[iQ] +
          fTRef[2,iDoFFeT,iQ] * cu2[iQ])
      end
    end
    for iDoFFeT = 1 : FeT.DoF
      ind = FeT.Glob[iDoFFeT,iF]
      Cross[ind] += CrossLoc[iDoFFeT]
    end
  end
end
