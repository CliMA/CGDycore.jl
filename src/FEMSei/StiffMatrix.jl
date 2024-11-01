function VortCrossVel!(backend,FTB,Rhs,u,uFe::HDivConfElement,q,qFe::ScalarElement,FeT::HDivConfElement,
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
  @inbounds for iF = 1 : Grid.NumFaces
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

function DivMatrix(backend,FTB,FeF::HDivConfElement,FeT::ScalarElement,Grid,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(FeT.Type,QuadOrd)
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

  @inbounds for iF = 1 : Grid.NumFaces
    DivLoc .= 0
    for i = 1 : length(Weights)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      DivLoc += Grid.Faces[iF].Orientation * Weights[i] * (fTRef[:,:,i]' * fFRef[:,:,i])
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

function LaplMatrix(backend,FTB,FeF::ScalarElement,FeT::ScalarElement,Grid,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(Grid.Type,QuadOrd)
  fFRef  = zeros(2,FeF.DoF,length(Weights))
  fTRef  = zeros(2,FeT.DoF,length(Weights))

  for i = 1 : length(Weights)
    for iComp = 1 : 2
      for iD = 1 : FeT.DoF
        fTRef[iComp,iD,i] = FeT.Gradphi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
    for iComp = 1 : 2
      for iD = 1 : FeF.DoF
        fFRef[iComp,iD,i] = FeF.Gradphi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  RowInd = Int64[]
  ColInd = Int64[]
  Val = Float64[]
  LaplLoc = zeros(FeT.DoF,FeF.DoF)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)

  @inbounds for iF = 1 : Grid.NumFaces
    LaplLoc .= 0
    for i = 1 : length(Weights)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      LaplLoc +=  -detDFLoc * Weights[i] * ((pinvDF * fTRef[:,:,i])' * 
        (pinvDF * fFRef[:,:,i]))
    end
    for j = 1 : size(LaplLoc,2)
      for i = 1 : size(LaplLoc,1)
        push!(RowInd,FeT.Glob[i,iF])
        push!(ColInd,FeF.Glob[j,iF])
        push!(Val,LaplLoc[i,j])
      end
    end
  end
  Lapl = sparse(RowInd, ColInd, Val)
  return Lapl
end

function DivRhs!(backend,FTB,Div,u,uFeF::HDivConfElement,FeT::ScalarElement,Grid,ElemType::Grids.ElementType,
  QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  uFRef  = zeros(FeT.Comp,uFeF.DoF,length(Weights))
  fTRef  = zeros(FeT.Comp,FeT.DoF,length(Weights))

  for i = 1 : length(Weights)
    for iComp = 1 : FeT.Comp
      for iDoF = 1 : FeT.DoF
        fTRef[iComp,iDoF,i] = FeT.phi[iDoF,iComp](Points[i,1],Points[i,2])
      end
    end
    for iComp = 1 : FeT.Comp
      for iDoF = 1 : uFeF.DoF
        uFRef[iComp,iDoF,i] = uFeF.Divphi[iDoF,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  DivLoc = zeros(FeT.DoF)
  uLoc = zeros(uFeF.DoF)

  @inbounds for iF = 1 : Grid.NumFaces
    DivLoc .= 0
    for iDoF  = 1 : uFeF.DoF
      ind = uFeF.Glob[iDoF,iF]
      uLoc[iDoF] = u[ind]
    end
    for iQ = 1 : length(Weights)
      uFRefLoc = 0.0
      for iDoF = 1 : uFeF.DoF
        uFRefLoc += uFRef[1,iDoF,iQ] * uLoc[iDoF]  
      end  
      for iDoF = 1 : FeT.DoF
        DivLoc[iDoF] -= Grid.Faces[iF].Orientation * Weights[iQ] * fTRef[1,iDoF,iQ] * uFRefLoc 
      end  
    end
    for iDoF  = 1 : FeT.DoF
      ind = FeT.Glob[iDoF,iF]
      Div[ind] += DivLoc[iDoF]
    end
  end
end

function GradRhs!(backend,FTB,Grad,h,hFeF::ScalarElement,FeT::HDivConfElement,Grid,ElemType::Grids.ElementType,
  QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  fFRef  = zeros(hFeF.Comp,hFeF.DoF,length(Weights))
  fTRef  = zeros(hFeF.Comp,FeT.DoF,length(Weights))

  for i = 1 : length(Weights)
    for iComp = 1 : hFeF.Comp
      for iDoF = 1 : FeT.DoF
        fTRef[iComp,iDoF,i] = FeT.Divphi[iDoF,iComp](Points[i,1],Points[i,2])
      end
    end
    for iComp = 1 : hFeF.Comp
      for iDoF = 1 : hFeF.DoF
        fFRef[iComp,iDoF,i] = hFeF.phi[iDoF,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  GradLoc = zeros(FeT.DoF)
  hLoc = zeros(hFeF.DoF)

  @inbounds for iF = 1 : Grid.NumFaces
    GradLoc .= 0
    for iDoF  = 1 : hFeF.DoF
      ind = hFeF.Glob[iDoF,iF]
      hLoc[iDoF] = h[ind]
    end
    for iQ = 1 : length(Weights)
      fFRefLoc = 0.0
      for iDoF = 1 : hFeF.DoF
        fFRefLoc += fFRef[1,iDoF,iQ] * hLoc[iDoF]  
      end  
      for iDoF = 1 : FeT.DoF
        GradLoc[iDoF] -= Grid.Faces[iF].Orientation * Weights[iQ] * fTRef[1,iDoF,iQ] * fFRefLoc 
      end  
    end
    for iDoF  = 1 : FeT.DoF
      ind = FeT.Glob[iDoF,iF]
      Grad[ind] += GradLoc[iDoF]
    end
  end
end
#HDiv->HCurl
function CurlMatrix(backend,FTB,FeF::HCurlConfElement,FeT::ScalarElement,Grid,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(FeT.Type,QuadOrd)
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

  @inbounds for iF = 1 : Grid.NumFaces
    CurlLoc .= 0
    for i = 1 : length(Weights)
      CurlLoc += Grid.Faces[iF].Orientation * Weights[i] * (fTRef[:,:,i]' * fFRef[:,:,i])
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
  NumQuad, Weights, Points = QuadRule(FeT.Type,QuadOrd)
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

  @inbounds for iF = 1 : Grid.NumFaces
    DivLoc .= 0
    for i = 1 : length(Weights)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      DivLoc += Grid.Faces[iF].Orientation * Weights[i] * (fTRef[:,:,i]' * fFRef[:,:,i])
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

function DivRhs!(backend,FTB,Div,u,uFeF::HDivKiteDElement,FeT::ScalarElement,Grid,ElemType::Grids.ElementType,
  QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  uFRef  = zeros(FeT.Comp,uFeF.DoF,length(Weights))
  fTRef  = zeros(FeT.Comp,FeT.DoF,length(Weights))

  for i = 1 : length(Weights)
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeT.DoF
        fTRef[iComp,iD,i] = FeT.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
    for iComp = 1 : FeT.Comp
      for iD = 1 : uFeF.DoF
        uFRef[iComp,iD,i] = uFeF.Divphi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  uFRefX  = zeros(FeT.Comp,uFeF.DoF,length(WeightsL))
  uFRefY  = zeros(FeT.Comp,uFeF.DoF,length(WeightsL))
  fTRefX  = zeros(FeT.Comp,FeT.DoF,length(WeightsL))
  fTRefY  = zeros(FeT.Comp,FeT.DoF,length(WeightsL))
  for i = 1 : length(WeightsL)
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeT.DoF
        fTRefX[iComp,iD,i] = FeT.phi[iD,iComp](-1.0,PointsL[i])
        fTRefY[iComp,iD,i] = FeT.phi[iD,iComp](PointsL[i],-1.0)
      end
    end
    for iD = 1 : uFeF.DoF
      uFRefX[1,iD,i] = uFeF.phi[iD,1](-1.0,PointsL[i])
      uFRefY[1,iD,i] = uFeF.phi[iD,2](PointsL[i],-1.0)
    end
  end

  DivLoc = zeros(FeT.DoF)
  uF = zeros(uFeF.DoF)

  @inbounds for iF = 1 : Grid.NumFaces
    DivLoc .= 0
    for iDoF = 1 : uFeF.DoF
      ind = uFeF.Glob[iDoF,iF]  
      uF[iDoF] = u[ind]  
    end  
    for i = 1 : length(Weights)
      uLoc = 0.0
      for iDoF = 1 : uFeF.DoF
        uLoc += uFRef[1,iDoF,i] * uF[iDoF]  
      end  
      for iDoF = 1 : FeT.DoF
        DivLoc[iDoF] -= Grid.Faces[iF].Orientation * Weights[i] * fTRef[1,iDoF,i] * uLoc
      end  
    end
    for i = 1 : length(WeightsL)
       uLocX = 0.0 
       uLocY = 0.0 
       for iDoF = 1 : uFeF.DoF
         uLocX += uFRefX[1,iDoF,i] * uF[iDoF]  
         uLocY += uFRefY[1,iDoF,i] * uF[iDoF]  
       end  
      for iDoF = 1 : FeT.DoF
         DivLoc[iDoF] -= WeightsL[i] * (fTRefX[1,iDoF,i] * uLocX +
         fTRefY[1,iDoF,i] * uLocY)
      end   
    end   
    for iDoF = 1 : FeT.DoF
      ind = FeT.Glob[iDoF,iF]  
      Div[ind] += DivLoc[iDoF]  
    end
  end
end

function CurlMatrix(backend,FTB,FeF::HCurlKiteDElement,FeT::ScalarElement,Grid,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(FeT.Type,QuadOrd)
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

  @inbounds for iF = 1 : Grid.NumFaces
    CurlLoc .= 0
    for i = 1 : length(Weights)
      CurlLoc += Weights[i] * (fTRef[:,:,i]' * fFRef[:,:,i])
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

function GradMatrix(backend,FTB,FeF::ScalarElement,FeT::HDivConfElement,Grid,QuadOrd,Jacobi)
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

  @inbounds for iF = 1 : Grid.NumFaces
    GradLoc .= 0
    for i = 1 : length(Weights)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      GradLoc += Grid.Faces[iF].Orientation * Weights[i] * (fTRef[:,:,i]' * fFRef[:,:,i])
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

  @inbounds for iF = 1 : Grid.NumFaces
    GradLoc .= 0
    for i = 1 : length(Weights)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      GradLoc += Grid.Faces[iF].Orientation * Weights[i] * (fTRef[:,:,i]' * fFRef[:,:,i])
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

function GradKinHeight!(backend,FTB,Rhs,h,hFeF::ScalarElement,u,uFeF::HDivConfElement,
  FeT::HDivConfElement,Grid,ElemType::Grids.ElementType,QuadOrd,Jacobi)

  Grav =  9.81
  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  hFRef  = zeros(hFeF.Comp,hFeF.DoF,length(Weights))
  uFRef  = zeros(FeT.Comp,FeT.DoF,length(Weights))
  DivfTRef  = zeros(hFeF.Comp,FeT.DoF,length(Weights))


  for iQ = 1 : NumQuad
    for iComp = 1 : hFeF.Comp
      for iD = 1 : FeT.DoF
        DivfTRef[iComp,iD,iQ] = FeT.Divphi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    for iComp = 1 : uFeF.Comp
      for iD = 1 : uFeF.DoF
        uFRef[iComp,iD,iQ] = uFeF.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    for iComp = 1 : hFeF.Comp
      for iD = 1 : hFeF.DoF
        hFRef[iComp,iD,iQ] = hFeF.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
  end
  GradLoc = zeros(FeT.DoF)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  uLoc = zeros(uFeF.DoF)
  hLoc = zeros(hFeF.DoF)

  @inbounds for iF = 1 : Grid.NumFaces
    GradLoc .= 0
    for iDoF = 1 : uFeF.DoF
      ind = uFeF.Glob[iDoF,iF]  
      uLoc[iDoF] = u[ind]
    end  
    for iDoF = 1 : hFeF.DoF
      ind = hFeF.Glob[iDoF,iF]  
      hLoc[iDoF] = h[ind]
    end  
    for iQ = 1 : length(Weights)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[iQ,1],Points[iQ,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      u1 = 0.0
      u2 = 0.0
      for iDoF = 1 : uFeF.DoF
        u1 += uFRef[1,iDoF,iQ] * uLoc[iDoF]
        u2 += uFRef[2,iDoF,iQ] * uLoc[iDoF]
      end  
      uLoc1 = 1 / detDFLoc * (DF[1,1] * u1 + DF[1,2] * u2)
      uLoc2 = 1 / detDFLoc * (DF[2,1] * u1 + DF[2,2] * u2)
      uLoc3 = 1 / detDFLoc * (DF[3,1] * u1 + DF[3,2] * u2)
      hhLoc = 0.5 * (uLoc1 * uLoc1 + uLoc2 * uLoc2 + uLoc3 * uLoc3)
      for iDoF = 1 : hFeF.DoF
        hhLoc += hFRef[1,iDoF,iQ] * Grav * hLoc[iDoF]
      end  
      for iDoF = 1 : FeT.DoF
        GradLoc[iDoF] += Grid.Faces[iF].Orientation * Weights[iQ] * DivfTRef[1,iDoF,iQ] * hhLoc
      end  
    end
    for iDoF = 1 : FeT.DoF
      ind = FeT.Glob[iDoF,iF]  
      Rhs[ind] += GradLoc[iDoF]
    end  
  end
end
#! heißt drin berechnet
function GradKinHeight!(backend,FTB,Rhs,h,hFeF::ScalarElement,u,uFeF::HDivKiteDElement,
  FeT::HDivKiteDElement,Grid,ElemType::Grids.ElementType,QuadOrd,Jacobi)

  Grav =  9.81
  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  hFRef  = zeros(hFeF.Comp,hFeF.DoF,length(Weights))
  fTRef  = zeros(hFeF.Comp,FeT.DoF,length(Weights))
  uFRef  = zeros(uFeF.Comp,uFeF.DoF,length(Weights))

  for iQ = 1 : length(Weights)
    for iComp = 1 : hFeF.Comp
      for iDoF = 1 : FeT.DoF
        fTRef[iComp,iDoF,iQ] = FeT.Divphi[iDoF,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    for iComp = 1 : hFeF.Comp
      for iDoF = 1 : hFeF.DoF
        hFRef[iComp,iDoF,iQ] = hFeF.phi[iDoF,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    for iComp = 1 : uFeF.Comp
      for iDoF = 1 : uFeF.DoF
        uFRef[iComp,iDoF,iQ] = uFeF.phi[iDoF,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
  end
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  fTRefX  = zeros(hFeF.Comp,FeT.DoF,length(WeightsL))
  fTRefY  = zeros(hFeF.Comp,FeT.DoF,length(WeightsL))
  hFRefX  = zeros(hFeF.Comp,hFeF.DoF,length(WeightsL))
  hFRefY  = zeros(hFeF.Comp,hFeF.DoF,length(WeightsL))
  uFRefX  = zeros(uFeF.Comp,uFeF.DoF,length(WeightsL))
  uFRefY  = zeros(uFeF.Comp,uFeF.DoF,length(WeightsL))
  for iQ = 1 : length(WeightsL)
    for iComp = 1 : hFeF.Comp
      for iDoF = 1 : hFeF.DoF
        hFRefX[iComp,iDoF,iQ] = hFeF.phi[iDoF,iComp](-1.0,PointsL[iQ])
        hFRefY[iComp,iDoF,iQ] = hFeF.phi[iDoF,iComp](PointsL[iQ],-1.0)
      end
    end
    for iDoF = 1 : uFeF.DoF
      uFRefX[1,iDoF,iQ] = uFeF.phi[iDoF,1](-1.0,PointsL[iQ])
      uFRefX[2,iDoF,iQ] = uFeF.phi[iDoF,2](-1.0,PointsL[iQ])
      uFRefY[1,iDoF,iQ] = uFeF.phi[iDoF,1](PointsL[iQ],-1.0)
      uFRefY[2,iDoF,iQ] = uFeF.phi[iDoF,2](PointsL[iQ],-1.0)
    end
    for iDoF = 1 : FeT.DoF
      fTRefX[1,iDoF,iQ] = FeT.phi[iDoF,1](-1.0,PointsL[iQ])
      fTRefY[1,iDoF,iQ] = FeT.phi[iDoF,2](PointsL[iQ],-1.0)
    end
  end
  GradLoc = zeros(FeT.DoF)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  uF = zeros(uFeF.DoF)
  hF = zeros(hFeF.DoF)

  @inbounds for iF = 1 : Grid.NumFaces
    GradLoc .= 0
    for iDoF = 1 : uFeF.DoF
      ind = uFeF.Glob[iDoF,iF]  
      uF[iDoF] = u[ind]
    end  
    for iDoF = 1 : hFeF.DoF
      ind = hFeF.Glob[iDoF,iF]  
      hF[iDoF] = h[ind]
    end  
    for iQ = 1 : length(Weights)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[iQ,1],Points[iQ,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      u1 = 0.0
      u2 = 0.0
      for iDoF = 1 : uFeF.DoF
        u1 += uFRef[1,iDoF,iQ] * uF[iDoF]
        u2 += uFRef[2,iDoF,iQ] * uF[iDoF]
      end
      uLoc1 = 1 / detDFLoc * (DF[1,1] * u1 + DF[1,2] * u2)
      uLoc2 = 1 / detDFLoc * (DF[2,1] * u1 + DF[2,2] * u2)
      uLoc3 = 1 / detDFLoc * (DF[3,1] * u1 + DF[3,2] * u2)
      hLoc = 0.5 * (uLoc1 * uLoc1 + uLoc2 * uLoc2 + uLoc3 * uLoc3)
      for iDoF = 1 : hFeF.DoF 
        hLoc += hFRef[1,iDoF,iQ] * Grav * hF[iDoF]
      end  
      for iDoF = 1 : FeT.DoF
        GradLoc[iDoF] += Grid.Faces[iF].Orientation * Weights[iQ] * fTRef[1,iDoF,iQ] * hLoc
      end
    end
    for iQ = 1 : length(WeightsL)
      hLocX = 0.0 
      hLocY = 0.0 
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,-1.0,PointsL[iQ],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      u1X = 0.0
      u2X = 0.0
      for iDoF = 1 : uFeF.DoF
        u1X += uFRefX[1,iDoF,iQ] * uF[iDoF]
        u2X += uFRefX[2,iDoF,iQ] * uF[iDoF]
      end
      uLoc1X = 1 / detDFLoc * (DF[1,1] * u1X + DF[1,2] * u2X)
      uLoc2X = 1 / detDFLoc * (DF[2,1] * u1X + DF[2,2] * u2X)
      uLoc3X = 1 / detDFLoc * (DF[3,1] * u1X + DF[3,2] * u2X)
      hLocX = 0.5 * (uLoc1X * uLoc1X + uLoc2X * uLoc2X + uLoc3X * uLoc3X)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,PointsL[iQ],-1.0,Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      u1Y = 0.0
      u2Y = 0.0
      for iDoF = 1 : uFeF.DoF
        u1Y += uFRefY[1,iDoF,iQ] * uF[iDoF]
        u2Y += uFRefY[2,iDoF,iQ] * uF[iDoF]
      end
      uLoc1Y = 1 / detDFLoc * (DF[1,1] * u1Y + DF[1,2] * u2Y)
      uLoc2Y = 1 / detDFLoc * (DF[2,1] * u1Y + DF[2,2] * u2Y)
      uLoc3Y = 1 / detDFLoc * (DF[3,1] * u1Y + DF[3,2] * u2Y)
      hLocY = 0.5 * (uLoc1Y * uLoc1Y + uLoc2Y * uLoc2Y + uLoc3Y * uLoc3Y)
      for iDoF = 1 : hFeF.DoF
        hLocX += hFRefX[1,iDoF,iQ] * Grav * hF[iDoF]
        hLocY += hFRefY[1,iDoF,iQ] * Grav * hF[iDoF]
      end
      for iDoF = 1 : FeT.DoF
        GradLoc[iDoF] += WeightsL[iQ] * (fTRefX[1,iDoF,iQ] * hLocX + fTRefY[1,iDoF,iQ] * hLocY)
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
        fTRefX[iComp,iD,iQ] = FeT.phi[iD,iComp](1.0,PointsL[iQ])
        fTRefY[iComp,iD,iQ] = FeT.phi[iD,iComp](PointsL[iQ],1.0)
      end
    end
    for iD = 1 : uFeF.DoF
      uFFRefX[1,iD,iQ] = uFeF.phi[iD,1](1.0,PointsL[iQ])
      uFFRefY[1,iD,iQ] = uFeF.phi[iD,2](PointsL[iQ],1.0)
    end
    for iD = 1 : hFeF.DoF
      hFFRefX[1,iD,iQ] = hFeF.phi[iD,1](1.0,PointsL[iQ])
      hFFRefY[1,iD,iQ] = hFeF.phi[iD,1](PointsL[iQ],1.0)
    end
  end

  DivLoc = zeros(FeT.DoF)
  uFLoc = zeros(NumQuad)
  uFLocX = zeros(NumQuadL)
  uFLocY = zeros(NumQuadL)
  hFLoc = zeros(NumQuad)
  hFLocX = zeros(NumQuadL)
  hFLocY = zeros(NumQuadL)

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
      for iDoFFeT = 1 : FeT.DoF
        DivLoc[iDoFFeT] +=  Weights[iQ] * fTRef[1,iDoFFeT,iQ] * uFLoc[iQ] * hFLoc[iQ]
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

  @inbounds for iF = 1 : Grid.NumFaces
    @. GradLoc = 0
    @. hFLoc = 0
    @. hFLocX = 0
    @. hFLocY = 0
    for iDoFuFeF = 1 : uFeF.DoF
      ind = uFeF.Glob[iDoFuFeF,iF]  
      uLoc[iDoFuFeF] = u[ind]
    end  
    for iQ = 1 : NumQuad
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[iQ,1],Points[iQ,2],Grid.Faces[iF], Grid)
      detDFLoc[iQ] = detDF[1]
      @. uFLoc = 0
      for iDoFuFeF = 1 : uFeF.DoF
        uFLoc[1] += uFFRef[1,iDoFuFeF,iQ] * uLoc[iDoFuFeF]
        uFLoc[2] += uFFRef[2,iDoFuFeF,iQ] * uLoc[iDoFuFeF]
      end
      u1 = (DF[1,1] * uFLoc[1] + DF[1,2] * uFLoc[2]) / detDFLoc[iQ]
      u2 = (DF[2,1] * uFLoc[1] + DF[2,2] * uFLoc[2]) / detDFLoc[iQ]
      KLoc[iQ] = 0.5 * (u1 * u1 + u2 * u2)
    end  
    for iQ = 1 : NumQuadL
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,-1.0,Points[iQ,2],Grid.Faces[iF], Grid)
      @. uFLoc = 0
      for iDoFuFeF = 1 : uFeF.DoF
        uFLoc[1] += uFFRefX[1,iDoFuFeF,iQ] * uLoc[iDoFuFeF]
        uFLoc[2] += uFFRefX[2,iDoFuFeF,iQ] * uLoc[iDoFuFeF]
      end
      u1 = (DF[1,1] * uFLoc[1] + DF[1,2] * uFLoc[2]) / detDF[1]
      u2 = (DF[2,1] * uFLoc[1] + DF[2,2] * uFLoc[2]) / detDF[1]
      KLocX[iQ] = 0.5 * (u1 * u1 + u2 * u2)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[iQ,1],-1.0,Grid.Faces[iF], Grid)
      @. uFLoc = 0
      for iDoFuFeF = 1 : uFeF.DoF
        uFLoc[1] += uFFRefY[1,iDoFuFeF,iQ] * uLoc[iDoFuFeF]
        uFLoc[2] += uFFRefY[2,iDoFuFeF,iQ] * uLoc[iDoFuFeF]
      end
      u1 = (DF[1,1] * uFLoc[1] + DF[1,2] * uFLoc[2]) / detDF[1]
      u2 = (DF[2,1] * uFLoc[1] + DF[2,2] * uFLoc[2]) / detDF[1]
      KLocY[iQ] = 0.5 * (u1 * u1 + u2 * u2)
    end  
      
    for iDoFhFeF = 1 : hFeF.DoF
      ind = hFeF.Glob[iDoFhFeF,iF]  
      for iQ = 1 : NumQuad
        hFLoc[iQ] += hFFRef[1,iDoFhFeF,iQ] * h[ind]
      end  
      for iQ = 1 : NumQuadL
        hFLocX[iQ] += hFFRefX[1,iDoFhFeF,iQ] * h[ind]
        hFLocY[iQ] += hFFRefY[1,iDoFhFeF,iQ] * h[ind]
      end  
    end   
    for iQ = 1 : NumQuad
      for iDoFFeT = 1 : FeT.DoF
        GradLoc[iDoFFeT] +=  Grid.Faces[iF].Orientation * Weights[iQ] * fTRef[1,iDoFFeT,iQ] * (hFLoc[iQ] + KLoc[iQ])
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

function CrossRhs!(backend,FTB,Cross,q,qFeF::ScalarElement,u,uFeF::HDivElement,FeT::HDivElement,Grid,
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
  qLoc = zeros(qFeF.DoF)
  cu1 = zeros(NumQuad)
  cu2 = zeros(NumQuad)
  DF = zeros(3,2)
  detDF = zeros(1)
  detDFLoc = zeros(NumQuad)
  pinvDF = zeros(3,2)
  X = zeros(3)
  Omega = 2 * pi / 24.0 / 3600.0

  @inbounds for iF = 1 : Grid.NumFaces
    @. CrossLoc = 0
    for iDoF = 1 : uFeF.DoF
      ind = uFeF.Glob[iDoF,iF]  
      uLoc[iDoF] = u[ind]
    end  
    for iDoF = 1 : qFeF.DoF
      ind = qFeF.Glob[iDoF,iF]  
      qLoc[iDoF] = q[ind]
    end  
    for iQ = 1 : NumQuad
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[iQ,1],Points[iQ,2],Grid.Faces[iF], Grid)
      detDFLoc[iQ] = detDF[1]
      uFLoc1 = 0.0
      uFLoc2 = 0.0
      for iDoF = 1 : uFeF.DoF
        uFLoc1 += uFFRef[1,iDoF,iQ] * uLoc[iDoF]
        uFLoc2 += uFFRef[2,iDoF,iQ] * uLoc[iDoF]
      end
      u1 = (DF[1,1] * uFLoc1 + DF[1,2] * uFLoc2) 
      u2 = (DF[2,1] * uFLoc1 + DF[2,2] * uFLoc2)
      u3 = (DF[3,1] * uFLoc1 + DF[3,2] * uFLoc2)
      Rad = sqrt(X[1]^2 + X[2]^2 + X[3]^2)
      k1 = X[1] / Rad
      k2 = X[2] / Rad
      k3 = X[3] / Rad
      kcru1 = k2 * u3 - k3 * u2
      kcru2 = -(k1 * u3 - k3 * u1)
      kcru3 = k1 * u2 - k2 * u1
      cu1 = (DF[1,1] * kcru1 + DF[2,1] * kcru2 + DF[3,1] * kcru3)
      cu2 = (DF[1,2] * kcru1 + DF[2,2] * kcru2 + DF[3,2] * kcru3)
      sinlat = k3
      qFLoc = 2 * Omega * sinlat
      for iDoF = 1 : qFeF.DoF
        qFLoc += qFFRef[1,iDoF,iQ] * qLoc[iDoF]
      end
      for iDoF = 1 : FeT.DoF
        CrossLoc[iDoF] -=  Weights[iQ] * qFLoc * (fTRef[1,iDoF,iQ] * cu1 +
          fTRef[2,iDoF,iQ] * cu2) / detDFLoc[iQ]
      end    
    end  
    for iDoF = 1 : FeT.DoF
      ind = FeT.Glob[iDoF,iF]
      Cross[ind] += CrossLoc[iDoF]
    end
  end
end

#vector field
function DivMomentumVector!(backend,FTB,Rhs,uHDiv,FeHDiv::HDivElement,uVecDG,FeVecDG::VectorElement,FeTHDiv::HDivElement,Grid,
  ElemType::Grids.ElementType,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  fuHDiv  = zeros(FeHDiv.DoF,FeTHDiv.Comp,NumQuad)
  fuVecDG  = zeros(FeVecDG.DoF,3,NumQuad)
  fuTHDiv  = zeros(FeTHDiv.DoF,FeTHDiv.Comp,NumQuad)
  fuDivHDiv = zeros(FeHDiv.DoF,NumQuad)
  fuGradVecDG = zeros(FeVecDG.DoF,3,2,NumQuad)

  #computation of the ansatz functions in the quadrature points
  for iQ = 1 : NumQuad
    for iComp = 1 : FeTHDiv.Comp
      for iD = 1 : FeTHDiv.DoF
        fuTHDiv[iD,iComp,iQ] = FeTHDiv.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    for iComp = 1 : FeTHDiv.Comp
      for iD = 1 : FeHDiv.DoF
        fuHDiv[iD,iComp,iQ] = FeHDiv.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    for iD = 1 : FeVecDG.DoF
      fuVecDG[iD,1,iQ] = FeVecDG.phi[iD,1](Points[iQ,1],Points[iQ,2])
      fuVecDG[iD,2,iQ] = FeVecDG.phi[iD,2](Points[iQ,1],Points[iQ,2])
      fuVecDG[iD,3,iQ] = FeVecDG.phi[iD,3](Points[iQ,1],Points[iQ,2])
    end
    for iD = 1 : FeHDiv.DoF
      fuDivHDiv[iD,iQ] = FeHDiv.Divphi[iD,1](Points[iQ,1],Points[iQ,2])
    end
    
    for iD = 1 : FeVecDG.DoF
      fuGradVecDG[iD,1,1,iQ] = FeVecDG.Gradphi[iD,1,1](Points[iQ,1],Points[iQ,2])
      fuGradVecDG[iD,1,2,iQ] = FeVecDG.Gradphi[iD,1,2](Points[iQ,1],Points[iQ,2])
      fuGradVecDG[iD,2,1,iQ] = FeVecDG.Gradphi[iD,2,1](Points[iQ,1],Points[iQ,2])
      fuGradVecDG[iD,2,2,iQ] = FeVecDG.Gradphi[iD,2,2](Points[iQ,1],Points[iQ,2])
      fuGradVecDG[iD,3,1,iQ] = FeVecDG.Gradphi[iD,3,1](Points[iQ,1],Points[iQ,2])
      fuGradVecDG[iD,3,2,iQ] = FeVecDG.Gradphi[iD,3,2](Points[iQ,1],Points[iQ,2])
    end
  end


  #local variables
  uVecDGLoc = zeros(FeVecDG.DoF) #c
  uHDivLoc = zeros(FeHDiv.DoF)
  RhsLoc = zeros(FeTHDiv.DoF)

  uHDivLocLeft = zeros(FeHDiv.DoF)
  uHDivLocRight = zeros(FeHDiv.DoF)
  uVecDGLocLeft = zeros(FeDG.DoF)
  uVecDGLocRight = zeros(FeDG.DoF)

  uuVecDGLoc = zeros(3,1)
  uuGradVecDGLoc = zeros(3,2)
  uuHDivLoc = zeros(2,1)
  uuDivHDivLoc = 0.0


  DFL = zeros(3,2)
  detDFL = zeros(1)
  pinvDFL = zeros(3,2)
  XL = zeros(3)

  DFR = zeros(3,2)
  detDFR = zeros(1)
  pinvDFR = zeros(3,2)
  XR = zeros(3)

  #copy from global variables to local variables  
  for iF = 1 : Grid.NumFaces
    for iD = 1 : FeVecDG.DoF
      ind = FeVecDG.Glob[iD,iF]  
      uVecDGLoc[iD] = uVecDG[ind]
    end  
    for iD = 1 : FeHDiv.DoF
      ind = FeHDiv.Glob[iD,iF]  
      uHDivLoc[iD] = uHDiv[ind]
    end  
    @. RhsLoc = 0.0

    for iQ = 1 : NumQuad
      uuVecDGLoc = 0.0
      @. uuGradVecDGLoc = 0.0
      #computation of local variables in a quadrature point
      for iD = 1 : FeVecDG.DoF
          uuVecDGLoc += uVecDGLoc[iD] * fuVecDG[iD,1,iQ] 
          @views @. uuGradVecDGLoc[:,:] += fuGradVecDG[iD,1,:,iQ] 
      end
      uuDivHDivLoc = 0.0
      @. uuHDivLoc = 0.0
      for iD = 1 : FeHDiv.DoF
          uuDivHDivLoc += uHDivLoc[iD] * fuDivHDiv[iD,iQ] 
          uuHDivLoc[1] += uHDivLoc[iD] * fuHDiv[iD,1,iQ]
          uuHDivLoc[2] += uHDivLoc[iD] * fuHDiv[iD,2,iQ]  
      end
      GradDGHDiv = uuGradDGLoc' * uuHDivLoc
      DGDivHDiv = uuDivHDivLoc * uuDGLoc
      innersum = Grid.Faces[iF].Orientation * (DGDivHDiv + GradDGHDiv)
      #product incoming functions and test function
      for iD = 1 : FeT.DoF
        product = innersum * fuT[iD,1,iQ]
        RhsLoc[iD] += Weights[iQ] * product
      end 
    end

    #save the new without upwind 
    for iD = 1 : FeT.DoF
      ind = FeT.Glob[iD,iF] 
      Rhs[ind] += RhsLoc[iD]
    end  
  end 

  #UPWIND on Edges
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  uFFRef  = zeros(FeHDiv.DoF, FeHDiv.Comp,NumQuadL,4)
  hFFRef  = zeros(FeVecDG.DoF,FeVecDG.DoF,NumQuadL,4)
  fTRef  = zeros(FeTHDiv.DoF,FeTHDiv.Comp,NumQuadL,4)
  uFFLocL = zeros(FeHDiv.DoF)
  uFFLocR = zeros(FeHDiv.DoF)
  
  for iQ = 1 : NumQuadL 
    for iComp = 1 : FeTHDiv.Comp
      for iD = 1 : FeTHDiv.DoF
        fTRef[iD,iComp,iQ,1] = FeTHDiv.phi[iD,iComp](PointsL[iQ],-1.0)
        fTRef[iD,iComp,iQ,2] = FeTHDiv.phi[iD,iComp](1.0,PointsL[iQ])
        fTRef[iD,iComp,iQ,3] = FeTHDiv.phi[iD,iComp](PointsL[iQ],1.0)
        fTRef[iD,iComp,iQ,4] = FeTHDiv.phi[iD,iComp](-1.0,PointsL[iQ])
      end
    end
    #uf RT UHDiv
    for iComp = 1 : FeHDiv.Comp
      for iD = 1 : FeHDiv.DoF
        uFFRef[iD,iComp,iQ,1] = FeHDiv.phi[iD,iComp](PointsL[iQ],-1.0)
        uFFRef[iD,iComp,iQ,2] = FeHDiv.phi[iD,iComp](1.0,PointsL[iQ])
        uFFRef[iD,iComp,iQ,3] = FeHDiv.phi[iD,iComp](PointsL[iQ],1.0)
        uFFRef[iD,iComp,iQ,4] = FeHDiv.phi[iD,iComp](-1.0,PointsL[iQ])
      end
    end
    #hf
    for iComp = 1 : FeVecDG.Comp
      for iD = 1 : FeVecDG.DoF
        hFFRef[iD,iComp,iQ,1] = FeVecDG.phi[iD,iComp](PointsL[iQ],-1.0)
        hFFRef[iD,iComp,iQ,2] = FeVecDG.phi[iD,iComp](1.0,PointsL[iQ])
        hFFRef[iD,iComp,iQ,3] = FeVecDG.phi[iD,iComp](PointsL[iQ],1.0)
        hFFRef[iD,iComp,iQ,4] = FeVecDG.phi[iD,iComp](-1.0,PointsL[iQ])
      end
    end
  end
  for iE in 1:Grid.NumEdges
    Edge = Grid.Edges[iE]
    if length(Edge.F) > 1
      iFL = Edge.F[1]
      EdgeTypeL = Edge.FE[1]
      iFR = Edge.F[2]
      EdgeTypeR = Edge.FE[2]
      #computation normales of edges
      nBarLocL = Grid.nBar[:, EdgeTypeL]
      nBarLocR = Grid.nBar[:, EdgeTypeR]
      #gamma upwind value
      uE = 0.0 
      for iD = 1 : FeHDiv.DoF
        ind = FeHDiv.Glob[iD,iFL]  
        uHDivLocLeft[iD] = uHDiv[ind]
        ind = FeHDiv.Glob[iD,iFR]  
        uHDivLocRight[iD] = uHDiv[ind]
      end 
      for iD = 1 : FeVecDG.DoF
        ind = FeVecDG.Glob[iD,iFL]  
        uVecDGLocLeft[iD] = uVecDG[ind]
        ind = FeVecDG.Glob[iD,iFR]  
        uVecDGLocRight[iD] = uVecDG[ind]
      end 
      for iQ = 1:length(WeightsL)
        for iD = 1 : FeHDiv.DoF
          uE += WeightsL[iQ]*(uHDivLocLeft[iD]*(uFFRef[iD,:,iQ,EdgeTypeL]'*nBarLocL)+uHDivLocRight[iD]*(uFFRef[iD,:,iQ,EdgeTypeR]'*nBarLocR))
        end
      end
      #upwind value
      gammaU = 0.5
      gammaLoc = uE > 0 ? gammaU : -gammaU

      #allocation Mloc
      MLoc11 = zeros(FeTHDiv.DoF, FeHDiv.DoF)
      MLoc12 = zeros(FeTHDiv.DoF, FeHDiv.DoF)
      MLoc21 = zeros(FeTHDiv.DoF, FeHDiv.DoF)
      MLoc22 = zeros(FeTHDiv.DoF, FeHDiv.DoF)

      for iQ in 1:length(WeightsL)
        #computation of Jacobi EdgetypeL
        Jacobi!(DFL,detDFL,pinvDFL,XL,Grid.Type,Points[iQ,1],Points[iQ,2],Grid.Faces[iFL], Grid)
        detDFLLoc = detDFL[1]
        #EdgeTypeR
        Jacobi!(DFR,detDFR,pinvDFR,XR,Grid.Type,Points[iQ,1],Points[iQ,2],Grid.Faces[iFR], Grid)
        detDFLRoc = detDFR[1]

        fTLocL = (1/detDFL) * DFL * fTRef[:,:,iQ,EdgeTypeL]
        fTLocR = (1/detDFR) * DFR * fTRef[:,:,iQ,EdgeTypeR]
        for iD in 1:FeHDiv.DoF
          uFFLocL[iD] = uFFRef[iD,:,iQ,EdgeTypeL]' * nBarLocL
          uFFLocR[iD] = uFFRef[iD,:,iQ,EdgeTypeR]' * nBarLocR
        end
        hFFL = 0.0
        hFFR = 0.0
        for iD in 1:FeVecDG.DoF
          hFFL += (1/detDFL) * DFL * hFFRef[iD,1,iQ,EdgeTypeL] * cDGLocLeft[iD]
          hFFR += (1/detDFR) * DFR * hFFRef[iD,1,iQ,EdgeTypeR] * cDGLocRight[iD]
        end
        MLoc11 += Weights[iQ] * (fTRef[:,1,iQ,EdgeTypeL] * hFFL) * uFFLocL'
        MLoc12 += Weights[iQ] * (fTRef[:,1,iQ,EdgeTypeR] * hFFR) * uFFLocR'
        MLoc21 += Weights[iQ] * (fTRef[:,1,iQ,EdgeTypeL] * hFFL) * uFFLocL'
        MLoc22 += Weights[iQ] * (fTRef[:,1,iQ,EdgeTypeR] * hFFR) * uFFLocR'
      end
      
      #all together over edges
      for iD = 1 : FeTHDiv.DoF
        indL = FeTHDiv.Glob[iD,iFL] 
        Rhs[indL] += (+0.5 - gammaLoc) * MLoc11[iD,:]' * uHDivLocLeft
        Rhs[indL] += (-0.5 + gammaLoc) * MLoc12[iD,:]' * uHDivLocRight
        indR = FeTHDiv.Glob[iD,iFR] 
        Rhs[indR] += (+0.5 + gammaLoc) * MLoc21[iD,:]' * uHDivLocLeft
        Rhs[indR] += (-0.5 - gammaLoc) * MLoc22[iD,:]' * uHDivLocRight
      end
    end
  end
end

#scalar variant
function DivMomentumScalar!(backend,FTB,Rhs,uHDiv,FeHDiv::HDivElement,cDG,FeDG::ScalarElement,FeT::ScalarElement,Grid,
  ElemType::Grids.ElementType,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  fuHDiv  = zeros(FeHDiv.DoF,FeHDiv.Comp,NumQuad) #4x1
  fuDG  = zeros(FeDG.DoF,FeDG.DoF,NumQuad) #1x1
  fuT  = zeros(FeT.DoF,FeT.Comp,NumQuad) #1x1
  fuDivHDiv = zeros(FeHDiv.DoF,NumQuad) 
  fuGradDG = zeros(FeDG.DoF,FeDG.DoF,2,NumQuad) #1x2

  #computation of the ansatz functions in the quadrature points
  for iQ = 1 : NumQuad
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeT.DoF
        fuT[iD,iComp,iQ] = FeT.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeHDiv.DoF
        fuHDiv[iD,iComp,iQ] = FeHDiv.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    for iD = 1 : FeDG.DoF
      fuGradDG[iD,1,1,iQ] = FeDG.Gradphi[iD,1,1](Points[iQ,1],Points[iQ,2])
      fuGradDG[iD,1,2,iQ] = FeDG.Gradphi[iD,1,2](Points[iQ,1],Points[iQ,2])
    end
    for iD = 1 : FeDG.DoF
      fuDG[iD,1,iQ] = FeDG.phi[iD,1](Points[iQ,1],Points[iQ,2])
    end
    for iD = 1 : FeHDiv.DoF
      fuDivHDiv[iD,iQ] = FeHDiv.Divphi[iD,1](Points[iQ,1],Points[iQ,2])
    end
  end

  #local variables
  cDGLoc = zeros(FeDG.DoF)
  uHDivLoc = zeros(FeHDiv.DoF)
  RhsLoc = zeros(FeT.DoF)
  uHDivLocLeft = zeros(FeHDiv.DoF)
  uHDivLocRight = zeros(FeHDiv.DoF)
  cDGLocLeft = zeros(FeDG.DoF)
  cDGLocRight = zeros(FeDG.DoF)

  uuGradDGLoc = zeros(2)
  uuHDivLoc = zeros(2)
 
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
 
  #copy from global variables to local variables  
  for iF = 1 : Grid.NumFaces
    for iD = 1 : FeDG.DoF
      ind = FeDG.Glob[iD,iF]  
      cDGLoc[iD] = cDG[ind]
    end 
    for iD = 1 : FeHDiv.DoF
      ind = FeHDiv.Glob[iD,iF]  
      uHDivLoc[iD] = uHDiv[ind]
    end  
    @. RhsLoc = 0.0

    for iQ = 1 : NumQuad
      #computation of local variables in a quadrature point
      uuDGLoc = 0.0
      @. uuGradDGLoc = 0.0
      for iD = 1 : FeDG.DoF
        uuDGLoc += cDGLoc[iD] * fuDG[iD,1,iQ] 
        @views @. uuGradDGLoc += cDGLoc[iD] * fuGradDG[iD,1,:,iQ] 
      end
      uuDivHDivLoc = 0.0
      @. uuHDivLoc = 0.0
      for iD = 1 : FeHDiv.DoF
        uuDivHDivLoc += uHDivLoc[iD] * fuDivHDiv[iD,iQ] 
        uuHDivLoc[1] += uHDivLoc[iD] * fuHDiv[iD,1,iQ]
        uuHDivLoc[2] += uHDivLoc[iD] * fuHDiv[iD,2,iQ] 
      end
      GradDGHDiv = uuGradDGLoc' * uuHDivLoc
      DGDivHDiv = uuDivHDivLoc * uuDGLoc
      innersum = Grid.Faces[iF].Orientation * (DGDivHDiv + GradDGHDiv)
      #product incoming functions and test function
      for iD = 1 : FeT.DoF
        product = innersum * fuT[iD,1,iQ]
        RhsLoc[iD] += Weights[iQ] * product
      end 
    end

    #save the new without upwind 
    for iD = 1 : FeT.DoF
      ind = FeT.Glob[iD,iF] 
      Rhs[ind] += RhsLoc[iD]
    end  
  end 

  #UPWIND on Edges
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  uFFRef  = zeros(FeHDiv.DoF, FeHDiv.Comp,NumQuadL,4)
  hFFRef  = zeros(FeDG.DoF,FeDG.DoF,NumQuadL,4)
  fTRef  = zeros(FeT.DoF,FeT.Comp,NumQuadL,4)
  uFFLocL = zeros(FeHDiv.DoF)
  uFFLocR = zeros(FeHDiv.DoF)

  #testfunctions DG scalar 
  for iQ = 1 : NumQuadL 
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeT.DoF
        fTRef[iD,iComp,iQ,1] = FeT.phi[iD,iComp](PointsL[iQ],-1.0)
        fTRef[iD,iComp,iQ,2] = FeT.phi[iD,iComp](1.0,PointsL[iQ])
        fTRef[iD,iComp,iQ,3] = FeT.phi[iD,iComp](PointsL[iQ],1.0)
        fTRef[iD,iComp,iQ,4] = FeT.phi[iD,iComp](-1.0,PointsL[iQ])
      end
    end
    #uf RT UHDiv
    for iComp = 1 : FeHDiv.Comp
      for iD = 1 : FeHDiv.DoF
        uFFRef[iD,iComp,iQ,1] = FeHDiv.phi[iD,iComp](PointsL[iQ],-1.0)
        uFFRef[iD,iComp,iQ,2] = FeHDiv.phi[iD,iComp](1.0,PointsL[iQ])
        uFFRef[iD,iComp,iQ,3] = FeHDiv.phi[iD,iComp](PointsL[iQ],1.0)
        uFFRef[iD,iComp,iQ,4] = FeHDiv.phi[iD,iComp](-1.0,PointsL[iQ])
      end
    end
    #hf
    for iComp = 1 : FeDG.Comp
      for iD = 1 : FeDG.DoF
        hFFRef[iD,iComp,iQ,1] = FeDG.phi[iD,iComp](PointsL[iQ],-1.0)
        hFFRef[iD,iComp,iQ,2] = FeDG.phi[iD,iComp](1.0,PointsL[iQ])
        hFFRef[iD,iComp,iQ,3] = FeDG.phi[iD,iComp](PointsL[iQ],1.0)
        hFFRef[iD,iComp,iQ,4] = FeDG.phi[iD,iComp](-1.0,PointsL[iQ])
      end
    end
  end
  for iE in 1:Grid.NumEdges
    Edge = Grid.Edges[iE]
    if length(Edge.F) > 1
      iFL = Edge.F[1]
      EdgeTypeL = Edge.FE[1]
      iFR = Edge.F[2]
      EdgeTypeR = Edge.FE[2]
      #computation normales of edges
      nBarLocL = Grid.nBar[:, EdgeTypeL]
      nBarLocR = Grid.nBar[:, EdgeTypeR]
      #gamma upwind value
      uE = 0.0 
      for iD = 1 : FeHDiv.DoF
        ind = FeHDiv.Glob[iD,iFL]  
        uHDivLocLeft[iD] = uHDiv[ind]
        ind = FeHDiv.Glob[iD,iFR]  
        uHDivLocRight[iD] = uHDiv[ind]
      end 
      for iD = 1 : FeDG.DoF
        ind = FeDG.Glob[iD,iFL]  
        cDGLocLeft[iD] = cDG[ind]
        ind = FeDG.Glob[iD,iFR]  
        cDGLocRight[iD] = cDG[ind]
      end 
      for iQ = 1:length(WeightsL)
        for iD = 1 : FeHDiv.DoF
          uE += WeightsL[iQ]*(uHDivLocLeft[iD]*(uFFRef[iD,:,iQ,EdgeTypeL]'*nBarLocL))
        end
      end
      #upwind value
      gammaU = 0.5
      gammaLoc = uE > 0 ? gammaU : -gammaU

      #allocation Mloc
      MLoc11 = zeros(FeT.DoF, FeHDiv.DoF)
      MLoc12 = zeros(FeT.DoF, FeHDiv.DoF)
      MLoc21 = zeros(FeT.DoF, FeHDiv.DoF)
      MLoc22 = zeros(FeT.DoF, FeHDiv.DoF)

      for iQ in 1:length(WeightsL)
        fTLocL = fTRef[:,:,iQ,EdgeTypeL]
        fTLocR = fTRef[:,:,iQ,EdgeTypeR]
        for iD in 1:FeHDiv.DoF
          uFFLocL[iD] = uFFRef[iD,:,iQ,EdgeTypeL]' * nBarLocL
          uFFLocR[iD] = uFFRef[iD,:,iQ,EdgeTypeR]' * nBarLocR
        end
        hFFL = 0.0
        hFFR = 0.0
        for iD in 1:FeDG.DoF
          hFFL += hFFRef[iD,1,iQ,EdgeTypeL] * cDGLocLeft[iD]
          hFFR += hFFRef[iD,1,iQ,EdgeTypeR] * cDGLocRight[iD]
        end
        MLoc11 += Weights[iQ] * (fTRef[:,1,iQ,EdgeTypeL] * hFFL) * uFFLocL'
        MLoc12 += Weights[iQ] * (fTRef[:,1,iQ,EdgeTypeR] * hFFR) * uFFLocR'
        MLoc21 += Weights[iQ] * (fTRef[:,1,iQ,EdgeTypeL] * hFFL) * uFFLocL'
        MLoc22 += Weights[iQ] * (fTRef[:,1,iQ,EdgeTypeR] * hFFR) * uFFLocR'
      end
      
      #all together over edges
      for iD = 1 : FeT.DoF
        indL = FeT.Glob[iD,iFL] 
        Rhs[indL] += (+0.5 - gammaLoc) * MLoc11[iD,:]' * uHDivLocLeft
        Rhs[indL] += (-0.5 + gammaLoc) * MLoc12[iD,:]' * uHDivLocRight
        indR = FeT.Glob[iD,iFR] 
        Rhs[indR] += (+0.5 + gammaLoc) * MLoc21[iD,:]' * uHDivLocLeft
        Rhs[indR] += (-0.5 - gammaLoc) * MLoc22[iD,:]' * uHDivLocRight
      end
    end
  end
end
