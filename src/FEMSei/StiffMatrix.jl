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


