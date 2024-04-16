function VortCrossVel!(backend,FTB,Rhs,u,uFe::HDivElement,q,qFe::ScalarElement,FeT::HDivElement,
  Grid,QuadOrd,Jacobi)
  QQ = FEMSei.QuadRule{FTB}(Grid.Type,backend,QuadOrd)
  Weights = QQ.Weights
  Points = QQ.Points
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
  for iF = 1 : Grid.NumFaces
    RhsLoc .= 0
    for i = 1 : length(Weights)
      DF, detJ = Jacobi(Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      qLoc = qfFRef[1,1,i] * q[qFe.Glob[1,iF]]
      for j = 2 : qFe.DoF 
        qLoc += qfFRef[1,j,i] * q[qFe.Glob[j,iF]]
      end  
      uLoc .= [-ufFRef[2,1,i],ufFRef[1,1,i]] * u[uFe.Glob[1,iF]]
      for j = 2 : FeT.DoF 
        uLoc .+= [-ufFRef[2,j,i],ufFRef[1,j,i]] * u[uFe.Glob[j,iF]]
      end  
      RhsLoc += 1 / detJ * Weights[i] * qLoc * ((DF * uLoc)' * (DF * fTRef[:,:,i]))'
    end
    @views @. Rhs[FeT.Glob[:,iF]] += RhsLoc
  end
end  

function GradKinHeight!(backend,FTB,Rhs,u,uFe::HDivElement,h,hFe::ScalarElement,FeT::HDivElement,
  Grid,QuadOrd,Jacobi)
  QQ = FEMSei.QuadRule{FTB}(Grid.Type,backend,QuadOrd)
  Weights = QQ.Weights
  Points = QQ.Points
  fTRef  = zeros(hFe.Comp,FeT.DoF,length(Weights))
  hfFRef  = zeros(hFe.Comp,hFe.DoF,length(Weights))
  ufFRef  = zeros(uFe.Comp,uFe.DoF,length(Weights))

  for i = 1 : length(Weights)
    for iComp = 1 : hFe.Comp
      for iD = 1 : FeT.DoF
        fTRef[iComp,iD,i] = FeT.Divphi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
    for iComp = 1 : hFe.Comp
      for iD = 1 : hFe.DoF
        hfFRef[iComp,iD,i] = hFe.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
    for iComp = 1 : uFe.Comp
      for iD = 1 : uFe.DoF
        ufFRef[iComp,iD,i] = uFe.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  RhsLoc = zeros(FeT.DoF)
  uLoc = zeros(uFe.Comp)
  for iF = 1 : Grid.NumFaces
    RhsLoc .= 0
    for i = 1 : length(Weights)
      DF, detJ = Jacobi(Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      hLoc = hfFRef[1,1,i] * h[hFe.Glob[1,iF]]
      for j = 2 : hFe.DoF 
        hLoc = hfFRef[1,j,i] * h[hFe.Glob[j,iF]]
      end  
      uLoc .= ufFRef[:,1,i] * u[uFe.Glob[1,iF]]
      for j = 2 : FeT.DoF 
        uLoc .+= ufFRef[:,j,i] * u[uFe.Glob[j,iF]]
      end  
      uLoc3 = 1 / detJ * DF * uLoc
      RhsLoc += Weights[i] * (uLoc3'*uLoc3 + hLoc) * fTRef[:,:,i]'
    end
    @views @. Rhs[FeT.Glob[:,iF]] += RhsLoc
  end
end  

function DivHeight!(backend,FTB,Rhs,u,uFe::HDivElement,h,hFe::ScalarElement,FeT::ScalarElement,
  Grid,QuadOrd,Jacobi)
  QQ = FEMSei.QuadRule{FTB}(Grid.Type,backend,QuadOrd)
  Weights = QQ.Weights
  Points = QQ.Points
  fTRef  = zeros(hFe.Comp,FeT.DoF,length(Weights))
  ufFRef  = zeros(FeT.Comp,uFe.DoF,length(Weights))

  for i = 1 : length(Weights)
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeT.DoF
        fTRef[iComp,iD,i] = FeT.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
    for iComp = 1 : FeT.Comp
      for iD = 1 : uFe.DoF
        ufFRef[iComp,iD,i] = uFe.Divphi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  hfFRef = fTRef
  RhsLoc = zeros(FeT.DoF)
  for iF = 1 : Grid.NumFaces
    RhsLoc .= 0
    for i = 1 : length(Weights)
      hLoc = hfFRef[1,1,i] * h[hFe.Glob[1,iF]]
      for j = 2 : hFe.DoF 
        hLoc += hfFRef[1,j,i] * h[hFe.Glob[j,iF]]
      end  
      uLoc = ufFRef[1,1,i] * u[uFe.Glob[1,iF]]
      for j = 2 : uFe.DoF 
        uLoc += ufFRef[1,j,i] * u[uFe.Glob[j,iF]]
      end  
      RhsLoc += Weights[i] * uLoc * hLoc * fTRef[:,:,i]'
    end
    @views @. Rhs[FeT.Glob[:,iF]] += RhsLoc
  end
end  

function DivMatrix(backend,FTB,FeF::HDivElement,FeT::ScalarElement,Grid,QuadOrd,Jacobi)
  QQ = FEMSei.QuadRule{FTB}(Grid.Type,backend,QuadOrd)
  Weights = QQ.Weights
  Points = QQ.Points
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

  for iF = 1 : Grid.NumFaces
    DivLoc .= 0
    for i = 1 : length(Weights)
      DF, detJ = Jacobi(Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      DivLoc += sign(detJ) * Weights[i] * (fTRef[:,:,i]' * fFRef[:,:,i])
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

function DivMatrix(backend,FTB,FeF::HDivKiteDElement,FeT::ScalarElement,Grid,QuadOrd,Jacobi)
  QQ = FEMSei.QuadRule{FTB}(Grid.Type,backend,QuadOrd)
  Weights = QQ.Weights
  Points = QQ.Points
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
  QQL = FEMSei.QuadRule{FTB}(Grids.Line(),backend,QuadOrd)
  WeightsL = QQL.Weights
  PointsL = QQL.Points
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

  for iF = 1 : Grid.NumFaces
    DivLoc .= 0
    for i = 1 : length(Weights)
      DF, detJ = Jacobi(Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      DivLoc += sign(detJ) * Weights[i] * (fTRef[:,:,i]' * fFRef[:,:,i])
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
  QQ = FEMSei.QuadRule{FTB}(Grid.Type,backend,QuadOrd)
  Weights = QQ.Weights
  Points = QQ.Points
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
  QQL = FEMSei.QuadRule{FTB}(Grids.Line(),backend,QuadOrd)
  WeightsL = QQL.Weights
  PointsL = QQL.Points
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

  for iF = 1 : Grid.NumFaces
    CurlLoc .= 0
    for i = 1 : length(Weights)
      _, detJ = Jacobi(Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      CurlLoc += -sign(detJ) * Weights[i] * (fTRef[:,:,i]' * fFRef[:,:,i])
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
  QQ = FEMSei.QuadRule{FTB}(Grid.Type,backend,QuadOrd)
  Weights = QQ.Weights
  Points = QQ.Points
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

  for iF = 1 : Grid.NumFaces
    GradLoc .= 0
    for i = 1 : length(Weights)
      DF, detJ = Jacobi(Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      GradLoc += sign(detJ) * Weights[i] * (fTRef[:,:,i]' * fFRef[:,:,i])
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
  QQ = FEMSei.QuadRule{FTB}(Grid.Type,backend,QuadOrd)
  Weights = QQ.Weights
  Points = QQ.Points
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
  QQL = FEMSei.QuadRule{FTB}(Grids.Line(),backend,QuadOrd)
  WeightsL = QQL.Weights
  PointsL = QQL.Points
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

  for iF = 1 : Grid.NumFaces
    GradLoc .= 0
    for i = 1 : length(Weights)
      DF, detJ = Jacobi(Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      GradLoc += sign(detJ) * Weights[i] * (fTRef[:,:,i]' * fFRef[:,:,i])
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

    
