function MassMatrix(backend,FTB,Fe::HCurlElement,Grid,QuadOrd,Jacobi)
  QQ = FEMSei.QuadRule{FTB}(Fe.Type,backend,QuadOrd)
  Weights = QQ.Weights
  Points = QQ.Points
  fRef  = zeros(Fe.Comp,Fe.DoF,length(Weights))
  DF  = zeros(Fe.Comp,Fe.DoF,length(Weights))

  for i = 1 : length(Weights)
    for iComp = 1 : Fe.Comp
      for iD = 1 : Fe.DoF
        fRef[iComp,iD,i] = Fe.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end

  RowInd = Int64[]
  ColInd = Int64[]
  Val = Float64[]

  MLoc = zeros(Fe.DoF,Fe.DoF)
  for iF = 1 : Grid.NumFaces
    MLoc .= 0  
    for i = 1 : length(Weights)
      DF, detJ, invPDF = Jacobi(Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      fLoc = invPDF * fRef[:, :, i]
      MLoc = MLoc + abs(detJ) * Weights[i]*(fLoc'*fLoc)
    end
    for j = 1 : size(MLoc,2)
      for i = 1 : size(MLoc,1)
        if abs(MLoc[i,j]) > 1.e-6 
          push!(RowInd,Fe.Glob[i,iF])
          push!(ColInd,Fe.Glob[j,iF])
          push!(Val,MLoc[i,j])
        end  
      end
    end
  end
  M = sparse(RowInd, ColInd, Val)
  return M
end

function MassMatrix(backend,FTB,Fe::HDivElement,Grid,QuadOrd,Jacobi)
  QQ = FEMSei.QuadRule{FTB}(Fe.Type,backend,QuadOrd)
  Weights = QQ.Weights
  Points = QQ.Points
  fRef  = zeros(Fe.Comp,Fe.DoF,length(Weights))
  DF  = zeros(Fe.Comp,Fe.DoF,length(Weights))

  for i = 1 : length(Weights)
    for iComp = 1 : Fe.Comp
      for iD = 1 : Fe.DoF
        fRef[iComp,iD,i] = Fe.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end

  RowInd = Int64[]
  ColInd = Int64[]
  Val = Float64[]

  MLoc = zeros(Fe.DoF,Fe.DoF)
  for iF = 1 : Grid.NumFaces
    MLoc .= 0  
    for i = 1 : length(Weights)
      DF, detJ = Jacobi(Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      fLoc = DF * fRef[:, :, i]
      MLoc = MLoc + 1 / abs(detJ) * Weights[i]*(fLoc'*fLoc)
    end
    for j = 1 : size(MLoc,2)
      for i = 1 : size(MLoc,1)
        if abs(MLoc[i,j]) > 1.e-6 
          push!(RowInd,Fe.Glob[i,iF])
          push!(ColInd,Fe.Glob[j,iF])
          push!(Val,MLoc[i,j])
        end  
      end
    end
  end
  M = sparse(RowInd, ColInd, Val)
  return M
end

function MassMatrix(backend,FTB,Fe::ScalarElement,Grid,QuadOrd,Jacobi)
  QQ = FEMSei.QuadRule{FTB}(Fe.Type,backend,QuadOrd)
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

  RowInd = Int64[]
  ColInd = Int64[]
  Val = Float64[]

  for iF = 1 : Grid.NumFaces
    MLoc = zeros(Fe.DoF,Fe.DoF)
      for i = 1 : length(Weights)
        _,detJ = Jacobi(Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF],Grid)
        fLoc = fRef[:, :, i]
        MLoc = MLoc + abs(detJ) * Weights[i] *(fLoc' * fLoc)
      end
      for j = 1 : size(MLoc,2)
        for i = 1 : size(MLoc,1)
          if abs(MLoc[i,j]) > 0  
            push!(RowInd,Fe.Glob[i,iF])
            push!(ColInd,Fe.Glob[j,iF])
            push!(Val,MLoc[i,j])
          end
        end
      end
  end
  M = sparse(RowInd, ColInd, Val)
  return M
end


