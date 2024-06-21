function MassMatrix(backend,FTB,Fe::HCurlElement,Grid,QuadOrd,Jacobi)
  NumQuad,Weights,Points = FEMSei.QuadRule(Fe.Type,QuadOrd)
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
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  for iF = 1 : Grid.NumFaces
    MLoc .= 0  
    for i = 1 : length(Weights)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      fLoc = pinvDF * fRef[:, :, i]
      MLoc = MLoc + abs(detDFLoc) * Weights[i]*(fLoc'*fLoc)
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
  NumQuad,Weights,Points = FEMSei.QuadRule(Fe.Type,QuadOrd)
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
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  for iF = 1 : Grid.NumFaces
    MLoc .= 0  
    for i = 1 : length(Weights)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      fLoc = DF * fRef[:, :, i]
      MLoc = MLoc + 1 / abs(detDFLoc) * Weights[i]*(fLoc'*fLoc)
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
  NumQuad,Weights,Points = FEMSei.QuadRule(Fe.Type,QuadOrd)
  fRef  = zeros(Fe.Comp,Fe.DoF,length(Weights))
  @show Weights
  @show Points
  @show NumQuad

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
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)

  for iF = 1 : Grid.NumFaces
    MLoc = zeros(Fe.DoF,Fe.DoF)
      for i = 1 : length(Weights)
        Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF],Grid)
        detDFLoc = detDF[1]
        fLoc = fRef[:, :, i]
        MLoc = MLoc + abs(detDFLoc) * Weights[i] *(fLoc' * fLoc)
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


