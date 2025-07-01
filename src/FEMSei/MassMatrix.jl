"""
  MassMatrix

Assembles the global mass matrix for a finite element discretization using H(curl) elements.

# Arguments
- `backend`: Computational backend to use (e.g., CPU, GPU).
- `FTB`: Finite element basis or transformation object.
- `Fe::HCurlElement`: The finite element object containing basis functions and related data for H(curl) elements.
- `Grid`: The computational grid or mesh structure containing faces and geometry information.
- `QuadOrd`: The order of the quadrature rule to use for numerical integration.
- `Jacobi`: Function that computes the Jacobian, its determinant, and pseudo-inverse for mapping reference to physical elements.

# Returns
- `M`: The assembled sparse global mass matrix as a `SparseMatrixCSC`.

# Description
For each face in the grid, the function computes the local mass matrix by integrating the inner product of the mapped basis functions over the reference element using a specified quadrature rule. The local contributions are then assembled into the global sparse mass matrix. Only non-negligible entries (greater than `1.e-6` in magnitude) are included in the global matrix.

# Notes
- The function assumes that the basis functions and their mapping are provided in the `Fe` object.
- The Jacobian and its pseudo-inverse are used to map basis functions from the reference to the physical element.
- The function is tailored for H(curl) elements and may require adaptation for other element types.
"""
function MassMatrix(backend,FTB,Fe::HCurlElement,Grid,QuadOrd,Jacobi)
  NumQuad,Weights,Points = FEMSei.QuadRule(Fe.Type,QuadOrd)
  fRef  = zeros(Fe.Comp,Fe.DoF,NumQuad)
  DF  = zeros(Fe.Comp,Fe.DoF,NumQuad)

  for iQ = 1 : NumQuad
    for iComp = 1 : Fe.Comp
      for iD = 1 : Fe.DoF
        fRef[iComp,iD,iQ] = Fe.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
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
    for iQ = 1 : NumQuad
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,Points[iQ,1],Points[iQ,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      fLoc = pinvDF * fRef[:, :, iQ]
      MLoc = MLoc + abs(detDFLoc) * Weights[iQ]*(fLoc'*fLoc)
    end
    for j = 1 : Fe.DoF
      for i = 1 : Fe.DoF
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
  fRef  = zeros(Fe.Comp,Fe.DoF,NumQuad)
  DF  = zeros(Fe.Comp,Fe.DoF,NumQuad)

  for i = 1 : NumQuad
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
    for i = 1 : NumQuad
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      fLoc = DF * fRef[:, :, i]
      MLoc = MLoc + 1 / abs(detDFLoc) * Weights[i]*(fLoc'*fLoc)
    end
    for j = 1 : size(MLoc,2)
      for i = 1 : size(MLoc,1)
        if abs(MLoc[i,j]) > 1.e-12 
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
  fRef  = zeros(Fe.Comp,Fe.DoF,NumQuad)

  for i = 1 : NumQuad
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
      for i = 1 : NumQuad
        Jacobi(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF],Grid)
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

function MassMatrix(backend,FTB,Fe::ScalarElement,w,wFe::ScalarElement,Grid,QuadOrd,Jacobi)
  NumQuad,Weights,Points = FEMSei.QuadRule(Fe.Type,QuadOrd)
  fRef  = zeros(Fe.Comp,Fe.DoF,NumQuad)
  fwRef  = zeros(wFe.Comp,wFe.DoF,NumQuad)

  for i = 1 : NumQuad
    for iComp = 1 : Fe.Comp
      for iD = 1 : Fe.DoF
        fRef[iComp,iD,i] = Fe.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
    for iComp = 1 : wFe.Comp
      for iD = 1 : wFe.DoF
        fwRef[iComp,iD,i] = wFe.phi[iD,iComp](Points[i,1],Points[i,2])
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
  wLoc = zeros(wFe.DoF)

  for iF = 1 : Grid.NumFaces
    MLoc = zeros(Fe.DoF,Fe.DoF)
    for iDoF = 1 : wFE.DoF
      ind = wFe.Glob[iDoF,iF]
      wLoc[iDoF] = w[ind]
    end  
    for i = 1 : NumQuad
      wwLoc = 0.0  
      for iDoF = 1 : wFE.DoF
        wwLoc += wLoc[iDoF] * fwRef[1,iD,i]
      end  
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF],Grid)
      detDFLoc = detDF[1]
      fLoc = fRef[:, :, i]
      MLoc = MLoc + abs(detDFLoc) * Weights[i] * wwLoc * (fLoc' * fLoc)
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

function MassMatrix(backend,FTB,Fe::VectorElement,Grid,QuadOrd,Jacobi)
  NumQuad,Weights,Points = FEMSei.QuadRule(Fe.Type,QuadOrd)
  fRef  = zeros(Fe.Comp,Fe.DoF,NumQuad)

  for i = 1 : NumQuad
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
      for i = 1 : NumQuad
        Jacobi(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF],Grid)
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




