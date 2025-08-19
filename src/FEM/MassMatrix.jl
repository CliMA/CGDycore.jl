"""
  MassMatrix H(curl) elements

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

"""
  MassMatrix H(div) elements

Assembles the global mass matrix for a finite element space defined by `Fe` over the mesh `Grid` using a specified quadrature order `QuadOrd` and a Jacobian computation function `Jacobi`.

# Arguments
- `backend`: Computational backend to use (not used directly in this function).
- `FTB`: Finite element type or basis (not used directly in this function).
- `Fe::HDivElement`: Finite element object containing basis functions, degrees of freedom, and global indexing.
- `Grid`: Mesh/grid object containing face information and geometry.
- `QuadOrd`: Integer specifying the order of the quadrature rule.
- `Jacobi`: Function to compute the Jacobian, its determinant, and pseudo-inverse for mapping reference to physical elements.

# Returns
- `M`: Sparse matrix representing the assembled global mass matrix.

# Details
- For each face in the mesh, the function computes the local mass matrix using numerical quadrature.
- Basis functions are evaluated at quadrature points and mapped to the physical element using the Jacobian.
- The local mass matrix contributions are assembled into the global sparse matrix using the global degree of freedom indices.
- Entries with absolute value less than `1e-12` are ignored for sparsity.

# Notes
- The function assumes that the finite element basis functions are provided as callable objects in `Fe.phi`.
- The function is tailored for `H(div)`-conforming elements, where basis functions are vector-valued.
"""
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

"""
  MassMatrix Scalar elements

Assembles the global mass matrix for a finite element mesh.

# Arguments
- `backend`: Computational backend to use (e.g., CPU, GPU).
- `FTB`: Finite element type or basis information.
- `Fe::ScalarElement`: Scalar finite element object containing basis functions and related data.
- `Grid`: Grid or mesh structure containing information about faces and connectivity.
- `QuadOrd`: Quadrature order for numerical integration.
- `Jacobi`: Function to compute the Jacobian, its determinant, and its pseudo-inverse for mapping reference to physical elements.

# Returns
- `M`: Sparse matrix representing the assembled global mass matrix.

# Description
This function computes the mass matrix by looping over all faces in the mesh, evaluating basis functions at quadrature points, and assembling local contributions into the global sparse matrix. The Jacobian is used to map reference element integrals to the physical domain. Only nonzero entries are stored in the sparse matrix.

# Notes
- The function assumes that `Fe.phi` contains callable basis functions.
- The global degree of freedom mapping is provided by `Fe.Glob`.
- The quadrature rule is determined by `FEMSei.QuadRule`.
"""
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

"""
  MassMatrix Scalar elements with weights

Assembles the mass matrix for a finite element method (FEM) problem using the provided scalar elements, quadrature order, and grid information.

# Arguments
- `backend`: Computational backend to use (e.g., CPU, GPU).
- `FTB`: Finite element basis or related data structure.
- `Fe::ScalarElement`: Finite element object representing the trial space.
- `w`: Vector of weights or coefficients for the weighted mass matrix.
- `wFe::ScalarElement`: Finite element object representing the weighting/test space.
- `Grid`: Grid or mesh data structure containing faces and connectivity.
- `QuadOrd`: Integer specifying the quadrature order for numerical integration.
- `Jacobi`: Function to compute the Jacobian, its determinant, and pseudo-inverse for mapping reference to physical elements.

# Returns
- `M`: Sparse matrix representing the assembled mass matrix.

# Description
The function computes the local mass matrix for each face in the grid using numerical quadrature, basis functions, and the provided weights. It then assembles these local matrices into a global sparse mass matrix. The Jacobian is used to map quadrature points from the reference element to the physical element, and the determinant is used to scale the integration weights.

# Notes
- The function assumes that the basis functions (`phi`) are provided as callable objects in `Fe` and `wFe`.
- The global indices for assembling the sparse matrix are taken from `Fe.Glob` and `wFe.Glob`.
- The function supports assembling weighted mass matrices by incorporating the vector `w` and its associated basis functions.
"""
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

"""
  MassMatrix Vector elements

Assembles the global mass matrix for a finite element discretization.

# Arguments
- `backend`: Computational backend to use (e.g., CPU, GPU).
- `FTB`: Finite element type or basis information.
- `Fe::VectorElement`: Finite element object containing basis functions, degrees of freedom, and global mapping.
- `Grid`: Grid or mesh structure containing face and connectivity information.
- `QuadOrd`: Quadrature order for numerical integration.
- `Jacobi`: Function to compute the Jacobian, its determinant, and pseudo-inverse for mapping reference to physical elements.

# Returns
- `M`: Sparse matrix representing the assembled global mass matrix.

# Description
For each face in the grid, this function computes the local mass matrix using numerical quadrature and the provided basis functions. The local contributions are then assembled into the global sparse mass matrix using the global degree of freedom mapping.

# Notes
- The function assumes 2D elements and vector-valued basis functions.
- The quadrature rule and basis function evaluations are performed on the reference element and mapped to the physical element using the Jacobian.
- Only nonzero entries are stored in the sparse matrix.
"""
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




