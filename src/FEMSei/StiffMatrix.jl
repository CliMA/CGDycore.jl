"""
  VortCrossVel

Assembles the right-hand side (Rhs) vector for the vorticity-cross-velocity term in a finite element method (FEM) context, specifically for mixed HDiv and scalar elements.

# Arguments
- `backend`: Computational backend or context (e.g., CPU/GPU).
- `FTB`: (Unused in this function) Placeholder for future extensions or backend-specific data.
- `Rhs`: The right-hand side vector to be assembled (modified in-place).
- `u`: Coefficient vector for the velocity field.
- `uFe::HDivConfElement`: HDiv-conforming finite element description for the velocity field.
- `q`: Coefficient vector for the scalar field (e.g., vorticity or pressure).
- `qFe::ScalarElement`: Scalar finite element description for the scalar field.
- `FeT::HDivConfElement`: HDiv-conforming finite element description for the test functions.
- `Grid`: Grid or mesh data structure containing geometry and connectivity.
- `QuadOrd`: Quadrature order for numerical integration.
- `Jacobi`: Function to compute the Jacobian, its determinant, and pseudo-inverse for mapping reference to physical elements.

# Description
For each face in the grid, this function:
- Evaluates basis functions and their values at quadrature points.
- Computes the local contributions to the right-hand side vector by integrating the vorticity-cross-velocity term using the provided quadrature rule.
- Applies the transformation from reference to physical element using the Jacobian.
- Accumulates the local contributions into the global `Rhs` vector.

# Notes
- Assumes that the basis functions and their global-to-local mappings are provided in the finite element objects.
- Uses in-place operations and preallocated arrays for efficiency.
- The function is intended for use in mixed finite element formulations involving HDiv and scalar elements.

# See also
- `HDivConfElement`
- `ScalarElement`
- `QuadRule`
- `Jacobi`
"""

function VortCrossVel!(backend,FTB,Rhs,u,uFe::HDivConfElement,q,qFe::ScalarElement,FeT::HDivConfElement,
  Grid,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(Grid.Type,QuadOrd)
  fTRef  = zeros(FeT.Comp,FeT.DoF,NumQuad)
  qfFRef  = zeros(qFe.Comp,qFe.DoF,NumQuad)

  @inbounds for i = 1 : NumQuad
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : FeT.DoF
        fTRef[iComp,iD,i] = FeT.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
    @inbounds for iComp = 1 : qFe.Comp
      @inbounds for iD = 1 : qFe.DoF
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
    @inbounds for i = 1 : NumQuad
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      qLoc = qfFRef[1,1,i] * q[qFe.Glob[1,iF]]
      @inbounds for j = 2 : qFe.DoF 
        qLoc += qfFRef[1,j,i] * q[qFe.Glob[j,iF]]
      end  
      uLoc .= [-ufFRef[2,1,i],ufFRef[1,1,i]] * u[uFe.Glob[1,iF]]
      @inbounds for j = 2 : FeT.DoF 
        uLoc .+= [-ufFRef[2,j,i],ufFRef[1,j,i]] * u[uFe.Glob[j,iF]]
      end  
      RhsLoc += 1 / detDFLoc * Weights[i] * qLoc * ((DF * uLoc)' * (DF * fTRef[:,:,i]))'
    end
    @views @. Rhs[FeT.Glob[:,iF]] += RhsLoc
  end
end  

"""
  DivMatrix(backend, FTB, FeF::HDivConfElement, FeT::ScalarElement, Grid, QuadOrd, Jacobi)

Assembles the divergence matrix for a finite element method (FEM) problem, coupling an `HDivConfElement` (face element) and a `ScalarElement` (cell element) over a given grid.

# Arguments
- `backend`: Computational backend to use (e.g., CPU, GPU).
- `FTB`: Additional backend or finite element type information.
- `FeF::HDivConfElement`: The face finite element with divergence basis functions.
- `FeT::ScalarElement`: The scalar finite element with basis functions.
- `Grid`: The grid structure containing mesh and face information.
- `QuadOrd`: Quadrature order for numerical integration.
- `Jacobi`: Function to compute the Jacobian, its determinant, and pseudo-inverse for mapping reference to physical elements.

# Returns
- `Div`: A sparse matrix representing the divergence operator assembled over the grid.

# Details
- Uses quadrature to integrate the product of test and trial functions (and their divergences) over each face in the grid.
- Handles mapping from reference to physical elements using the provided `Jacobi` function.
- Assembles the global sparse matrix using local contributions from each face.

# Notes
- Assumes that `FeT.phi` and `FeF.Divphi` are arrays of functions for evaluating basis functions and their divergences at quadrature points.
- The orientation of each face is taken into account via `Grid.Faces[iF].Orientation`.
"""

function DivMatrix(backend,FTB,FeF::HDivConfElement,FeT::ScalarElement,Grid,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(FeT.Type,QuadOrd)
  fFRef  = zeros(FeT.Comp,FeF.DoF,NumQuad)
  fTRef  = zeros(FeT.Comp,FeT.DoF,NumQuad)

  @inbounds for i = 1 : NumQuad
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : FeT.DoF
        fTRef[iComp,iD,i] = FeT.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : FeF.DoF
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
    @inbounds for i = 1 : NumQuad
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      DivLoc += Grid.Faces[iF].Orientation * Weights[i] * (fTRef[:,:,i]' * fFRef[:,:,i])
    end
    @inbounds for j = 1 : size(DivLoc,2)
      @inbounds for i = 1 : size(DivLoc,1)
        push!(RowInd,FeT.Glob[i,iF])
        push!(ColInd,FeF.Glob[j,iF])
        push!(Val,DivLoc[i,j])
      end
    end
  end
  Div = sparse(RowInd, ColInd, Val)
  return Div
end

"""
  LaplMatrix(backend, FTB, FeF::ScalarElement, FeT::ScalarElement, Grid, QuadOrd, Jacobi)

Assembles the global Laplacian (stiffness) matrix for a finite element method (FEM) discretization.

# Arguments
- `backend`: Computational backend to use (not used directly in this function).
- `FTB`: Additional backend or transformation information (not used directly in this function).
- `FeF::ScalarElement`: Finite element object for the test space (provides basis function gradients and global DoF mapping).
- `FeT::ScalarElement`: Finite element object for the trial space (provides basis function gradients and global DoF mapping).
- `Grid`: Grid or mesh object containing information about the domain, faces, and geometry.
- `QuadOrd`: Quadrature order to use for numerical integration.
- `Jacobi`: Function that computes the Jacobian matrix, its determinant, and its pseudo-inverse for a given face and reference point.

# Returns
- `Lapl`: Sparse matrix representing the assembled global Laplacian (stiffness) matrix.

# Description
This function computes the global Laplacian matrix by looping over all faces of the mesh, performing numerical integration using the specified quadrature rule, and assembling local element matrices into the global sparse matrix. The gradients of the basis functions are mapped from the reference element to the physical element using the Jacobian and its pseudo-inverse.

# Notes
- The function assumes 2D elements and uses two components for gradients.
- The local-to-global mapping for degrees of freedom is provided by `FeT.Glob` and `FeF.Glob`.
- The function is optimized for performance using `@inbounds` and preallocated arrays.
"""

function LaplMatrix(backend,FTB,FeF::ScalarElement,FeT::ScalarElement,Grid,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(Grid.Type,QuadOrd)
  fFRef  = zeros(2,FeF.DoF,NumQuad)
  fTRef  = zeros(2,FeT.DoF,NumQuad)

  @inbounds for i = 1 : NumQuad
    @inbounds for iComp = 1 : 2
      @inbounds for iD = 1 : FeT.DoF
        fTRef[iComp,iD,i] = FeT.Gradphi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
    @inbounds for iComp = 1 : 2
      @inbounds for iD = 1 : FeF.DoF
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
    @inbounds for i = 1 : NumQuad
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      LaplLoc +=  -detDFLoc * Weights[i] * ((pinvDF * fTRef[:,:,i])' * 
        (pinvDF * fFRef[:,:,i]))
    end
    @inbounds for j = 1 : size(LaplLoc,2)
      @inbounds for i = 1 : size(LaplLoc,1)
        push!(RowInd,FeT.Glob[i,iF])
        push!(ColInd,FeF.Glob[j,iF])
        push!(Val,LaplLoc[i,j])
      end
    end
  end
  Lapl = sparse(RowInd, ColInd, Val)
  return Lapl
end

"""
  DivRhs!(backend, FTB, Div, FeT::ScalarElement, u, uFeF::HDivConfElement,
      Grid, ElemType::Grids.ElementType, QuadOrd, Jacobi)

Assembles the right-hand side vector for the divergence operator in a finite element method (FEM) context, specifically for scalar test functions and H(div)-conforming trial functions.

# Arguments
- `backend`: Computational backend to use (e.g., CPU, GPU).
- `FTB`: Finite element basis or related data structure.
- `Div`: Output vector to accumulate the divergence right-hand side contributions (modified in-place).
- `FeT::ScalarElement`: Scalar finite element type for the test space.
- `u`: Global vector of degrees of freedom for the H(div)-conforming field.
- `uFeF::HDivConfElement`: H(div)-conforming finite element type for the trial space.
- `Grid`: Grid or mesh data structure containing face and orientation information.
- `ElemType::Grids.ElementType`: Element type descriptor for quadrature rule selection.
- `QuadOrd`: Quadrature order for numerical integration.
- `Jacobi`: Jacobian or transformation data (not used directly in this function).

# Description
For each face in the grid, this function:
1. Evaluates basis functions and their divergences at quadrature points.
2. Computes local contributions to the divergence right-hand side by integrating the product of test and trial functions over each face, taking into account face orientation and quadrature weights.
3. Assembles these local contributions into the global `Div` vector.

# Notes
- The function assumes 2D reference elements and that basis function evaluations are provided as callable objects.
- The function is intended for use in mixed finite element methods, such as those involving H(div) spaces (e.g., Raviart-Thomas elements).
- All operations are performed in-place for efficiency.

# In-place Behavior
- The `Div` vector is updated in-place with the assembled right-hand side contributions.

# Performance
- Uses `@inbounds` for performance; ensure input data is valid to avoid out-of-bounds errors.

"""

function DivRhs!(backend,FTB,Div,FeT::ScalarElement,u,uFeF::HDivConfElement,
  Grid,ElemType::Grids.ElementType,QuadOrd,Jacobi)

  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  uFRef  = zeros(FeT.Comp,uFeF.DoF,NumQuad)
  fTRef  = zeros(FeT.Comp,FeT.DoF,NumQuad)

  @inbounds for iQ = 1 : NumQuad
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iDoF = 1 : FeT.DoF
        fTRef[iComp,iDoF,iQ] = FeT.phi[iDoF,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iDoF = 1 : uFeF.DoF
        uFRef[iComp,iDoF,iQ] = uFeF.Divphi[iDoF,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
  end
  DivLoc = zeros(FeT.DoF)
  uLoc = zeros(uFeF.DoF)

  
  @inbounds for iF = 1 : Grid.NumFaces
    DivLoc .= 0
    @inbounds for iDoF  = 1 : uFeF.DoF
      ind = uFeF.Glob[iDoF,iF]
      uLoc[iDoF] = u[ind]
    end
    @inbounds for iQ = 1 : NumQuad
      uFRefLoc = 0.0
      @inbounds for iDoF = 1 : uFeF.DoF
        uFRefLoc += uFRef[1,iDoF,iQ] * uLoc[iDoF]  
      end  
      @inbounds for iDoF = 1 : FeT.DoF
        DivLoc[iDoF] -= Grid.Faces[iF].Orientation * Weights[iQ] * fTRef[1,iDoF,iQ] * uFRefLoc 
      end  
    end
    @inbounds for iDoF  = 1 : FeT.DoF
      ind = FeT.Glob[iDoF,iF]
      Div[ind] += DivLoc[iDoF]
    end
  end
end

"""
  GradRhs!(backend, FTB, Grad, h, hFeF::ScalarElement, FeT::HDivConfElement, Grid, ElemType::Grids.ElementType,
       QuadOrd, Jacobi)

Assembles the right-hand side (RHS) contribution for the gradient operator in a finite element method (FEM) context, specifically for mixed formulations involving scalar and HDiv-conforming elements.

# Arguments
- `backend`: Computational backend to use (e.g., CPU, GPU).
- `FTB`: (Purpose not specified in the code snippet; likely a data structure for finite element assembly).
- `Grad`: Array to accumulate the assembled gradient RHS contributions (modified in-place).
- `h`: Array of scalar field degrees of freedom (DoFs).
- `hFeF::ScalarElement`: Scalar finite element descriptor for the field `h`.
- `FeT::HDivConfElement`: HDiv-conforming finite element descriptor for the test space.
- `Grid`: Grid or mesh data structure containing face and orientation information.
- `ElemType::Grids.ElementType`: Element type descriptor (e.g., triangle, quadrilateral).
- `QuadOrd`: Quadrature order for numerical integration.
- `Jacobi`: (Purpose not specified in the code snippet; possibly a flag or structure for Jacobian handling).

# Description
For each face in the grid, this function:
- Evaluates basis functions and their divergences at quadrature points.
- Computes local contributions to the gradient RHS using quadrature.
- Assembles these local contributions into the global `Grad` array.

The function is optimized with `@inbounds` for performance and assumes that the basis function arrays and grid structures are properly initialized.

# Notes
- The function modifies `Grad` in-place.
- The function assumes that the basis function arrays (`phi`, `Divphi`) are callable with reference coordinates.
- The purpose of `backend`, `FTB`, and `Jacobi` is not fully specified in the provided code.

# See Also
- `ScalarElement`
- `HDivConfElement`
- `QuadRule`
"""

function GradRhs!(backend,FTB,Grad,h,hFeF::ScalarElement,FeT::HDivConfElement,Grid,ElemType::Grids.ElementType,
  QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  fFRef  = zeros(hFeF.Comp,hFeF.DoF,NumQuad)
  fTRef  = zeros(hFeF.Comp,FeT.DoF,NumQuad)

  @inbounds for i = 1 : NumQuad
    @inbounds for iComp = 1 : hFeF.Comp
      @inbounds for iDoF = 1 : FeT.DoF
        fTRef[iComp,iDoF,i] = FeT.Divphi[iDoF,iComp](Points[i,1],Points[i,2])
      end
    end
    @inbounds for iComp = 1 : hFeF.Comp
      @inbounds for iDoF = 1 : hFeF.DoF
        fFRef[iComp,iDoF,i] = hFeF.phi[iDoF,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  GradLoc = zeros(FeT.DoF)
  hLoc = zeros(hFeF.DoF)

  @inbounds for iF = 1 : Grid.NumFaces
    GradLoc .= 0
    @inbounds for iDoF  = 1 : hFeF.DoF
      ind = hFeF.Glob[iDoF,iF]
      hLoc[iDoF] = h[ind]
    end
    @inbounds for iQ = 1 : NumQuad
      fFRefLoc = 0.0
      @inbounds for iDoF = 1 : hFeF.DoF
        fFRefLoc += fFRef[1,iDoF,iQ] * hLoc[iDoF]  
      end  
      @inbounds for iDoF = 1 : FeT.DoF
        GradLoc[iDoF] += Grid.Faces[iF].Orientation * Weights[iQ] * fTRef[1,iDoF,iQ] * fFRefLoc 
      end  
    end
    @inbounds for iDoF  = 1 : FeT.DoF
      ind = FeT.Glob[iDoF,iF]
      Grad[ind] += GradLoc[iDoF]
    end
  end
end

"""
  CurlMatrix(backend, FTB, FeF::HCurlConfElement, FeT::ScalarElement, Grid, QuadOrd, Jacobi)

Assembles the curl matrix for a finite element method (FEM) problem involving an `HCurlConfElement` (typically a vector-valued element with curl-conforming basis functions) and a `ScalarElement` (typically a scalar-valued element). The resulting sparse matrix represents the action of the curl operator in the FEM discretization.

# Arguments
- `backend`: Computational backend to use (e.g., CPU, GPU, etc.).
- `FTB`: Additional backend or finite element type information (usage depends on context).
- `FeF::HCurlConfElement`: The finite element object representing the curl-conforming (vector-valued) basis functions.
- `FeT::ScalarElement`: The finite element object representing the scalar-valued basis functions.
- `Grid`: The grid or mesh structure containing information about faces, orientations, and global DoF mappings.
- `QuadOrd`: The order of the quadrature rule to use for numerical integration.
- `Jacobi`: Jacobi data or transformation (usage depends on context).

# Returns
- `Curl`: A sparse matrix representing the assembled curl operator between the given finite element spaces.

# Details
- Uses numerical quadrature to integrate the product of test and trial functions (and their curls) over each face of the grid.
- Handles orientation of faces and global degree-of-freedom (DoF) mapping.
- Efficiently assembles the matrix in sparse format.

# Notes
- Assumes that `FeT.phi` and `FeF.Curlphi` are arrays of function handles for evaluating basis functions and their curls at quadrature points.
- The function is performance-optimized using `@inbounds` and preallocated arrays.
"""

#HDiv->HCurl
function CurlMatrix(backend,FTB,FeF::HCurlConfElement,FeT::ScalarElement,Grid,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(FeT.Type,QuadOrd)
  fFRef  = zeros(FeT.Comp,FeF.DoF,NumQuad)
  fTRef  = zeros(FeT.Comp,FeT.DoF,NumQuad)

  @inbounds for i = 1 : NumQuad
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : FeT.DoF
        fTRef[iComp,iD,i] = FeT.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : FeF.DoF
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
    @inbounds for i = 1 : NumQuad
      CurlLoc += Grid.Faces[iF].Orientation * Weights[i] * (fTRef[:,:,i]' * fFRef[:,:,i])
    end
    @inbounds for j = 1 : size(CurlLoc,2)
      @inbounds for i = 1 : size(CurlLoc,1)
        push!(RowInd,FeT.Glob[i,iF])
        push!(ColInd,FeF.Glob[j,iF])
        push!(Val,CurlLoc[i,j])
      end
    end
  end
  Curl = sparse(RowInd, ColInd, Val)
  return Curl
end

"""
  DivMatrix(backend, FTB, FeF::HDivKiteDElement, FeT::ScalarElement, Grid, QuadOrd, Jacobi)

Assembles the divergence matrix for a finite element method (FEM) discretization using HDiv-conforming kite elements (`FeF`) and scalar elements (`FeT`) on a given grid.

# Arguments
- `backend`: Computational backend to use (e.g., for parallelization or hardware selection).
- `FTB`: Additional backend or finite element type information (usage context-dependent).
- `FeF::HDivKiteDElement`: The HDiv-conforming finite element for the faces (vector-valued).
- `FeT::ScalarElement`: The scalar finite element for the target space (scalar-valued).
- `Grid`: The computational grid containing geometry and topology information.
- `QuadOrd`: The order of the quadrature rule to use for numerical integration.
- `Jacobi`: Function to compute the Jacobian, its determinant, and pseudo-inverse for mapping reference to physical elements.

# Returns
- `Div`: A sparse matrix representing the divergence operator assembled over all faces of the grid.

# Details
- Computes basis function values and their divergences at quadrature points on both the reference element and its edges.
- Integrates the contributions using quadrature rules for both the element interior and its boundaries.
- Assembles the local divergence matrices into a global sparse matrix using global degree-of-freedom (DoF) mappings.
- Handles orientation and geometry transformations via the provided `Jacobi` function.

# Notes
- The function assumes that the basis functions and their divergences are provided as callable objects in `FeT.phi`, `FeF.phi`, and `FeF.Divphi`.
- The grid and element structures must provide necessary geometric and topological information, including face orientations and global DoF mappings.
"""

function DivMatrix(backend,FTB,FeF::HDivKiteDElement,FeT::ScalarElement,Grid,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(FeT.Type,QuadOrd)
  fFRef  = zeros(FeT.Comp,FeF.DoF,NumQuad)
  fTRef  = zeros(FeT.Comp,FeT.DoF,NumQuad)

  @inbounds for i = 1 : NumQuad
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : FeT.DoF
        fTRef[iComp,iD,i] = FeT.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : FeF.DoF
        fFRef[iComp,iD,i] = FeF.Divphi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  fFRefX  = zeros(FeT.Comp,FeF.DoF,NumQuadL)
  fFRefY  = zeros(FeT.Comp,FeF.DoF,NumQuadL)
  fTRefX  = zeros(FeT.Comp,FeT.DoF,NumQuadL)
  fTRefY  = zeros(FeT.Comp,FeT.DoF,NumQuadL)
  @inbounds for i = 1 : NumQuadL
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : FeT.DoF
        fTRefX[iComp,iD,i] = FeT.phi[iD,iComp](-1.0,PointsL[i])
        fTRefY[iComp,iD,i] = FeT.phi[iD,iComp](PointsL[i],-1.0)
      end
    end
    @inbounds for iD = 1 : FeF.DoF
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
    @inbounds for i = 1 : NumQuad
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      DivLoc += Grid.Faces[iF].Orientation * Weights[i] * (fTRef[:,:,i]' * fFRef[:,:,i])
    end
    @inbounds for i = 1 : NumQuadL
       DivLoc += WeightsL[i] * (fTRefX[:,:,i]' * fFRefX[:,:,i] +  
         fTRefY[:,:,i]' * fFRefY[:,:,i])  
    end   
    @inbounds for j = 1 : size(DivLoc,2)
      @inbounds for i = 1 : size(DivLoc,1)
        push!(RowInd,FeT.Glob[i,iF])
        push!(ColInd,FeF.Glob[j,iF])
        push!(Val,DivLoc[i,j])
      end
    end
  end
  Div = sparse(RowInd, ColInd, Val)
  return Div
end

"""
  DivRhs!(backend, FTB, Div, u, uFeF::HDivKiteDElement, FeT::ScalarElement, Grid, ElemType::Grids.ElementType,
      QuadOrd, Jacobi)

Assembles the right-hand side vector for the divergence operator in a finite element method (FEM) context, specifically for mixed formulations involving H(div) elements on quadrilateral (kite) meshes.

# Arguments
- `backend`: Computational backend to use (e.g., CPU, GPU).
- `FTB`: (Purpose not specified in the code snippet; likely a finite element transformation or buffer).
- `Div`: Output vector to accumulate the divergence contributions (modified in-place).
- `u`: Input vector of degrees of freedom for the H(div) field.
- `uFeF::HDivKiteDElement`: H(div) finite element definition for the face elements.
- `FeT::ScalarElement`: Scalar finite element definition for the test space.
- `Grid`: Grid structure containing mesh and face information.
- `ElemType::Grids.ElementType`: Element type descriptor for quadrature rule selection.
- `QuadOrd`: Quadrature order for numerical integration.
- `Jacobi`: (Purpose not specified in the code snippet; possibly a flag or structure for Jacobian transformations).

# Description
- Precomputes basis function values and their divergences at quadrature points for both the element interior and its edges.
- Loops over all faces in the grid, assembling local contributions to the divergence right-hand side vector.
- Integrates both over the element interior and its boundaries using appropriate quadrature rules.
- Handles orientation of faces and global-to-local degree of freedom mappings.
- Accumulates the computed local divergence contributions into the global `Div` vector.

# Notes
- The function is intended for use in mixed FEM formulations, such as those arising in incompressible flow or electromagnetics.
- Assumes that the finite element types and grid structures provide the necessary basis function evaluation and mapping utilities.
- Uses in-place modification for performance.

# See also
- [`QuadRule`](@ref): For quadrature rule generation.
- [`HDivKiteDElement`](@ref), [`ScalarElement`](@ref): For finite element definitions.
"""

function DivRhs!(backend,FTB,Div,u,uFeF::HDivKiteDElement,FeT::ScalarElement,Grid,ElemType::Grids.ElementType,
  QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  uFRef  = zeros(FeT.Comp,uFeF.DoF,NumQuad)
  fTRef  = zeros(FeT.Comp,FeT.DoF,NumQuad)

  @inbounds for i = 1 : NumQuad
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : FeT.DoF
        fTRef[iComp,iD,i] = FeT.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : uFeF.DoF
        uFRef[iComp,iD,i] = uFeF.Divphi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  uFRefX  = zeros(FeT.Comp,uFeF.DoF,NumQuadL)
  uFRefY  = zeros(FeT.Comp,uFeF.DoF,NumQuadL)
  fTRefX  = zeros(FeT.Comp,FeT.DoF,NumQuadL)
  fTRefY  = zeros(FeT.Comp,FeT.DoF,NumQuadL)
  @inbounds for i = 1 : NumQuadL
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : FeT.DoF
        fTRefX[iComp,iD,i] = FeT.phi[iD,iComp](-1.0,PointsL[i])
        fTRefY[iComp,iD,i] = FeT.phi[iD,iComp](PointsL[i],-1.0)
      end
    end
    @inbounds for iD = 1 : uFeF.DoF
      uFRefX[1,iD,i] = uFeF.phi[iD,1](-1.0,PointsL[i])
      uFRefY[1,iD,i] = uFeF.phi[iD,2](PointsL[i],-1.0)
    end
  end

  DivLoc = zeros(FeT.DoF)
  uF = zeros(uFeF.DoF)

  @inbounds for iF = 1 : Grid.NumFaces
    DivLoc .= 0
    @inbounds for iDoF = 1 : uFeF.DoF
      ind = uFeF.Glob[iDoF,iF]  
      uF[iDoF] = u[ind]  
    end  
    @inbounds for i = 1 : NumQuad
      uLoc = 0.0
      @inbounds for iDoF = 1 : uFeF.DoF
        uLoc += uFRef[1,iDoF,i] * uF[iDoF]  
      end  
      @inbounds for iDoF = 1 : FeT.DoF
        DivLoc[iDoF] -= Grid.Faces[iF].Orientation * Weights[i] * fTRef[1,iDoF,i] * uLoc
      end  
    end
    @inbounds for i = 1 : NumQuadL
       uLocX = 0.0 
       uLocY = 0.0 
       @inbounds for iDoF = 1 : uFeF.DoF
         uLocX += uFRefX[1,iDoF,i] * uF[iDoF]  
         uLocY += uFRefY[1,iDoF,i] * uF[iDoF]  
       end  
      @inbounds for iDoF = 1 : FeT.DoF
         DivLoc[iDoF] -= WeightsL[i] * (fTRefX[1,iDoF,i] * uLocX +
         fTRefY[1,iDoF,i] * uLocY)
      end   
    end   
    @inbounds for iDoF = 1 : FeT.DoF
      ind = FeT.Glob[iDoF,iF]  
      Div[ind] += DivLoc[iDoF]  
    end
  end
end

"""
  CurlMatrix(backend, FTB, FeF::HCurlKiteDElement, FeT::ScalarElement, Grid, QuadOrd, Jacobi)

Assembles the curl matrix for a finite element method (FEM) discretization involving an `HCurlKiteDElement` (vector-valued, H(curl)-conforming element) and a `ScalarElement` (scalar-valued element) over a given grid.

# Arguments
- `backend`: Computational backend to use (e.g., CPU, GPU).
- `FTB`: Finite element type or basis information (purpose depends on context).
- `FeF::HCurlKiteDElement`: The H(curl)-conforming finite element object for the face (vector-valued).
- `FeT::ScalarElement`: The scalar-valued finite element object for the test function.
- `Grid`: The computational grid or mesh structure.
- `QuadOrd`: Quadrature order for numerical integration.
- `Jacobi`: Jacobian or transformation information for mapping reference to physical elements.

# Returns
- `Curl`: A sparse matrix representing the assembled curl operator coupling the given finite element spaces.

# Details
- Computes basis and curl-basis functions at quadrature points on both the reference element and its edges.
- Integrates the contributions using quadrature rules for both the element interior and its boundaries.
- Assembles the local contributions into a global sparse matrix using the global degree-of-freedom (DoF) mappings.
- The resulting matrix can be used in FEM formulations involving curl operators, such as in electromagnetics or fluid dynamics.

# Notes
- Assumes that `FeF` and `FeT` provide callable basis and curl-basis functions, as well as global DoF mappings.
- The function is performance-optimized using `@inbounds` and preallocated arrays.
"""

function CurlMatrix(backend,FTB,FeF::HCurlKiteDElement,FeT::ScalarElement,Grid,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(FeT.Type,QuadOrd)
  fFRef  = zeros(FeF.Comp,FeF.DoF,NumQuad)
  fTRef  = zeros(FeF.Comp,FeT.DoF,NumQuad)

  @inbounds for i = 1 : NumQuad
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : FeT.DoF
        fTRef[iComp,iD,i] = FeT.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : FeF.DoF
        fFRef[iComp,iD,i] = FeF.Curlphi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  fFRefX  = zeros(FeT.Comp,FeF.DoF,NumQuadL)
  fFRefY  = zeros(FeT.Comp,FeF.DoF,NumQuadL)
  fTRefX  = zeros(FeT.Comp,FeT.DoF,NumQuadL)
  fTRefY  = zeros(FeT.Comp,FeT.DoF,NumQuadL)
  @inbounds for i = 1 : NumQuadL
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : FeT.DoF
        fTRefX[iComp,iD,i] = FeT.phi[iD,iComp](-1.0,PointsL[i])
        fTRefY[iComp,iD,i] = FeT.phi[iD,iComp](PointsL[i],-1.0)
      end
    end
    @inbounds for iD = 1 : FeF.DoF
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
    @inbounds for i = 1 : NumQuad
      CurlLoc += Weights[i] * (fTRef[:,:,i]' * fFRef[:,:,i])
    end
    @inbounds for i = 1 : NumQuadL
      CurlLoc += WeightsL[i] * (fTRefX[:,:,i]' * fFRefX[:,:,i] +  
        fTRefY[:,:,i]' * fFRefY[:,:,i])  
    end   
    @inbounds for j = 1 : size(CurlLoc,2)
      @inbounds for i = 1 : size(CurlLoc,1)
        push!(RowInd,FeT.Glob[i,iF])
        push!(ColInd,FeF.Glob[j,iF])
        push!(Val,CurlLoc[i,j])
      end
    end
  end
  Curl = sparse(RowInd, ColInd, Val)
  return Curl
end

"""
  GradMatrix(backend, FTB, FeF::ScalarElement, FeT::HDivConfElement, Grid, QuadOrd, Jacobi)

Assembles the gradient matrix for a finite element method (FEM) problem involving a scalar element (`FeF`) and an H(div)-conforming element (`FeT`) on a given grid.

# Arguments
- `backend`: Computational backend to use (e.g., CPU, GPU).
- `FTB`: Finite element type or basis (usage context-dependent).
- `FeF::ScalarElement`: Scalar finite element object, providing basis functions and gradients.
- `FeT::HDivConfElement`: H(div)-conforming finite element object, providing vector-valued basis functions.
- `Grid`: Grid or mesh structure containing geometric and topological information.
- `QuadOrd`: Quadrature order for numerical integration.
- `Jacobi`: Function to compute the Jacobian, its determinant, and pseudo-inverse for mapping reference to physical elements.

# Returns
- `Grad`: Sparse matrix representing the assembled gradient operator between the given finite element spaces.

# Description
The function computes the local contributions to the gradient matrix for each face in the grid using numerical quadrature. It evaluates the basis functions and their gradients at quadrature points, applies the appropriate geometric transformations, and assembles the global sparse matrix.

# Notes
- The function assumes that `FeF` and `FeT` provide callable basis and gradient functions.
- The grid structure must provide face orientation and global degree-of-freedom (DoF) mappings.
- The assembly is performed in a memory-efficient way using sparse matrix construction.

"""

function GradMatrix(backend,FTB,FeF::ScalarElement,FeT::HDivConfElement,Grid,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(Grid.Type,QuadOrd)
  fFRef  = zeros(FeT.Comp,FeF.DoF,NumQuad)
  fTRef  = zeros(FeT.Comp,FeT.DoF,NumQuad)

  @inbounds for i = 1 : NumQuad
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : FeT.DoF
        fTRef[iComp,iD,i] = FeT.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : FeF.DoF
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
    @inbounds for i = 1 : NumQuad
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      GradLoc += Grid.Faces[iF].Orientation * Weights[i] * (fTRef[:,:,i]' * fFRef[:,:,i])
    end
    @inbounds for j = 1 : size(GradLoc,2)
      @inbounds for i = 1 : size(GradLoc,1)
        push!(RowInd,FeT.Glob[i,iF])
        push!(ColInd,FeF.Glob[j,iF])
        push!(Val,GradLoc[i,j])
      end
    end
  end
  Grad = sparse(RowInd, ColInd, Val)
  return Grad
end
    
"""
  GradMatrix(backend, FTB, FeF::ScalarElement, FeT::HDivKiteDElement, Grid, QuadOrd, Jacobi)

Assembles the gradient matrix for a finite element method (FEM) discretization involving a scalar element (`FeF`) and an H(div) kite element (`FeT`) on a given grid.

# Arguments
- `backend`: Computational backend to use for matrix assembly.
- `FTB`: Additional backend or transformation information (purpose depends on context).
- `FeF::ScalarElement`: Scalar finite element providing basis functions and gradients.
- `FeT::HDivKiteDElement`: H(div) finite element providing vector-valued basis functions.
- `Grid`: Grid structure containing mesh and face information.
- `QuadOrd`: Quadrature order for numerical integration.
- `Jacobi`: Function to compute the Jacobian, its determinant, and pseudo-inverse for mapping reference to physical elements.

# Returns
- `Grad`: Sparse matrix representing the assembled gradient operator.

# Description
The function computes the local contributions to the gradient matrix by evaluating basis functions and their gradients at quadrature points in the reference element and on its edges. It then maps these contributions to the global matrix using the grid's connectivity and orientation information. The resulting sparse matrix can be used in FEM solvers for problems involving mixed formulations or H(div) spaces.

# Notes
- The function assumes that `FeF` and `FeT` provide callable basis and gradient functions indexed by degree of freedom and component.
- The assembly uses in-place operations and preallocated arrays for efficiency.
- The function is specialized for 2D elements with two components and their corresponding edge integrals.
"""

function GradMatrix(backend,FTB,FeF::ScalarElement,FeT::HDivKiteDElement,Grid,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(Grid.Type,QuadOrd)
  fFRef  = zeros(FeT.Comp,FeF.DoF,NumQuad)
  fTRef  = zeros(FeT.Comp,FeT.DoF,NumQuad)

  @inbounds for i = 1 : NumQuad
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : FeT.DoF
        fTRef[iComp,iD,i] = FeT.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : FeF.DoF
        fFRef[iComp,iD,i] = FeF.Gradphi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  fTRefX  = zeros(FeF.Comp,FeT.DoF,NumQuadL)
  fTRefY  = zeros(FeF.Comp,FeT.DoF,NumQuadL)
  fFRefX  = zeros(FeF.Comp,FeF.DoF,NumQuadL)
  fFRefY  = zeros(FeF.Comp,FeF.DoF,NumQuadL)
  @inbounds for i = 1 : NumQuadL
    @inbounds for iComp = 1 : FeF.Comp
      @inbounds for iD = 1 : FeF.DoF
        fFRefX[iComp,iD,i] = FeF.phi[iD,iComp](1.0,PointsL[i])
        fFRefY[iComp,iD,i] = FeF.phi[iD,iComp](PointsL[i],1.0)
      end
    end
    @inbounds for iD = 1 : FeT.DoF
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
    @inbounds for i = 1 : NumQuad
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      GradLoc += Grid.Faces[iF].Orientation * Weights[i] * (fTRef[:,:,i]' * fFRef[:,:,i])
    end
    @inbounds for i = 1 : NumQuadL
       GradLoc -= WeightsL[i] * (fTRefX[:,:,i]' * fFRefX[:,:,i] +  
         fTRefY[:,:,i]' * fFRefY[:,:,i])  
    end   
    @inbounds for j = 1 : size(GradLoc,2)
      @inbounds for i = 1 : size(GradLoc,1)
        push!(RowInd,FeT.Glob[i,iF])
        push!(ColInd,FeF.Glob[j,iF])
        push!(Val,GradLoc[i,j])
      end
    end
  end
  Grad = sparse(RowInd, ColInd, Val)
  return Grad
end

"""
  GradHeightSquared!(backend, FTB, Rhs, FeT::HDivConfElement, h, hFeF::ScalarElement,
            Grid, ElemType::Grids.ElementType, QuadOrd, Jacobi)

Assembles the contribution of the gradient of the squared height field to the right-hand side vector `Rhs`
for a finite element method (FEM) discretization, typically used in shallow water or similar PDEs.

# Arguments
- `backend`: Computational backend or context (not used directly in this function).
- `FTB`: Additional backend or context parameter (not used directly in this function).
- `Rhs`: The right-hand side vector to be updated in-place.
- `FeT::HDivConfElement`: The finite element configuration for the vector-valued (H(div)) space.
- `h`: The global vector of height field values.
- `hFeF::ScalarElement`: The finite element configuration for the scalar (height) field.
- `Grid`: The grid or mesh structure containing face and connectivity information.
- `ElemType::Grids.ElementType`: The element type used for quadrature rule selection.
- `QuadOrd`: The order of the quadrature rule.
- `Jacobi`: The Jacobian or transformation information (not used directly in this function).

# Description
For each face in the grid, the function:
1. Interpolates the height field `h` at quadrature points using the scalar basis functions.
2. Computes the divergence of the vector basis functions at quadrature points.
3. Assembles the local contributions of the gradient of the squared height (`∇(h^2)`) to the global right-hand side vector `Rhs`, scaled by gravity.

# Notes
- The function assumes that the basis functions and their divergences are provided as callable objects.
- The gravity constant `Grav` is hardcoded as `9.80616`.
- The function operates in-place and is performance-optimized with `@inbounds` and preallocated arrays.

# Modifies
- `Rhs`: The right-hand side vector is updated in-place with the assembled contributions.

# See also
- `QuadRule`
- `HDivConfElement`
- `ScalarElement`
"""

function GradHeightSquared!(backend,FTB,Rhs,FeT::HDivConfElement,h,hFeF::ScalarElement,
  Grid,ElemType::Grids.ElementType,QuadOrd,Jacobi)

  Grav = 9.80616
  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  hFRef  = zeros(hFeF.Comp,hFeF.DoF,NumQuad)
  DivfTRef  = zeros(hFeF.Comp,FeT.DoF,NumQuad)


  @inbounds for iQ = 1 : NumQuad
    @inbounds for iComp = 1 : hFeF.Comp
      @inbounds for iD = 1 : FeT.DoF
        DivfTRef[iComp,iD,iQ] = FeT.Divphi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    @inbounds for iComp = 1 : hFeF.Comp
      @inbounds for iD = 1 : hFeF.DoF
        hFRef[iComp,iD,iQ] = hFeF.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
  end
  GradLoc = zeros(FeT.DoF)
  hLoc = zeros(hFeF.DoF)

  @inbounds for iF = 1 : Grid.NumFaces
    GradLoc .= 0
    @inbounds for iDoF = 1 : hFeF.DoF
      ind = hFeF.Glob[iDoF,iF]  
      hLoc[iDoF] = h[ind]
    end  
    @inbounds for iQ = 1 : NumQuad
      hhLoc = 0.0
      @inbounds for iDoF = 1 : hFeF.DoF
        hhLoc += hFRef[1,iDoF,iQ] * hLoc[iDoF]
      end  
      @inbounds for iDoF = 1 : FeT.DoF
        GradLoc[iDoF] += Grid.Faces[iF].Orientation * Weights[iQ] * DivfTRef[1,iDoF,iQ] * hhLoc^2
      end  
    end
    @inbounds for iDoF = 1 : FeT.DoF
      ind = FeT.Glob[iDoF,iF]  
      Rhs[ind] += 0.5 * Grav * GradLoc[iDoF] 
    end  
  end
end

"""
  GradKinHeight!(backend, FTB, Rhs, h, hFeF::ScalarElement, u, uFeF::HDivConfElement,
           FeT::HDivConfElement, Grid, ElemType::Grids.ElementType, QuadOrd, Jacobi)

Assembles the gradient of the kinetic height term for a finite element method (FEM) discretization,
and adds the local contributions to the global right-hand side vector `Rhs`.

# Arguments
- `backend`: Computational backend (not used directly in this function).
- `FTB`: (Purpose not specified in the code snippet).
- `Rhs`: Global right-hand side vector to be updated in-place.
- `h`: Vector of height values at global degrees of freedom.
- `hFeF::ScalarElement`: Scalar finite element descriptor for the height field.
- `u`: Vector of velocity values at global degrees of freedom.
- `uFeF::HDivConfElement`: H(div)-conforming finite element descriptor for the velocity field.
- `FeT::HDivConfElement`: H(div)-conforming finite element descriptor for the test functions.
- `Grid`: Grid structure containing mesh and face information.
- `ElemType::Grids.ElementType`: Element type for quadrature rule selection.
- `QuadOrd`: Quadrature order for numerical integration.
- `Jacobi`: Function to compute the Jacobian, its determinant, and pseudo-inverse for mapping reference to physical elements.

# Description
For each face in the grid, this function:
1. Evaluates basis functions and their divergences at quadrature points.
2. Maps local degrees of freedom to global indices and gathers local solution values.
3. Computes the local kinetic energy and potential energy contributions at quadrature points.
4. Integrates these contributions using the quadrature rule and assembles them into the global right-hand side vector `Rhs`.

The gravitational constant is hardcoded as `Grav = 9.80616`.

# Notes
- The function operates in-place on `Rhs`.
- Assumes that the finite element descriptors (`hFeF`, `uFeF`, `FeT`) and the grid structure are properly initialized and compatible.
- The function is performance-oriented, using `@inbounds` to skip bounds checking in loops.

# Returns
Nothing. The function updates `Rhs` in-place.
"""

function GradKinHeight!(backend,FTB,Rhs,h,hFeF::ScalarElement,u,uFeF::HDivConfElement,
  FeT::HDivConfElement,Grid,ElemType::Grids.ElementType,QuadOrd,Jacobi)

  Grav =  9.80616
  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  hFRef  = zeros(hFeF.Comp,hFeF.DoF,NumQuad)
  uFRef  = zeros(FeT.Comp,FeT.DoF,NumQuad)
  DivfTRef  = zeros(hFeF.Comp,FeT.DoF,NumQuad)


  @inbounds for iQ = 1 : NumQuad
    @inbounds for iComp = 1 : hFeF.Comp
      @inbounds for iD = 1 : FeT.DoF
        DivfTRef[iComp,iD,iQ] = FeT.Divphi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    @inbounds for iComp = 1 : hFeF.Comp
      @inbounds for iD = 1 : hFeF.DoF
        hFRef[iComp,iD,iQ] = hFeF.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    @inbounds for iComp = 1 : uFeF.Comp
      @inbounds for iD = 1 : uFeF.DoF
        uFRef[iComp,iD,iQ] = uFeF.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
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
    @inbounds for iDoF = 1 : uFeF.DoF
      ind = uFeF.Glob[iDoF,iF]  
      uLoc[iDoF] = u[ind]
    end  
    @inbounds for iDoF = 1 : hFeF.DoF
      ind = hFeF.Glob[iDoF,iF]  
      hLoc[iDoF] = Grav * h[ind]
    end  
    @inbounds for iQ = 1 : NumQuad
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,Points[iQ,1],Points[iQ,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      u1 = 0.0
      u2 = 0.0
      @inbounds for iDoF = 1 : uFeF.DoF
        u1 += uFRef[1,iDoF,iQ] * uLoc[iDoF]
        u2 += uFRef[2,iDoF,iQ] * uLoc[iDoF]
      end  
      uLoc1 = 1 / detDFLoc * (DF[1,1] * u1 + DF[1,2] * u2)
      uLoc2 = 1 / detDFLoc * (DF[2,1] * u1 + DF[2,2] * u2)
      uLoc3 = 1 / detDFLoc * (DF[3,1] * u1 + DF[3,2] * u2)
      hhLoc = 0.5 * (uLoc1 * uLoc1 + uLoc2 * uLoc2 + uLoc3 * uLoc3)
      @inbounds for iDoF = 1 : hFeF.DoF
        hhLoc += hFRef[1,iDoF,iQ] * hLoc[iDoF] 
      end  
      @inbounds for iDoF = 1 : FeT.DoF
        GradLoc[iDoF] += Grid.Faces[iF].Orientation * Weights[iQ] * DivfTRef[1,iDoF,iQ] * hhLoc
      end  
    end
    @inbounds for iDoF = 1 : FeT.DoF
      ind = FeT.Glob[iDoF,iF]  
      Rhs[ind] += GradLoc[iDoF]
    end  
  end
end

"""
  GradKinHeightInter!(backend, FTB, Rhs, h, hFeF::ScalarElement, u, uFeF::HDivConfElement,
             FeT::HDivConfElement, Grid, ElemType::Grids.ElementType, QuadOrd, Jacobi)

Assembles the gradient of the kinetic height interaction term for a finite element method (FEM) discretization, and adds its contribution to the right-hand side vector `Rhs`.

# Arguments
- `backend`: Computational backend to use (e.g., CPU, GPU).
- `FTB`: (Purpose not specified in the code snippet; likely a data structure related to the finite element assembly).
- `Rhs`: Right-hand side vector to be updated in-place.
- `h`: Vector of height (scalar field) degrees of freedom.
- `hFeF::ScalarElement`: Scalar finite element descriptor for the height field.
- `u`: Vector of velocity (vector field) degrees of freedom.
- `uFeF::HDivConfElement`: H(div)-conforming finite element descriptor for the velocity field.
- `FeT::HDivConfElement`: H(div)-conforming finite element descriptor for the test functions.
- `Grid`: Grid or mesh data structure containing face and element information.
- `ElemType::Grids.ElementType`: Type of the finite element (e.g., triangle, quadrilateral).
- `QuadOrd`: Quadrature order for numerical integration.
- `Jacobi`: Function to compute the Jacobian, its determinant, and pseudo-inverse for the mapping from reference to physical element.

# Description
This function computes the local contributions of the gradient of the kinetic height interaction term for each face in the grid, using numerical quadrature. It evaluates basis functions and their divergences at quadrature points, transforms velocity components, and assembles the contributions into the global right-hand side vector `Rhs`. The gravitational constant is hardcoded as `Grav = 9.80616`.

# Notes
- The function operates in-place and assumes all input arrays are properly sized and initialized.
- The function is performance-optimized using `@inbounds` to skip array bounds checking.
- The function is intended for use in the assembly phase of a mixed finite element method for shallow water or similar PDEs.

# See also
- [`ScalarElement`](@ref)
- [`HDivConfElement`](@ref)
- [`QuadRule`](@ref)
"""

function GradKinHeightInter!(backend,FTB,Rhs,h,hFeF::ScalarElement,u,uFeF::HDivConfElement,
  FeT::HDivConfElement,Grid,ElemType::Grids.ElementType,QuadOrd,Jacobi)

  Grav =  9.80616
  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  hFRef  = zeros(hFeF.Comp,hFeF.DoF,NumQuad)
  PointsI = hFeF.points
  NumPointsI = size(PointsI,1)
  uFRef  = zeros(FeT.Comp,FeT.DoF,NumPointsI)
  DivfTRef  = zeros(hFeF.Comp,FeT.DoF,NumQuad)


  @inbounds for iQ = 1 : NumQuad
    @inbounds for iComp = 1 : hFeF.Comp
      @inbounds for iD = 1 : FeT.DoF
        DivfTRef[iComp,iD,iQ] = FeT.Divphi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    @inbounds for iComp = 1 : hFeF.Comp
      @inbounds for iD = 1 : hFeF.DoF
        hFRef[iComp,iD,iQ] = hFeF.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
  end
  @inbounds for iP = 1 : NumPointsI
    @inbounds for iComp = 1 : uFeF.Comp
      @inbounds for iD = 1 : uFeF.DoF
        uFRef[iComp,iD,iP] = uFeF.phi[iD,iComp](PointsI[iP,1],PointsI[iP,2])
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
    @inbounds for iDoF = 1 : uFeF.DoF
      ind = uFeF.Glob[iDoF,iF]  
      uLoc[iDoF] = u[ind]
    end  
    @inbounds for iDoF = 1 : hFeF.DoF
      ind = hFeF.Glob[iDoF,iF]  
      hLoc[iDoF] = Grav * h[ind]
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,PointsI[iDoF,1],PointsI[iDoF,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      u1 = 0.0
      u2 = 0.0
      @inbounds for jDoF = 1 : uFeF.DoF
        u1 += uFRef[1,jDoF,iDoF] * uLoc[jDoF]
        u2 += uFRef[2,jDoF,iDoF] * uLoc[jDoF]
      end  
      uLoc1 = 1 / detDFLoc * (DF[1,1] * u1 + DF[1,2] * u2)
      uLoc2 = 1 / detDFLoc * (DF[2,1] * u1 + DF[2,2] * u2)
      uLoc3 = 1 / detDFLoc * (DF[3,1] * u1 + DF[3,2] * u2)
      hLoc[iDoF] += 0.5 * (uLoc1 * uLoc1 + uLoc2 * uLoc2 + uLoc3 * uLoc3)
    end  
    @inbounds for iQ = 1 : NumQuad
      hhLoc = 0.0
      @inbounds for iDoF = 1 : hFeF.DoF
        hhLoc += hFRef[1,iDoF,iQ] * hLoc[iDoF]
      end  
      @inbounds for iDoF = 1 : FeT.DoF
        GradLoc[iDoF] += Grid.Faces[iF].Orientation * Weights[iQ] * DivfTRef[1,iDoF,iQ] * hhLoc
      end  
    end
    @inbounds for iDoF = 1 : FeT.DoF
      ind = FeT.Glob[iDoF,iF]  
      Rhs[ind] += GradLoc[iDoF]
    end  
  end
end

"""
  GradKinHeight!(backend, FTB, Rhs, h, hFeF::ScalarElement, u, uFeF::HDivKiteDElement,
           FeT::HDivKiteDElement, Grid, ElemType::Grids.ElementType, QuadOrd, Jacobi)

Assembles the gradient of the kinetic height term for a finite element method (FEM) discretization,
and adds its contribution to the right-hand side vector `Rhs`. This function is typically used in the
context of shallow water or similar equations, where the kinetic and potential energy terms are
projected onto a finite element basis.

# Arguments
- `backend`: Computational backend (not used directly in this function, but may be required for compatibility).
- `FTB`: Placeholder or auxiliary argument (usage context-dependent).
- `Rhs`: The right-hand side vector to which the local contributions will be added (modified in-place).
- `h`: Vector of height (or potential) degrees of freedom.
- `hFeF::ScalarElement`: Scalar finite element descriptor for the height field.
- `u`: Vector of velocity degrees of freedom.
- `uFeF::HDivKiteDElement`: HDiv-conforming finite element descriptor for the velocity field.
- `FeT::HDivKiteDElement`: HDiv-conforming finite element descriptor for the test functions.
- `Grid`: Grid or mesh data structure containing face and element information.
- `ElemType::Grids.ElementType`: Type of the finite element (e.g., quadrilateral, triangle).
- `QuadOrd`: Quadrature order for numerical integration.
- `Jacobi`: Function to compute the Jacobian, its determinant, and pseudo-inverse for mapping reference to physical elements.

# Details
- Computes local projections of the kinetic and potential energy terms at quadrature points.
- Integrates these projections over each face of the grid using the specified quadrature rule.
- Handles both volume and edge (face) integrals, accounting for orientation and mapping via the Jacobian.
- The gravitational constant is hardcoded as `Grav = 9.80616`.
- The function is performance-optimized using `@inbounds` and preallocated arrays.

# Modifies
- `Rhs`: The right-hand side vector is updated in-place with the assembled local contributions.

# Notes
- Assumes that the finite element descriptors (`hFeF`, `uFeF`, `FeT`) and the grid structure provide all necessary shape function evaluations and global indexing.
- The function is intended for internal use in FEM assembly routines and may require adaptation for other element types or problem settings.

# See Also
- `QuadRule`
- `Jacobi`
- `ScalarElement`
- `HDivKiteDElement`
"""

#! heißt drin berechnet
function GradKinHeight!(backend,FTB,Rhs,h,hFeF::ScalarElement,u,uFeF::HDivKiteDElement,
  FeT::HDivKiteDElement,Grid,ElemType::Grids.ElementType,QuadOrd,Jacobi)

  Grav =  9.80616
  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  hFRef  = zeros(hFeF.Comp,hFeF.DoF,NumQuad)
  fTRef  = zeros(hFeF.Comp,FeT.DoF,NumQuad)
  uFRef  = zeros(uFeF.Comp,uFeF.DoF,NumQuad)

  @inbounds for iQ = 1 : NumQuad
    @inbounds for iComp = 1 : hFeF.Comp
      @inbounds for iDoF = 1 : FeT.DoF
        fTRef[iComp,iDoF,iQ] = FeT.Divphi[iDoF,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    @inbounds for iComp = 1 : hFeF.Comp
      @inbounds for iDoF = 1 : hFeF.DoF
        hFRef[iComp,iDoF,iQ] = hFeF.phi[iDoF,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    @inbounds for iComp = 1 : uFeF.Comp
      @inbounds for iDoF = 1 : uFeF.DoF
        uFRef[iComp,iDoF,iQ] = uFeF.phi[iDoF,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
  end
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  fTRefX  = zeros(hFeF.Comp,FeT.DoF,NumQuadL)
  fTRefY  = zeros(hFeF.Comp,FeT.DoF,NumQuadL)
  hFRefX  = zeros(hFeF.Comp,hFeF.DoF,NumQuadL)
  hFRefY  = zeros(hFeF.Comp,hFeF.DoF,NumQuadL)
  uFRefX  = zeros(uFeF.Comp,uFeF.DoF,NumQuadL)
  uFRefY  = zeros(uFeF.Comp,uFeF.DoF,NumQuadL)
  @inbounds for iQ = 1 : NumQuadL
    @inbounds for iComp = 1 : hFeF.Comp
      @inbounds for iDoF = 1 : hFeF.DoF
        hFRefX[iComp,iDoF,iQ] = hFeF.phi[iDoF,iComp](-1.0,PointsL[iQ])
        hFRefY[iComp,iDoF,iQ] = hFeF.phi[iDoF,iComp](PointsL[iQ],-1.0)
      end
    end
    @inbounds for iDoF = 1 : uFeF.DoF
      uFRefX[1,iDoF,iQ] = uFeF.phi[iDoF,1](-1.0,PointsL[iQ])
      uFRefX[2,iDoF,iQ] = uFeF.phi[iDoF,2](-1.0,PointsL[iQ])
      uFRefY[1,iDoF,iQ] = uFeF.phi[iDoF,1](PointsL[iQ],-1.0)
      uFRefY[2,iDoF,iQ] = uFeF.phi[iDoF,2](PointsL[iQ],-1.0)
    end
    @inbounds for iDoF = 1 : FeT.DoF
      fTRefX[1,iDoF,iQ] = FeT.phi[iDoF,1](-1.0,PointsL[iQ])
      fTRefY[1,iDoF,iQ] = FeT.phi[iDoF,2](PointsL[iQ],-1.0)
    end
  end
  GradLoc = zeros(FeT.DoF)
  DF = zeros(3,2) * Grav
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  uF = zeros(uFeF.DoF)
  hF = zeros(hFeF.DoF)

  @inbounds for iF = 1 : Grid.NumFaces
    GradLoc .= 0
    @inbounds for iDoF = 1 : uFeF.DoF
      ind = uFeF.Glob[iDoF,iF]  
      uF[iDoF] = u[ind]
    end  
    @inbounds for iDoF = 1 : hFeF.DoF
      ind = hFeF.Glob[iDoF,iF]  
      hF[iDoF] = h[ind]
    end  
    @inbounds for iQ = 1 : NumQuad
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,Points[iQ,1],Points[iQ,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      u1 = 0.0
      u2 = 0.0
      @inbounds for iDoF = 1 : uFeF.DoF
        u1 += uFRef[1,iDoF,iQ] * uF[iDoF]
        u2 += uFRef[2,iDoF,iQ] * uF[iDoF]
      end
      uLoc1 = 1 / detDFLoc * (DF[1,1] * u1 + DF[1,2] * u2)
      uLoc2 = 1 / detDFLoc * (DF[2,1] * u1 + DF[2,2] * u2)
      uLoc3 = 1 / detDFLoc * (DF[3,1] * u1 + DF[3,2] * u2)
      hLoc = 0.5 * (uLoc1 * uLoc1 + uLoc2 * uLoc2 + uLoc3 * uLoc3)
      @inbounds for iDoF = 1 : hFeF.DoF 
        hLoc += hFRef[1,iDoF,iQ] * Grav * hF[iDoF]
      end  
      @inbounds for iDoF = 1 : FeT.DoF
        GradLoc[iDoF] += Grid.Faces[iF].Orientation * Weights[iQ] * fTRef[1,iDoF,iQ] * hLoc
      end
    end
    @inbounds for iQ = 1 : NumQuadL
      hLocX = 0.0 
      hLocY = 0.0 
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,-1.0,PointsL[iQ],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      u1X = 0.0
      u2X = 0.0
      @inbounds for iDoF = 1 : uFeF.DoF
        u1X += uFRefX[1,iDoF,iQ] * uF[iDoF]
        u2X += uFRefX[2,iDoF,iQ] * uF[iDoF]
      end
      uLoc1X = 1 / detDFLoc * (DF[1,1] * u1X + DF[1,2] * u2X)
      uLoc2X = 1 / detDFLoc * (DF[2,1] * u1X + DF[2,2] * u2X)
      uLoc3X = 1 / detDFLoc * (DF[3,1] * u1X + DF[3,2] * u2X)
      hLocX = 0.5 * (uLoc1X * uLoc1X + uLoc2X * uLoc2X + uLoc3X * uLoc3X)
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,PointsL[iQ],-1.0,Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      u1Y = 0.0
      u2Y = 0.0
      @inbounds for iDoF = 1 : uFeF.DoF
        u1Y += uFRefY[1,iDoF,iQ] * uF[iDoF]
        u2Y += uFRefY[2,iDoF,iQ] * uF[iDoF]
      end
      uLoc1Y = 1 / detDFLoc * (DF[1,1] * u1Y + DF[1,2] * u2Y)
      uLoc2Y = 1 / detDFLoc * (DF[2,1] * u1Y + DF[2,2] * u2Y)
      uLoc3Y = 1 / detDFLoc * (DF[3,1] * u1Y + DF[3,2] * u2Y)
      hLocY = 0.5 * (uLoc1Y * uLoc1Y + uLoc2Y * uLoc2Y + uLoc3Y * uLoc3Y)
      @inbounds for iDoF = 1 : hFeF.DoF
        hLocX += hFRefX[1,iDoF,iQ] * Grav * hF[iDoF]
        hLocY += hFRefY[1,iDoF,iQ] * Grav * hF[iDoF]
      end
      @inbounds for iDoF = 1 : FeT.DoF
        GradLoc[iDoF] += WeightsL[iQ] * (fTRefX[1,iDoF,iQ] * hLocX + fTRefY[1,iDoF,iQ] * hLocY)
      end
    end   
    @inbounds for iDoF = 1 : FeT.DoF
      ind = FeT.Glob[iDoF,iF]  
      Rhs[ind] += GradLoc[iDoF]
    end  
  end
end

"""
  DivRhs!(backend, FTB, Div, u, uFeF::HDivKiteDElement, h, hFeF::ScalarElement, FeT::ScalarElement, Grid,
      ElemType::Grids.ElementType, QuadOrd, Jacobi)

Assembles the right-hand side (RHS) vector for the divergence operator in a finite element method (FEM) context, specifically for mixed finite element spaces involving H(div) and scalar elements on a kite-shaped reference element.

# Arguments
- `backend`: Computational backend to use (e.g., CPU, GPU).
- `FTB`: (Purpose not specified; likely a finite element basis or transformation object.)
- `Div`: Output vector to accumulate the divergence RHS contributions (modified in-place).
- `u`: Vector of degrees of freedom (DoFs) for the H(div) field.
- `uFeF::HDivKiteDElement`: H(div) finite element object for the vector field.
- `h`: Vector of DoFs for the scalar field.
- `hFeF::ScalarElement`: Scalar finite element object for the scalar field.
- `FeT::ScalarElement`: Scalar finite element object for the test space.
- `Grid`: Grid or mesh structure containing element and face information.
- `ElemType::Grids.ElementType`: Type of the element (e.g., reference kite).
- `QuadOrd`: Quadrature order for numerical integration.
- `Jacobi`: (Purpose not specified; possibly a flag or structure for Jacobi preconditioning or integration.)

# Description
- Computes local contributions to the divergence RHS by evaluating basis functions and their divergences at quadrature points.
- Integrates over both the element interior and its faces/edges using appropriate quadrature rules.
- Accumulates the results into the global RHS vector `Div` using the global DoF mapping.

# Notes
- The function is intended for use in mixed FEM formulations, such as those arising in fluid dynamics or electromagnetics.
- Assumes that the basis functions and their divergences are provided as callable objects in the finite element structures.
- Uses in-place modification for performance.

# In-place
- Modifies the `Div` vector in-place.

# Performance
- Uses `@inbounds` and preallocated arrays for efficiency.
- Suitable for use in performance-critical assembly loops.

# See also
- `QuadRule`
- `HDivKiteDElement`
- `ScalarElement`
"""

function DivRhs!(backend,FTB,Div,u,uFeF::HDivKiteDElement,h,hFeF::ScalarElement,FeT::ScalarElement,Grid,
  ElemType::Grids.ElementType,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  uFFRef  = zeros(FeT.Comp,uFeF.DoF,NumQuad)
  hFFRef  = zeros(FeT.Comp,hFeF.DoF,NumQuad)
  fTRef  = zeros(FeT.Comp,FeT.DoF,NumQuad)

  @inbounds for iQ = 1 : NumQuad
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : FeT.DoF
        fTRef[iComp,iD,iQ] = FeT.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : uFeF.DoF
        uFFRef[iComp,iD,iQ] = uFeF.Divphi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : hFeF.DoF
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
  @inbounds for iQ = 1 : NumQuadL
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : FeT.DoF
        fTRefX[iComp,iD,iQ] = FeT.phi[iD,iComp](1.0,PointsL[iQ])
        fTRefY[iComp,iD,iQ] = FeT.phi[iD,iComp](PointsL[iQ],1.0)
      end
    end
    @inbounds for iD = 1 : uFeF.DoF
      uFFRefX[1,iD,iQ] = uFeF.phi[iD,1](1.0,PointsL[iQ])
      uFFRefY[1,iD,iQ] = uFeF.phi[iD,2](PointsL[iQ],1.0)
    end
    @inbounds for iD = 1 : hFeF.DoF
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

  @inbounds for iF = 1 : Grid.NumFaces
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
    @inbounds for iQ = 1 : NumQuad
      @inbounds for iDoFFeT = 1 : FeT.DoF
        DivLoc[iDoFFeT] +=  Weights[iQ] * fTRef[1,iDoFFeT,iQ] * uFLoc[iQ] * hFLoc[iQ]
      end
    end
    @inbounds for iQ = 1 : NumQuadL
      @inbounds for iDoFFeT = 1 : FeT.DoF
        DivLoc[iDoFFeT] +=  WeightsL[iQ] * (fTRefX[1,iDoFFeT,iQ] * uFLocX[iQ] * hFLocX[iQ] +
          fTRefY[1,iDoFFeT,iQ] * uFLocY[iQ] * hFLocY[iQ])
      end
    end   
    @inbounds for iDoFFeT = 1 : FeT.DoF
      ind = FeT.Glob[iDoFFeT,iF]
      Div[ind] += DivLoc[iDoFFeT]
    end
  end
end

"""
  GradRhs!(backend, FTB, Grad, h, hFeF::ScalarElement, u, uFeF::HDivKiteDElement, FeT::HDivKiteDElement, Grid,
       ElemType::Grids.ElementType, QuadOrd, Jacobi)

Assembles the gradient right-hand side (RHS) vector for a finite element method (FEM) problem involving scalar and vector fields on a grid with kite-shaped elements.

# Arguments
- `backend`: Computational backend to use (e.g., CPU, GPU).
- `FTB`: (Purpose not specified in code snippet; likely a context or buffer object).
- `Grad`: Array to accumulate the assembled gradient RHS.
- `h`: Array of scalar field values (e.g., height or pressure) at global degrees of freedom.
- `hFeF::ScalarElement`: Scalar finite element function space for `h`.
- `u`: Array of vector field values (e.g., velocity) at global degrees of freedom.
- `uFeF::HDivKiteDElement`: H(div) finite element function space for `u` on kite-shaped elements.
- `FeT::HDivKiteDElement`: H(div) finite element function space for test functions on kite-shaped elements.
- `Grid`: Grid structure containing mesh and face information.
- `ElemType::Grids.ElementType`: Type of element used in the grid.
- `QuadOrd`: Quadrature order for numerical integration.
- `Jacobi`: Function to compute the Jacobian, its determinant, and pseudo-inverse for coordinate transformations.

# Description
This function computes the local contributions to the gradient RHS vector for each face in the grid, using quadrature rules for both the element interior and its boundaries. It evaluates basis functions and their divergences at quadrature points, transforms field values to the physical domain, and assembles the contributions into the global RHS vector.

The function supports mixed finite element spaces (scalar and H(div) vector fields) and handles the mapping between reference and physical elements. It is designed for use in discontinuous Galerkin or mixed FEM formulations, particularly for problems involving conservation laws or fluid dynamics.

# Notes
- The function modifies `Grad` in-place.
- Assumes that the finite element spaces and grid structures are properly initialized and compatible.
- The function is performance-oriented, using `@inbounds` and preallocated arrays to minimize allocations and bounds checking.

# See Also
- `QuadRule`
- `ScalarElement`
- `HDivKiteDElement`
- `Grids.ElementType`
"""

function GradRhs!(backend,FTB,Grad,h,hFeF::ScalarElement,u,uFeF::HDivKiteDElement,FeT::HDivKiteDElement,Grid,
  ElemType::Grids.ElementType,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  uFFRef  = zeros(uFeF.Comp,uFeF.DoF,NumQuad)
  hFFRef  = zeros(hFeF.Comp,hFeF.DoF,NumQuad)
  fTRef  = zeros(hFeF.Comp,FeT.DoF,NumQuad)

  @inbounds for iQ = 1 : NumQuad
    @inbounds for iComp = 1 : hFeF.Comp
      @inbounds for iD = 1 : FeT.DoF
        fTRef[iComp,iD,iQ] = FeT.Divphi[iD,iComp](Points[iQ,1],Points[iQ,2])
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
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  uFFRefX  = zeros(uFeF.Comp,uFeF.DoF,NumQuadL)
  uFFRefY  = zeros(uFeF.Comp,uFeF.DoF,NumQuadL)
  hFFRefX  = zeros(hFeF.Comp,hFeF.DoF,NumQuadL)
  hFFRefY  = zeros(hFeF.Comp,hFeF.DoF,NumQuadL)
  fTRefX  = zeros(hFeF.Comp,FeT.DoF,NumQuadL)
  fTRefY  = zeros(hFeF.Comp,FeT.DoF,NumQuadL)
  @inbounds for iQ = 1 : NumQuadL
    @inbounds for iD = 1 : FeT.DoF
      fTRefX[1,iD,iQ] = FeT.phi[iD,1](-1.0,PointsL[iQ])
      fTRefY[1,iD,iQ] = FeT.phi[iD,2](PointsL[iQ],-1.0)
    end
    @inbounds for iComp = 1 : uFeF.Comp
      @inbounds for iD = 1 : uFeF.DoF
        uFFRefX[iComp,iD,iQ] = uFeF.phi[iD,iComp](-1.0,PointsL[iQ])
        uFFRefY[iComp,iD,iQ] = uFeF.phi[iD,iComp](PointsL[iQ],-1.0)
      end  
    end
    @inbounds for iD = 1 : hFeF.DoF
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
    @inbounds for iDoFuFeF = 1 : uFeF.DoF
      ind = uFeF.Glob[iDoFuFeF,iF]  
      uLoc[iDoFuFeF] = u[ind]
    end  
    @inbounds for iQ = 1 : NumQuad
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,Points[iQ,1],Points[iQ,2],Grid.Faces[iF], Grid)
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
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,-1.0,Points[iQ,2],Grid.Faces[iF], Grid)
      @. uFLoc = 0
      @inbounds for iDoFuFeF = 1 : uFeF.DoF
        uFLoc[1] += uFFRefX[1,iDoFuFeF,iQ] * uLoc[iDoFuFeF]
        uFLoc[2] += uFFRefX[2,iDoFuFeF,iQ] * uLoc[iDoFuFeF]
      end
      u1 = (DF[1,1] * uFLoc[1] + DF[1,2] * uFLoc[2]) / detDF[1]
      u2 = (DF[2,1] * uFLoc[1] + DF[2,2] * uFLoc[2]) / detDF[1]
      KLocX[iQ] = 0.5 * (u1 * u1 + u2 * u2)
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,Points[iQ,1],-1.0,Grid.Faces[iF], Grid)
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
    @inbounds for iQ = 1 : NumQuad
      @inbounds for iDoFFeT = 1 : FeT.DoF
        GradLoc[iDoFFeT] +=  Grid.Faces[iF].Orientation * Weights[iQ] * fTRef[1,iDoFFeT,iQ] * (hFLoc[iQ] + KLoc[iQ])
      end
    end
    @inbounds for iQ = 1 : NumQuadL
      @inbounds for iDoFFeT = 1 : FeT.DoF
        GradLoc[iDoFFeT] +=  WeightsL[iQ] * (fTRefX[1,iDoFFeT,iQ] * (hFLocX[iQ] + KLocX[iQ]) +
          fTRefY[1,iDoFFeT,iQ] * (hFLocY[iQ] + KLocY[iQ]))
      end
    end   
    @inbounds for iDoFFeT = 1 : FeT.DoF
      ind = FeT.Glob[iDoFFeT,iF]
      Grad[ind] += GradLoc[iDoFFeT]
    end
  end
end

"""
  CrossRhs!(backend, FTB, Cross, FeT::HDivElement, u, uFeF::HDivElement, Grid,
        ElemType::Grids.ElementType, QuadOrd, Jacobi)

Assembles the right-hand side vector contribution for the cross product (Coriolis) term in a finite element method (FEM) context, specifically for HDiv-conforming elements on a given grid.

# Arguments
- `backend`: Computational backend (e.g., for parallelization or hardware abstraction).
- `FTB`: (Unused in this function, possibly for future extension or interface compatibility).
- `Cross`: Vector to accumulate the assembled right-hand side contributions.
- `FeT::HDivElement`: Test finite element (HDiv-conforming) structure.
- `u`: Global solution vector (degrees of freedom).
- `uFeF::HDivElement`: Trial finite element (HDiv-conforming) structure.
- `Grid`: Grid structure containing mesh and face information.
- `ElemType::Grids.ElementType`: Element type descriptor for quadrature rule selection.
- `QuadOrd`: Quadrature order for numerical integration.
- `Jacobi`: Function to compute the Jacobian, its determinant, and pseudo-inverse for mapping reference to physical elements.

# Description
For each face in the grid, the function:
1. Evaluates basis functions at quadrature points.
2. Maps local degrees of freedom to global solution vector.
3. Computes the mapped solution at quadrature points.
4. Applies the Coriolis operator (cross product with a vector field `k`).
5. Integrates the result against test functions using quadrature.
6. Assembles the local contributions into the global right-hand side vector `Cross`.

# Notes
- Assumes 3D vector fields and HDiv-conforming elements.
- The Coriolis operator is computed via the `Coriolis` function, which returns a scalar field and a vector `k`.
- The function is performance-oriented, using `@inbounds` for loop optimization.

# See Also
- `Coriolis`
- `QuadRule`
- `HDivElement`
"""

function CrossRhs!(backend,FTB,Cross,FeT::HDivElement,u,uFeF::HDivElement,Grid,
  ElemType::Grids.ElementType,QuadOrd,Jacobi)

  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  fTRef  = zeros(FeT.Comp,FeT.DoF,NumQuad)

  @inbounds for iQ = 1 : NumQuad
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : FeT.DoF
        fTRef[iComp,iD,iQ] = FeT.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
  end
  uFFRef = fTRef

  CrossLoc = zeros(FeT.DoF)
  uLoc = zeros(uFeF.DoF)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)

  @inbounds for iF = 1 : Grid.NumFaces
    @. CrossLoc = 0
    @inbounds for iDoF = 1 : uFeF.DoF
      ind = uFeF.Glob[iDoF,iF]  
      uLoc[iDoF] = u[ind]
    end  
    @inbounds for iQ = 1 : NumQuad
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,Points[iQ,1],Points[iQ,2],Grid.Faces[iF], Grid)
      uFLoc1 = 0.0
      uFLoc2 = 0.0
      @inbounds for iDoF = 1 : uFeF.DoF
        uFLoc1 += uFFRef[1,iDoF,iQ] * uLoc[iDoF]
        uFLoc2 += uFFRef[2,iDoF,iQ] * uLoc[iDoF]
      end
      u1 = (DF[1,1] * uFLoc1 + DF[1,2] * uFLoc2) 
      u2 = (DF[2,1] * uFLoc1 + DF[2,2] * uFLoc2)
      u3 = (DF[3,1] * uFLoc1 + DF[3,2] * uFLoc2)

      qFLoc,k= Coriolis(X,Grid.Form)
      
      kcru1 = k[2] * u3 - k[3] * u2
      kcru2 = -(k[1] * u3 - k[3] * u1)
      kcru3 = k[1] * u2 - k[2] * u1
      cu1 = (DF[1,1] * kcru1 + DF[2,1] * kcru2 + DF[3,1] * kcru3)
      cu2 = (DF[1,2] * kcru1 + DF[2,2] * kcru2 + DF[3,2] * kcru3)
      @inbounds for iDoF = 1 : FeT.DoF
        CrossLoc[iDoF] -= Weights[iQ] * qFLoc * (fTRef[1,iDoF,iQ] * cu1 +
          fTRef[2,iDoF,iQ] * cu2) / detDF[1] 
      end    
    end  
    @inbounds for iDoF = 1 : FeT.DoF
      ind = FeT.Glob[iDoF,iF]
      Cross[ind] += CrossLoc[iDoF] 
    end
  end
end

"""
  Coriolis(X, Form)

Compute the local Coriolis parameter and unit vector for a given position.

# Arguments
- `X::AbstractVector{<:Real}`: A 3-element vector representing the position in Cartesian coordinates.
- `Form::String`: Specifies the geometry. Use `"Sphere"` for spherical geometry, otherwise a default planar value is used.

# Returns
- `qFLoc::Float64`: The local Coriolis parameter.
- `k::SVector{3, Float64}`: The local unit vector (vertical direction).

# Description
If `Form` is `"Sphere"`, computes the Coriolis parameter based on the latitude derived from the position vector `X` on a sphere (assuming Earth-like rotation). Otherwise, returns a constant Coriolis parameter and a vertical unit vector for planar geometry.
"""

function Coriolis(X,Form)
  if Form == "Sphere"
    Omega = 2 * pi / 24.0 / 3600.0
    Rad = sqrt(X[1]^2 + X[2]^2 + X[3]^2)
    k1 = X[1] / Rad
    k2 = X[2] / Rad
    k3 = X[3] / Rad
    sinlat = k3
    qFLoc = 2 * Omega * sinlat
    return qFLoc, SVector{3}(k1,k2,k3)
  else
    k1= 0.0
    k2= 0.0
    k3= 1.0
    return 1.e-4,SVector{3}(k1,k2,k3)
  end    
end

"""
  CrossRhs!(backend, FTB, Cross, q, qFeF::ScalarElement, u, uFeF::HDivElement, FeT::HDivElement, Grid,
        ElemType::Grids.ElementType, QuadOrd, Jacobi)

Assembles the right-hand side vector for the cross product term in a finite element method (FEM) context, specifically for mixed finite element spaces involving scalar and H(div) elements.

# Arguments
- `backend`: Computational backend to use (e.g., CPU, GPU).
- `FTB`: (Purpose not specified in code snippet; likely a workspace or buffer).
- `Cross`: Output vector to accumulate the cross product right-hand side contributions.
- `q`: Coefficient vector for the scalar field.
- `qFeF::ScalarElement`: Scalar finite element descriptor for `q`.
- `u`: Coefficient vector for the H(div) field.
- `uFeF::HDivElement`: H(div) finite element descriptor for `u`.
- `FeT::HDivElement`: H(div) finite element descriptor for the test space.
- `Grid`: Grid or mesh structure containing geometry and topology information.
- `ElemType::Grids.ElementType`: Element type descriptor for quadrature rule selection.
- `QuadOrd`: Quadrature order for numerical integration.
- `Jacobi`: Function to compute the Jacobian, its determinant, and pseudo-inverse for mapping reference to physical space.

# Description
For each face in the grid, this function:
- Evaluates basis functions and their values at quadrature points.
- Maps local degrees of freedom to global vectors.
- Computes the mapped vector field and its cross product with a Coriolis-like vector.
- Integrates the resulting expression using quadrature and accumulates the contributions to the global right-hand side vector.

# Notes
- The function is performance-critical and uses in-place operations and preallocated buffers.
- Assumes that the finite element descriptors (`qFeF`, `uFeF`, `FeT`) provide basis function evaluation and global-to-local mapping.
- The `Coriolis` function is expected to return a vector field and its associated vector at a given point.

# In-place
Modifies the `Cross` vector in-place.

# See also
- `QuadRule`
- `Coriolis`
- `Jacobi`
"""

function CrossRhs!(backend,FTB,Cross,q,qFeF::ScalarElement,u,uFeF::HDivElement,FeT::HDivElement,Grid,
  ElemType::Grids.ElementType,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  qFFRef  = zeros(qFeF.Comp,qFeF.DoF,NumQuad)
  fTRef  = zeros(FeT.Comp,FeT.DoF,NumQuad)

  @inbounds for iQ = 1 : NumQuad
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : FeT.DoF
        fTRef[iComp,iD,iQ] = FeT.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    @inbounds for iComp = 1 : qFeF.Comp
      @inbounds for iD = 1 : qFeF.DoF
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
  pinvDF = zeros(3,2)
  X = zeros(3)

  @inbounds for iF = 1 : Grid.NumFaces
    @. CrossLoc = 0
    @inbounds for iDoF = 1 : uFeF.DoF
      ind = uFeF.Glob[iDoF,iF]  
      uLoc[iDoF] = u[ind]
    end  
    @inbounds for iDoF = 1 : qFeF.DoF
      ind = qFeF.Glob[iDoF,iF]  
      qLoc[iDoF] = q[ind]
    end  
    @inbounds for iQ = 1 : NumQuad
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,Points[iQ,1],Points[iQ,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      uFLoc1 = 0.0
      uFLoc2 = 0.0
      @inbounds for iDoF = 1 : uFeF.DoF
        uFLoc1 += uFFRef[1,iDoF,iQ] * uLoc[iDoF]
        uFLoc2 += uFFRef[2,iDoF,iQ] * uLoc[iDoF]
      end
      u1 = (DF[1,1] * uFLoc1 + DF[1,2] * uFLoc2) 
      u2 = (DF[2,1] * uFLoc1 + DF[2,2] * uFLoc2)
      u3 = (DF[3,1] * uFLoc1 + DF[3,2] * uFLoc2)
      
      qFLoc,k= Coriolis(X,Grid.Form)

      kcru1 = k[2] * u3 - k[3] * u2
      kcru2 = -(k[1] * u3 - k[3] * u1)
      kcru3 = k[1] * u2 - k[2] * u1
      cu1 = (DF[1,1] * kcru1 + DF[2,1] * kcru2 + DF[3,1] * kcru3)
      cu2 = (DF[1,2] * kcru1 + DF[2,2] * kcru2 + DF[3,2] * kcru3)    
      @inbounds for iDoF = 1 : qFeF.DoF
        qFLoc += qFFRef[1,iDoF,iQ] * qLoc[iDoF]
      end
      @inbounds for iDoF = 1 : FeT.DoF
        CrossLoc[iDoF] -= Weights[iQ] * qFLoc * (fTRef[1,iDoF,iQ] * cu1 +
          fTRef[2,iDoF,iQ] * cu2) / detDFLoc
      end    
    end  
    @inbounds for iDoF = 1 : FeT.DoF
      ind = FeT.Glob[iDoF,iF]
      Cross[ind] += CrossLoc[iDoF]
    end
  end
end

"""
  DivMomentumVector!(backend, FTB, Rhs, FeTHDiv::HDivElement, uHDiv, FeHDiv::HDivElement,
            uVecDG, FeVecDG::VectorElement, Grid, ElemType::Grids.ElementType,
            QuadOrd, Jacobi)

Assembles the right-hand side (Rhs) vector for the divergence of the momentum equation in a finite element method (FEM) context, using mixed HDiv and DG elements.

# Arguments
- `backend`: Computational backend (e.g., CPU/GPU, not used directly in this function).
- `FTB`: (Unused in this function, possibly for future extension or interface compatibility).
- `Rhs`: Array to be updated in-place with the assembled right-hand side vector.
- `FeTHDiv::HDivElement`: Test function finite element space (H(div) element).
- `uHDiv`: Coefficient vector for the H(div) field.
- `FeHDiv::HDivElement`: Trial function finite element space (H(div) element).
- `uVecDG`: Coefficient vector for the DG vector field.
- `FeVecDG::VectorElement`: Finite element space for the DG vector field.
- `Grid`: Grid structure containing mesh and topology information.
- `ElemType::Grids.ElementType`: Element type (e.g., triangle or quadrilateral).
- `QuadOrd`: Quadrature order for numerical integration.
- `Jacobi`: Function to compute the Jacobian and related geometric quantities.

# Description
This function computes the local contributions to the right-hand side vector for the divergence of the momentum equation, including both element and edge (upwind) terms. It supports both triangular and quadrilateral elements and handles the assembly of contributions from volume integrals and edge integrals (for upwinding).

The function:
- Evaluates basis functions and their derivatives at quadrature points.
- Computes local variables and geometric transformations.
- Assembles local contributions to the global right-hand side vector, including upwind fluxes across edges.

# Notes
- The function assumes that the finite element spaces and grid structures are properly initialized and compatible.
- The function is performance-oriented, using in-place updates and preallocated arrays.
- The upwind flux computation distinguishes between left and right elements at each edge and applies the appropriate sign based on the upwind direction.

# See Also
- `QuadRule`
- `HDivElement`
- `VectorElement`
- `Grids.ElementType`
- `Jacobi`
"""

function DivMomentumVector!(backend,FTB,Rhs,FeTHDiv::HDivElement,uHDiv,FeHDiv::HDivElement,
  uVecDG,FeVecDG::VectorElement,Grid,ElemType::Grids.ElementType,QuadOrd,Jacobi)

  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  fuHDiv  = zeros(FeHDiv.DoF,FeHDiv.Comp,NumQuad)
  fuVecDG  = zeros(FeVecDG.DoF,FeVecDG.Comp,NumQuad)
  fuTHDiv  = zeros(FeTHDiv.DoF,FeTHDiv.Comp,NumQuad)
  fuDivHDiv = zeros(FeHDiv.DoF,NumQuad)
  fuGradVecDG = zeros(FeVecDG.DoF,FeVecDG.Comp,2,NumQuad)

  #computation of the ansatz functions in the quadrature points
  @inbounds for iQ = 1 : NumQuad
    @inbounds for iComp = 1 : FeTHDiv.Comp
      @inbounds for iD = 1 : FeTHDiv.DoF
        fuTHDiv[iD,iComp,iQ] = FeTHDiv.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    @inbounds for iComp = 1 : FeTHDiv.Comp
      @inbounds for iD = 1 : FeHDiv.DoF
        fuHDiv[iD,iComp,iQ] = FeHDiv.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    @inbounds for iD = 1 : FeVecDG.DoF
      fuVecDG[iD,1,iQ] = FeVecDG.phi[iD,1](Points[iQ,1],Points[iQ,2])
      fuVecDG[iD,2,iQ] = FeVecDG.phi[iD,2](Points[iQ,1],Points[iQ,2])
      fuVecDG[iD,3,iQ] = FeVecDG.phi[iD,3](Points[iQ,1],Points[iQ,2])
    end
    @inbounds for iD = 1 : FeHDiv.DoF
      fuDivHDiv[iD,iQ] = FeHDiv.Divphi[iD,1](Points[iQ,1],Points[iQ,2])
    end
    
    @inbounds for iD = 1 : FeVecDG.DoF
      fuGradVecDG[iD,1,1,iQ] = FeVecDG.Gradphi[iD,1,1](Points[iQ,1],Points[iQ,2])
      fuGradVecDG[iD,1,2,iQ] = FeVecDG.Gradphi[iD,1,2](Points[iQ,1],Points[iQ,2])
      fuGradVecDG[iD,2,1,iQ] = FeVecDG.Gradphi[iD,2,1](Points[iQ,1],Points[iQ,2])
      fuGradVecDG[iD,2,2,iQ] = FeVecDG.Gradphi[iD,2,2](Points[iQ,1],Points[iQ,2])
      fuGradVecDG[iD,3,1,iQ] = FeVecDG.Gradphi[iD,3,1](Points[iQ,1],Points[iQ,2])
      fuGradVecDG[iD,3,2,iQ] = FeVecDG.Gradphi[iD,3,2](Points[iQ,1],Points[iQ,2])
    end
  end

  #local variables
  uVecDGLoc = zeros(FeVecDG.DoF) 
  uuVecDGLoc = zeros(3) 
  uuGradVecDGLoc = zeros(3,2) 
  uHDivLoc = zeros(FeHDiv.DoF)
  uuHDivLoc = zeros(2)
  RhsLoc = zeros(FeTHDiv.DoF)

  GradVecDGHDiv = zeros(3)
  VecDGDivHDiv = zeros(3)
  
  #Jacobi-Allocation
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)

  DFL = zeros(3,2)
  detDFL = zeros(1)
  pinvDFL = zeros(3,2)
  XL = zeros(3)

  DFR = zeros(3,2)
  detDFR = zeros(1)
  pinvDFR = zeros(3,2)
  XR = zeros(3)
  innersum = zeros(3)

  @inbounds for iF = 1 : Grid.NumFaces
    @inbounds for iD = 1 : FeVecDG.DoF
      ind = FeVecDG.Glob[iD,iF]  
      uVecDGLoc[iD] = uVecDG[ind]
    end  
    @inbounds for iD = 1 : FeHDiv.DoF
      ind = FeHDiv.Glob[iD,iF]  
      uHDivLoc[iD] = uHDiv[ind]
    end  
    @. RhsLoc = 0.0

    @inbounds for iQ = 1 : NumQuad
      @. uuVecDGLoc = 0.0
      @. uuGradVecDGLoc = 0.0
      #computation of local variables in a quadrature point
      @inbounds for iD = 1 : FeVecDG.DoF
        @inbounds for i = 1 : 3
          uuVecDGLoc[i] += uVecDGLoc[iD] * fuVecDG[iD,i,iQ] 
          @inbounds for j = 1 : 2  
            uuGradVecDGLoc[i,j] += uVecDGLoc[iD] * fuGradVecDG[iD,i,j,iQ]  
          end
        end
      end
      uuDivHDivLoc = 0.0
      @. uuHDivLoc = 0.0
      @inbounds for iD = 1 : FeHDiv.DoF
        uuDivHDivLoc += uHDivLoc[iD] * fuDivHDiv[iD,iQ] 
        @inbounds for i = 1 : 2
          uuHDivLoc[i] += uHDivLoc[iD] * fuHDiv[iD,i,iQ]
        end 
      end

      @inbounds for i = 1 : 3
        GradVecDGHDiv[i] = 0.0  
        VecDGDivHDiv[i] = uuDivHDivLoc * uuVecDGLoc[i]
        @inbounds for j = 1 : 2  
          GradVecDGHDiv[i] += uuGradVecDGLoc[i,j] * uuHDivLoc[j]
        end
      end  
      @. innersum = (VecDGDivHDiv + GradVecDGHDiv)
      #computation of Jacobi
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,Points[iQ,1],Points[iQ,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      @inbounds for iD = 1 : FeTHDiv.DoF
        product = 0.0
        @inbounds for i = 1 : 3
          temp = 0.0  
          @inbounds for j = 1 : 2  
            temp +=  DF[i,j] * fuTHDiv[iD,j,iQ] / detDFLoc 
          end  
          product += innersum[i] * temp
        end  
        RhsLoc[iD] += product * Weights[iQ]
      end 
    end
    @inbounds for iD = 1 : FeTHDiv.DoF
      ind = FeTHDiv.Glob[iD,iF]
      Rhs[ind] -= RhsLoc[iD]
    end  
  end 

  #UPWIND on Edges
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  #uFFLocL = zeros(FeHDiv.DoF)
  #uFFLocR = zeros(FeHDiv.DoF)
  uHDivLocLeft = zeros(FeHDiv.DoF)
  uHDivLocRight = zeros(FeHDiv.DoF)
  uVecDGLocLeft = zeros(FeVecDG.DoF)
  uVecDGLocRight = zeros(FeVecDG.DoF)

  #distinction between Tri or Quad 
  if ElemType == Grids.Tri()
    nBar1 = [ 0 1 -1
             -1 1  0]
    PointsE = zeros(2,NumQuadL,3)
    @inbounds for iQ = 1 : NumQuadL
      PointsE[1,iQ,1] = PointsL[iQ]
      PointsE[2,iQ,1] = -1.0
      PointsE[1,iQ,2] = -PointsL[iQ]
      PointsE[2,iQ,2] = PointsL[iQ]
      PointsE[1,iQ,3] = -1.0
      PointsE[2,iQ,3] = PointsL[iQ]
    end
  elseif ElemType == Grids.Quad()
    nBar1 = [ 0 1 0 -1
             -1 0 1  0]
    PointsE = zeros(2,NumQuadL,4)
    @inbounds for iQ = 1 : NumQuadL
      PointsE[1,iQ,1] = PointsL[iQ]
      PointsE[2,iQ,1] = -1.0
      PointsE[1,iQ,2] = 1.0
      PointsE[2,iQ,2] = PointsL[iQ]
      PointsE[1,iQ,3] = PointsL[iQ]
      PointsE[2,iQ,3] = 1.0
      PointsE[1,iQ,4] = -1.0
      PointsE[2,iQ,4] = PointsL[iQ]
    end
  else 
    @show "wrong ElemType"
  end
  uFFRef  = zeros(FeHDiv.DoF, FeHDiv.Comp,NumQuadL,size(PointsE,3))
  uVecDGRef  = zeros(FeVecDG.DoF,FeVecDG.Comp,NumQuadL,size(PointsE,3))
  fTRef  = zeros(FeTHDiv.DoF,FeTHDiv.Comp,NumQuadL,size(PointsE,3))
  @inbounds for i = 1 : size(PointsE,3)
    @inbounds for iQ = 1 : NumQuadL 
      @inbounds for iComp = 1 : FeTHDiv.Comp
        @inbounds for iD = 1 : FeTHDiv.DoF
          fTRef[iD,iComp,iQ,i] = FeTHDiv.phi[iD,iComp](PointsE[1,iQ,i],PointsE[2,iQ,i])
        end
      end
      @inbounds for iComp = 1 : FeHDiv.Comp
        @inbounds for iD = 1 : FeHDiv.DoF
          uFFRef[iD,iComp,iQ,i] = FeHDiv.phi[iD,iComp](PointsE[1,iQ,i],PointsE[2,iQ,i])
        end
      end
      @inbounds for iComp = 1 : FeVecDG.Comp
        @inbounds for iD = 1 : FeVecDG.DoF
          uVecDGRef[iD,iComp,iQ,i] = FeVecDG.phi[iD,iComp](PointsE[1,iQ,i],PointsE[2,iQ,i])
        end
      end
    end
  end

  #allocation Mloc
  ILL = zeros(FeTHDiv.DoF)
  ILR = zeros(FeTHDiv.DoF)
  IRL = zeros(FeTHDiv.DoF)
  IRR = zeros(FeTHDiv.DoF)
  fTLocL = zeros(3,FeTHDiv.DoF)
  fTLocR = zeros(3,FeTHDiv.DoF)

  uL = zeros(3)
  uR = zeros(3)
  gammaU = 0.5
  @inbounds for iE in 1:Grid.NumEdges
    Edge = Grid.Edges[iE]
    if length(Edge.F) > 1
      iFL = Edge.F[1]
      EdgeTypeL = Edge.FE[1]
      iFR = Edge.F[2]
      EdgeTypeR = Edge.FE[2]
      #computation normales of edges
      @views nBarLocL = Grid.nBar[:, EdgeTypeL]
      @views nBarLocR = Grid.nBar[:, EdgeTypeR]
      @views nBarLocL = nBar1[:, EdgeTypeL] * Grid.Faces[iFL].Orientation
      @views nBarLocR = nBar1[:, EdgeTypeR] * Grid.Faces[iFR].Orientation
      #gamma upwind value
      uE = 0.0 
      uER = 0.0 
      @inbounds for iD = 1 : FeHDiv.DoF
        ind = FeHDiv.Glob[iD,iFL]  
        uHDivLocLeft[iD] = uHDiv[ind]
        ind = FeHDiv.Glob[iD,iFR]  
        uHDivLocRight[iD] = uHDiv[ind]
      end 
      @inbounds for iD = 1 : FeVecDG.DoF
        ind = FeVecDG.Glob[iD,iFL]  
        uVecDGLocLeft[iD] = uVecDG[ind]
        ind = FeVecDG.Glob[iD,iFR]  
        uVecDGLocRight[iD] = uVecDG[ind]
      end 
      @inbounds for iQ = 1:NumQuadL
        @inbounds for iD = 1 : FeHDiv.DoF
          @views uE += WeightsL[iQ]*(uHDivLocLeft[iD]*uFFRef[iD,:,iQ,EdgeTypeL]'*nBarLocL) 
        end
      end
      gammaLoc = uE > 0 ? gammaU : -gammaU

      @. ILL = 0
      @. ILR = 0
      @. IRL = 0
      @. IRR = 0

      @inbounds for iQ in 1: NumQuadL
        #computation of Jacobi EdgetypeL
        Jacobi(DFL,detDFL,pinvDFL,XL,Grid.Type,PointsE[1,iQ,EdgeTypeL],PointsE[2,iQ,EdgeTypeL],
          Grid.Faces[iFL], Grid)
        detDFLLoc = detDFL[1] * Grid.Faces[iFL].Orientation
        #EdgeTypeR
        Jacobi(DFR,detDFR,pinvDFR,XR,Grid.Type,PointsE[1,iQ,EdgeTypeR],PointsE[2,iQ,EdgeTypeR],
          Grid.Faces[iFR], Grid)
        detDFRLoc = detDFR[1] * Grid.Faces[iFR].Orientation

        @inbounds for iDoF  = 1 : FeTHDiv.DoF 
          @inbounds for j = 1 : 3
            fTLocL[j,iDoF] = (1/detDFLLoc) * (DFL[j,1] * fTRef[iDoF,1,iQ,EdgeTypeL] +
              DFL[j,2] * fTRef[iDoF,2,iQ,EdgeTypeL])
            fTLocR[j,iDoF] = (1/detDFRLoc) * (DFR[j,1] * fTRef[iDoF,1,iQ,EdgeTypeR] +
              DFR[j,2] * fTRef[iDoF,2,iQ,EdgeTypeR])
          end
        end
        uFFLocL = 0.0
        uFFLocR = 0.0
        @inbounds for iD in 1:FeHDiv.DoF
          @views uFFLocL += uHDivLocLeft[iD] * uFFRef[iD,:,iQ,EdgeTypeL]' * nBarLocL
          @views uFFLocR += uHDivLocRight[iD] * uFFRef[iD,:,iQ,EdgeTypeR]' * nBarLocR
        end

        @. uL = 0.0
        @. uR = 0.0
        @inbounds for iD in 1:FeVecDG.DoF 
          uL[1] +=  uVecDGRef[iD,1,iQ,EdgeTypeL] * uVecDGLocLeft[iD]          
          uL[2] +=  uVecDGRef[iD,2,iQ,EdgeTypeL] * uVecDGLocLeft[iD]          
          uL[3] +=  uVecDGRef[iD,3,iQ,EdgeTypeL] * uVecDGLocLeft[iD]          
          uR[1] +=  uVecDGRef[iD,1,iQ,EdgeTypeR] * uVecDGLocRight[iD]
          uR[2] +=  uVecDGRef[iD,2,iQ,EdgeTypeR] * uVecDGLocRight[iD]
          uR[3] +=  uVecDGRef[iD,3,iQ,EdgeTypeR] * uVecDGLocRight[iD]
        end
        @inbounds for iDoF  = 1 : FeTHDiv.DoF 
          @inbounds for k = 1 : 3  
            ILL[iDoF] += WeightsL[iQ] * (fTLocL[k,iDoF]' * uL[k]) * uFFLocL
            ILR[iDoF] += WeightsL[iQ] * (fTLocL[k,iDoF]' * uR[k]) * uFFLocL
            IRL[iDoF] += WeightsL[iQ] * (fTLocR[k,iDoF]' * uL[k]) * uFFLocR
            IRR[iDoF] += WeightsL[iQ] * (fTLocR[k,iDoF]' * uR[k]) * uFFLocR
          end
        end
      end
      #all together over edges
      @inbounds for iD = 1 : FeTHDiv.DoF
#       indR = FeTHDiv.Glob[iD,iFR]
#       Rhs[indR] = 0.5 * (IRR[iD] - IRL[iD])
#       indL = FeTHDiv.Glob[iD,iFL]
#       Rhs[indL] = 0.5 * (ILL[iD] - ILR[iD])
        
        if gammaLoc > 0.0
          indR = FeTHDiv.Glob[iD,iFR]
          Rhs[indR] += -1.0 * IRL[iD]
          Rhs[indR] += +1.0 * IRR[iD]
        else
          indL = FeTHDiv.Glob[iD,iFL]
          Rhs[indL] += +1.0  * ILL[iD]
          Rhs[indL] += -1.0  * ILR[iD]
        end  
      end
    end
  end
end

"""
  CurlVel!(q, FeT, u, uFe::HDivElement, QuadOrd, ElemType, Grid, Jacobi)

Assembles the right-hand side vector `q` for the weak form of the curl operator applied to a velocity field `u` in a finite element context. Specifically, computes the integral

  ∫ (curl u) ⋅ v dx = -∫ u ⋅ (rot v) dx

using numerical quadrature over all faces of the grid.

# Arguments
- `q`: Output vector to be assembled (modified in-place).
- `FeT`: Test function finite element type, containing basis and gradient information.
- `u`: Coefficient vector representing the velocity field in the HDiv space.
- `uFe::HDivElement`: Finite element type for the HDiv-conforming space, providing basis functions for `u`.
- `QuadOrd`: Quadrature order for numerical integration.
- `ElemType`: Element type descriptor (e.g., triangle, quadrilateral).
- `Grid`: Grid structure containing mesh and face connectivity information.
- `Jacobi`: Function to compute the Jacobian, its determinant, and pseudo-inverse for mapping reference to physical elements.

# Details
- Initializes quadrature rules and evaluates basis functions and their gradients at quadrature points.
- Loops over all faces in the grid, assembling local contributions to `q` using the provided velocity field `u`.
- Applies the Jacobian transformation to map reference element quantities to the physical element.
- Accumulates contributions using the quadrature weights and basis function gradients.
- Solves the local system using the LU factorization stored in `FeT.LUM`.

# Notes
- Assumes 2D or 3D vector fields, with appropriate mapping and transformation.
- The function modifies `q` in-place.
- The function is performance-optimized using `@inbounds` and preallocated arrays.

# See also
- [`QuadRule`](@ref)
- [`HDivElement`](@ref)
- [`ldiv!`](@ref)
"""
function CurlVel!(q,FeT,u,uFe::HDivElement,QuadOrd,ElemType,Grid,Jacobi)
#
#
# int q*v dx = int Curl u * v dx = - int u * rot v dx 
#
#
  @. q = 0

  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)  
  RotqRef  = zeros(2,FeT.DoF,NumQuad)
  uFRef  = zeros(uFe.Comp,uFe.DoF,NumQuad)
  @inbounds for iQ = 1 : NumQuad
    @inbounds for iDoF = 1 : uFe.DoF
      uFRef[1,iDoF,iQ] = uFe.phi[iDoF,1](Points[iQ,1],Points[iQ,2])
      uFRef[2,iDoF,iQ] = uFe.phi[iDoF,2](Points[iQ,1],Points[iQ,2])
    end  
    @inbounds for iDoF = 1 : FeT.DoF
      RotqRef[1,iDoF,iQ] = -FeT.Gradphi[iDoF,1,2](Points[iQ,1],Points[iQ,2])
      RotqRef[2,iDoF,iQ] = FeT.Gradphi[iDoF,1,1](Points[iQ,1],Points[iQ,2])
    end  
  end  
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  uLoc = zeros(uFe.DoF)
  uFLoc = zeros(2)
  uuFLoc = zeros(3)
  qLoc = zeros(FeT.DoF)

  @inbounds for iF = 1 : Grid.NumFaces
    @inbounds for iDoF = 1 : uFe.DoF
      ind = uFe.Glob[iDoF,iF]  
      uLoc[iDoF] = u[ind] 
    end  
    @. qLoc = 0
    @inbounds for iQ = 1 : NumQuad
      @. uFLoc = 0  
      @inbounds for iDoF = 1 : uFe.DoF  
        uFLoc[1] += uFRef[1,iDoF,iQ] * uLoc[iDoF]  
        uFLoc[2] += uFRef[2,iDoF,iQ] * uLoc[iDoF]  
      end   
      Jacobi(DF,detDF,pinvDF,X,ElemType,Points[iQ,1],Points[iQ,2],Grid.Faces[iF], Grid)
      uuFLoc[1] = (DF[1,1] * uFLoc[1] + DF[1,2] * uFLoc[2]) / detDF[1]
      uuFLoc[2] = (DF[2,1] * uFLoc[1] + DF[2,2] * uFLoc[2]) / detDF[1]
      uuFLoc[3] = (DF[3,1] * uFLoc[1] + DF[3,2] * uFLoc[2]) / detDF[1]
      uFLoc[1] = DF[1,1] * uuFLoc[1] + DF[2,1] * uuFLoc[2] + DF[3,1] * uuFLoc[3]
      uFLoc[2] = DF[1,2] * uuFLoc[1] + DF[2,2] * uuFLoc[2] + DF[3,2] * uuFLoc[3]
      @inbounds for iDoF = 1 : FeT.DoF  
        qLoc[iDoF] +=  -Weights[iQ] * (uFLoc[1] * RotqRef[1,iDoF,iQ] +
          uFLoc[2] * RotqRef[2,iDoF,iQ])  
      end  
    end  
    @inbounds for iDoF = 1 : FeT.DoF
      ind = FeT.Glob[iDoF,iF]
      q[ind] += qLoc[iDoF]
    end
  end    
  ldiv!(FeT.LUM,q)
end

"""
  CurlVel!(q, FeT, u, uFe::HDivElement, h, hFe::ScalarElement, QuadOrd, ElemType, Grid, Jacobi)

Assembles the right-hand side vector `q` for the weak form of the curl-velocity operator in a finite element method (FEM) context. Specifically, it computes the integral

  ∫ q * v dx = ∫ (curl u) * v dx = -∫ u * rot(v) dx

over all faces of the computational grid, using the provided finite element spaces and quadrature rules.

# Arguments
- `q`: Output vector to be assembled (modified in-place).
- `FeT`: Test function finite element type (provides basis functions and global indexing).
- `u`: Input vector field (typically velocity).
- `uFe::HDivElement`: Finite element description for the velocity field (H(div) space).
- `h`: Scalar field (e.g., height or pressure).
- `hFe::ScalarElement`: Finite element description for the scalar field.
- `QuadOrd`: Quadrature order for numerical integration.
- `ElemType`: Element type (e.g., triangle, quadrilateral).
- `Grid`: Grid structure containing mesh and face information.
- `Jacobi`: Function to compute the Jacobian, its determinant, and pseudo-inverse for mapping reference to physical elements.

# Details
- Loops over all faces in the grid, assembling local contributions to the global vector `q`.
- Uses quadrature to evaluate integrals on each face.
- Applies the Piola transformation to map basis functions from reference to physical elements.
- Handles both vector and scalar finite element spaces.
- The assembled vector `q` is solved in-place using the LU decomposition stored in `FeT.LUM`.

# Notes
- Assumes that all input structures (`FeT`, `uFe`, `hFe`, `Grid`) are properly initialized and compatible.
- The function is performance-critical and uses in-place operations and `@inbounds` for efficiency.
"""

function CurlVel!(q,FeT,u,uFe::HDivElement,h,hFe::ScalarElement,QuadOrd,ElemType,Grid,Jacobi)
#
#
# int q*v dx = int Curl u * v dx = - int u * rot v dx 
#
#
  @. q = 0

  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)  
  RotqRef  = zeros(2,FeT.DoF,NumQuad)
  uFRef  = zeros(uFe.Comp,uFe.DoF,NumQuad)
  hFRef  = zeros(hFe.DoF,NumQuad)
  @inbounds for iQ = 1 : NumQuad
    @inbounds for iDoF = 1 : uFe.DoF
      uFRef[1,iDoF,iQ] = uFe.phi[iDoF,1](Points[iQ,1],Points[iQ,2])
      uFRef[2,iDoF,iQ] = uFe.phi[iDoF,2](Points[iQ,1],Points[iQ,2])
    end  
    @inbounds for iDoF = 1 : FeT.DoF
      RotqRef[1,iDoF,iQ] = -FeT.Gradphi[iDoF,1,2](Points[iQ,1],Points[iQ,2])
      RotqRef[2,iDoF,iQ] = FeT.Gradphi[iDoF,1,1](Points[iQ,1],Points[iQ,2])
    end  
    @inbounds for iDoF = 1 : hFe.DoF
      hFRef[iDoF,iQ] = hFe.phi[iDoF](Points[iQ,1],Points[iQ,2])
    end  
  end  
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  uLoc = zeros(uFe.DoF)
  hLoc = zeros(hFe.DoF)
  uFLoc = zeros(2)
  uuFLoc = zeros(3)
  qLoc = zeros(FeT.DoF)

  @inbounds for iF = 1 : Grid.NumFaces
    @inbounds for iDoF = 1 : uFe.DoF
      ind = uFe.Glob[iDoF,iF]  
      uLoc[iDoF] = u[ind] 
    end  
    @inbounds for iDoF = 1 : hFe.DoF
      ind = hFe.Glob[iDoF,iF]  
      hLoc[iDoF] = h[ind] 
    end  
    @. qLoc = 0
    @inbounds for iQ = 1 : NumQuad
      @. uFLoc = 0  
      @inbounds for iDoF = 1 : uFe.DoF  
        uFLoc[1] += uFRef[1,iDoF,iQ] * uLoc[iDoF]  
        uFLoc[2] += uFRef[2,iDoF,iQ] * uLoc[iDoF]  
      end   
      hFLoc = 0.0 
      @inbounds for iDoF = 1 : hFe.DoF  
        hFLoc += hFRef[iDoF,iQ] * hLoc[iDoF]  
      end   
      Jacobi(DF,detDF,pinvDF,X,ElemType,Points[iQ,1],Points[iQ,2],Grid.Faces[iF], Grid)
      uuFLoc[1] = (DF[1,1] * uFLoc[1] + DF[1,2] * uFLoc[2]) / detDF[1]
      uuFLoc[2] = (DF[2,1] * uFLoc[1] + DF[2,2] * uFLoc[2]) / detDF[1]
      uuFLoc[3] = (DF[3,1] * uFLoc[1] + DF[3,2] * uFLoc[2]) / detDF[1]
      uFLoc[1] = DF[1,1] * uuFLoc[1] + DF[2,1] * uuFLoc[2] + DF[3,1] * uuFLoc[3]
      uFLoc[2] = DF[1,2] * uuFLoc[1] + DF[2,2] * uuFLoc[2] + DF[3,2] * uuFLoc[3]
      @inbounds for iDoF = 1 : FeT.DoF  
        qLoc[iDoF] +=  -Weights[iQ] * (uFLoc[1] * RotqRef[1,iDoF,iQ] +
          uFLoc[2] * RotqRef[2,iDoF,iQ]) / hFLoc 
      end  
    end  
    @inbounds for iDoF = 1 : FeT.DoF
      ind = FeT.Glob[iDoF,iF]
      q[ind] += qLoc[iDoF]
    end
  end    
  ldiv!(FeT.LUM,q)
end

"""
  DivMomentumVectorOld!(backend, FTB, Rhs, FeTHDiv::HDivElement, uHDiv, FeHDiv::HDivElement,
             uVecDG, FeVecDG::VectorElement, Grid, ElemType::Grids.ElementType,
             QuadOrd, Jacobi)

Assembles the right-hand side (Rhs) vector for the divergence of the momentum equation in a finite element method (FEM) context, using mixed HDiv and DG vector elements. This function computes local contributions on faces and edges of the mesh and accumulates them into the global Rhs vector.

# Arguments
- `backend`: Computational backend (e.g., CPU/GPU, not used directly in this function).
- `FTB`: (Purpose not specified in the code snippet; likely a data structure for finite element assembly).
- `Rhs`: The global right-hand side vector to be assembled (modified in-place).
- `FeTHDiv::HDivElement`: Test function finite element (H(div) space).
- `uHDiv`: Coefficient vector for H(div) solution.
- `FeHDiv::HDivElement`: Trial function finite element (H(div) space).
- `uVecDG`: Coefficient vector for DG vector solution.
- `FeVecDG::VectorElement`: DG vector finite element.
- `Grid`: Mesh/grid data structure containing faces, edges, and geometric information.
- `ElemType::Grids.ElementType`: Element type (e.g., triangle or quadrilateral).
- `QuadOrd`: Quadrature order for numerical integration.
- `Jacobi`: Function to compute the Jacobian, determinant, and related geometric quantities for mapping reference to physical elements.

# Description
- Evaluates basis functions and their derivatives at quadrature points for both H(div) and DG vector elements.
- Loops over all faces to compute local contributions to the Rhs vector using quadrature.
- Computes upwind fluxes and edge contributions for internal edges, handling both triangle and quadrilateral elements.
- Assembles all local contributions into the global Rhs vector, applying upwind stabilization on edges.

# Notes
- The function is performance-oriented, using `@inbounds` and preallocated arrays to minimize allocations.
- Assumes specific data structures for finite elements and grid, including global-to-local DoF mappings and geometric information.
- The function is tailored for 2D problems (triangles and quadrilaterals) and vector-valued fields.

# Returns
- Modifies `Rhs` in-place; does not return a value.

# See Also
- `QuadRule`
- `HDivElement`
- `VectorElement`
- `Grids.ElementType`
- `Jacobi`
"""

function DivMomentumVectorOld!(backend,FTB,Rhs,FeTHDiv::HDivElement,uHDiv,FeHDiv::HDivElement,
  uVecDG,FeVecDG::VectorElement,Grid,ElemType::Grids.ElementType,QuadOrd,Jacobi)

  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  fuHDiv  = zeros(FeHDiv.DoF,FeHDiv.Comp,NumQuad)
  fuVecDG  = zeros(FeVecDG.DoF,FeVecDG.DoF,NumQuad)
  fuTHDiv  = zeros(FeTHDiv.DoF,FeTHDiv.Comp,NumQuad)
  fuDivHDiv = zeros(FeHDiv.DoF,NumQuad)
  fuGradVecDG = zeros(FeVecDG.DoF,FeVecDG.DoF,2,NumQuad)

  #computation of the ansatz functions in the quadrature points
  @inbounds for iQ = 1 : NumQuad
    @inbounds for iComp = 1 : FeTHDiv.Comp
      @inbounds for iD = 1 : FeTHDiv.DoF
        fuTHDiv[iD,iComp,iQ] = FeTHDiv.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    @inbounds for iComp = 1 : FeTHDiv.Comp
      @inbounds for iD = 1 : FeHDiv.DoF
        fuHDiv[iD,iComp,iQ] = FeHDiv.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    @inbounds for iD = 1 : FeVecDG.DoF
      fuVecDG[iD,1,iQ] = FeVecDG.phi[iD,1](Points[iQ,1],Points[iQ,2])
      fuVecDG[iD,2,iQ] = FeVecDG.phi[iD,2](Points[iQ,1],Points[iQ,2])
      fuVecDG[iD,3,iQ] = FeVecDG.phi[iD,3](Points[iQ,1],Points[iQ,2])
    end
    @inbounds for iD = 1 : FeHDiv.DoF
      fuDivHDiv[iD,iQ] = FeHDiv.Divphi[iD,1](Points[iQ,1],Points[iQ,2])
    end
    
    @inbounds for iD = 1 : FeVecDG.DoF
      fuGradVecDG[iD,1,1,iQ] = FeVecDG.Gradphi[iD,1,1](Points[iQ,1],Points[iQ,2])
      fuGradVecDG[iD,1,2,iQ] = FeVecDG.Gradphi[iD,1,2](Points[iQ,1],Points[iQ,2])
      fuGradVecDG[iD,2,1,iQ] = FeVecDG.Gradphi[iD,2,1](Points[iQ,1],Points[iQ,2])
      fuGradVecDG[iD,2,2,iQ] = FeVecDG.Gradphi[iD,2,2](Points[iQ,1],Points[iQ,2])
      fuGradVecDG[iD,3,1,iQ] = FeVecDG.Gradphi[iD,3,1](Points[iQ,1],Points[iQ,2])
      fuGradVecDG[iD,3,2,iQ] = FeVecDG.Gradphi[iD,3,2](Points[iQ,1],Points[iQ,2])
    end
  end


  #local variables
  uVecDGLoc = zeros(FeVecDG.DoF) 
  uuVecDGLoc = zeros(3) 
  uuGradVecDGLoc = zeros(3,2) 
  uHDivLoc = zeros(FeHDiv.DoF)
  uuHDivLoc = zeros(2)
  RhsLoc = zeros(FeTHDiv.DoF)

  GradVecDGHDiv = zeros(3)
  VecDGDivHDiv = zeros(3)
  
  #Jacobi-Allocation
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)

  DFL = zeros(3,2)
  detDFL = zeros(1)
  pinvDFL = zeros(3,2)
  XL = zeros(3)

  DFR = zeros(3,2)
  detDFR = zeros(1)
  pinvDFR = zeros(3,2)
  XR = zeros(3)
  innersum = zeros(3)

  @inbounds for iF = 1 : Grid.NumFaces
    @inbounds for iD = 1 : FeVecDG.DoF
      ind = FeVecDG.Glob[iD,iF]  
      uVecDGLoc[iD] = uVecDG[ind]
    end  
    @inbounds for iD = 1 : FeHDiv.DoF
      ind = FeHDiv.Glob[iD,iF]  
      uHDivLoc[iD] = uHDiv[ind]
    end  
    @. RhsLoc = 0.0

    @inbounds for iQ = 1 : NumQuad
      @. uuVecDGLoc = 0.0
      @. uuGradVecDGLoc = 0.0
      #computation of local variables in a quadrature point
      @inbounds for iD = 1 : FeVecDG.DoF
        @inbounds for i = 1 : 3
          uuVecDGLoc[i] += uVecDGLoc[iD] * fuVecDG[iD,i,iQ] 
          @inbounds for j = 1 : 2  
            uuGradVecDGLoc[i,j] += uVecDGLoc[iD] * fuGradVecDG[iD,i,j,iQ]  
          end
        end
      end
      uuDivHDivLoc = 0.0
      @. uuHDivLoc = 0.0
      @inbounds for iD = 1 : FeHDiv.DoF
        uuDivHDivLoc += uHDivLoc[iD] * fuDivHDiv[iD,iQ] 
        @inbounds for i = 1 : 2
          uuHDivLoc[i] += uHDivLoc[iD] * fuHDiv[iD,i,iQ]
        end 
      end

      @inbounds for i = 1 : 3
        GradVecDGHDiv[i] = 0.0  
        VecDGDivHDiv[i] = uuDivHDivLoc * uuVecDGLoc[i]
        @inbounds for j = 1 : 2  
          GradVecDGHDiv[i] += uuGradVecDGLoc[i,j] * uuHDivLoc[j]
        end
      end  
#     @. innersum = Grid.Faces[iF].Orientation * (VecDGDivHDiv + GradVecDGHDiv)
      @. innersum = (VecDGDivHDiv + GradVecDGHDiv)
      #computation of Jacobi
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,Points[iQ,1],Points[iQ,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      @inbounds for iD = 1 : FeTHDiv.DoF
        product = 0.0
        @inbounds for i = 1 : 3
          temp = 0.0  
          @inbounds for j = 1 : 2  
            temp +=  DF[i,j] * fuTHDiv[iD,j,iQ] / detDFLoc 
          end  
          product += innersum[i] * temp
        end  
        RhsLoc[iD] += product * Weights[iQ]
      end 
    end
    @inbounds for iD = 1 : FeTHDiv.DoF
      ind = FeTHDiv.Glob[iD,iF]
      Rhs[ind] -= RhsLoc[iD]
    end  
  end 

  #UPWIND on Edges
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  uFFLocL = zeros(FeHDiv.DoF)
  uFFLocR = zeros(FeHDiv.DoF)
  uHDivLocLeft = zeros(FeHDiv.DoF)
  uHDivLocRight = zeros(FeHDiv.DoF)
  uVecDGLocLeft = zeros(FeVecDG.DoF)
  uVecDGLocRight = zeros(FeVecDG.DoF)

  #distinction between Tri or Quad 
  if ElemType == Grids.Tri()
    PointsE = zeros(2,NumQuadL,3)
    @inbounds for iQ = 1 : NumQuadL
      PointsE[1,iQ,1] = PointsL[iQ]
      PointsE[2,iQ,1] = -1.0
      PointsE[1,iQ,2] = -PointsL[iQ]
      PointsE[2,iQ,2] = PointsL[iQ]
      PointsE[1,iQ,3] = -1.0
      PointsE[2,iQ,3] = PointsL[iQ]
    end
  elseif ElemType == Grids.Quad()
    PointsE = zeros(2,NumQuadL,4)
    @inbounds for iQ = 1 : NumQuadL
      PointsE[1,iQ,1] = PointsL[iQ]
      PointsE[2,iQ,1] = -1.0
      PointsE[1,iQ,2] = 1.0
      PointsE[2,iQ,2] = PointsL[iQ]
      PointsE[1,iQ,3] = PointsL[iQ]
      PointsE[2,iQ,3] = 1.0
      PointsE[1,iQ,4] = -1.0
      PointsE[2,iQ,4] = PointsL[iQ]
    end
  else 
    @show "wrong ElemType"
  end
  uFFRef  = zeros(FeHDiv.DoF, FeHDiv.Comp,NumQuadL,size(PointsE,3))
  uVecDGRef  = zeros(FeVecDG.DoF,FeVecDG.Comp,NumQuadL,size(PointsE,3))
  fTRef  = zeros(FeTHDiv.DoF,FeTHDiv.Comp,NumQuadL,size(PointsE,3))
  @inbounds for i = 1 : size(PointsE,3)
    @inbounds for iQ = 1 : NumQuadL 
      @inbounds for iComp = 1 : FeTHDiv.Comp
        @inbounds for iD = 1 : FeTHDiv.DoF
          fTRef[iD,iComp,iQ,i] = FeTHDiv.phi[iD,iComp](PointsE[1,iQ,i],PointsE[2,iQ,i])
        end
      end
      @inbounds for iComp = 1 : FeHDiv.Comp
        @inbounds for iD = 1 : FeHDiv.DoF
          uFFRef[iD,iComp,iQ,i] = FeHDiv.phi[iD,iComp](PointsE[1,iQ,i],PointsE[2,iQ,i])
        end
      end
      @inbounds for iComp = 1 : FeVecDG.Comp
        @inbounds for iD = 1 : FeVecDG.DoF
          uVecDGRef[iD,iComp,iQ,i] = FeVecDG.phi[iD,iComp](PointsE[1,iQ,i],PointsE[2,iQ,i])
        end
      end
    end
  end

  #allocation Mloc
  MLoc11 = zeros(FeTHDiv.DoF, FeHDiv.DoF)
  MLoc12 = zeros(FeTHDiv.DoF, FeHDiv.DoF)
  MLoc21 = zeros(FeTHDiv.DoF, FeHDiv.DoF)
  MLoc22 = zeros(FeTHDiv.DoF, FeHDiv.DoF)
  fTLocL = zeros(3,FeTHDiv.DoF)
  fTLocR = zeros(3,FeTHDiv.DoF)

  uL = zeros(3)
  uR = zeros(3)
  @inbounds for iE in 1:Grid.NumEdges
    Edge = Grid.Edges[iE]
    if length(Edge.F) > 1
      iFL = Edge.F[1]
      EdgeTypeL = Edge.FE[1]
      iFR = Edge.F[2]
      EdgeTypeR = Edge.FE[2]
      #computation normales of edges
      @views nBarLocL = Grid.nBar[:, EdgeTypeL]
      @views nBarLocR = Grid.nBar[:, EdgeTypeR]
      #gamma upwind value
      uE = 0.0 
      uER = 0.0 
      @inbounds for iD = 1 : FeHDiv.DoF
        ind = FeHDiv.Glob[iD,iFL]  
        uHDivLocLeft[iD] = uHDiv[ind]
        ind = FeHDiv.Glob[iD,iFR]  
        uHDivLocRight[iD] = uHDiv[ind]
      end 
      @inbounds for iD = 1 : FeVecDG.DoF
        ind = FeVecDG.Glob[iD,iFL]  
        uVecDGLocLeft[iD] = uVecDG[ind]
        ind = FeVecDG.Glob[iD,iFR]  
        uVecDGLocRight[iD] = uVecDG[ind]
      end 
      @inbounds for iQ = 1:NumQuadL
        @inbounds for iD = 1 : FeHDiv.DoF
#         @views uE += WeightsL[iQ]*(uHDivLocLeft[iD]*(uFFRef[iD,:,iQ,EdgeTypeL]'*nBarLocL) +
#           uHDivLocRight[iD]*(uFFRef[iD,:,iQ,EdgeTypeR]'*nBarLocR))
          @views uE += WeightsL[iQ]*(uHDivLocLeft[iD]*uFFRef[iD,:,iQ,EdgeTypeL]'*nBarLocL) 
          @views uER += WeightsL[iQ]*uHDivLocRight[iD]*(uFFRef[iD,:,iQ,EdgeTypeR]'*nBarLocR)
        end
      end
      gammaU = 0.5
      gammaLoc = uE > 0 ? gammaU : -gammaU

      @. MLoc11 = 0
      @. MLoc12 = 0
      @. MLoc21 = 0
      @. MLoc22 = 0

      @inbounds for iQ in 1: NumQuadL
        #computation of Jacobi EdgetypeL
        Jacobi(DFL,detDFL,pinvDFL,XL,Grid.Type,PointsE[iQ,1,EdgeTypeL],PointsE[iQ,2,EdgeTypeL],Grid.Faces[iFL], Grid)
        detDFLLoc = detDFL[1] * Grid.Faces[iFL].Orientation
        #EdgeTypeR
        Jacobi(DFR,detDFR,pinvDFR,XR,Grid.Type,PointsE[iQ,1,EdgeTypeR],PointsE[iQ,2,EdgeTypeR],Grid.Faces[iFR], Grid)
        detDFRLoc = detDFR[1] * Grid.Faces[iFR].Orientation

        @inbounds for iDoF  = 1 : FeTHDiv.DoF 
          @inbounds for j = 1 : 3
            fTLocL[j,iDoF] = (1/detDFLLoc) * (DFL[j,1] * fTRef[iDoF,1,iQ,EdgeTypeL] +
              DFL[j,2] * fTRef[iDoF,2,iQ,EdgeTypeL])
            fTLocR[j,iDoF] = (1/detDFRLoc) * (DFR[j,1] * fTRef[iDoF,1,iQ,EdgeTypeR] +
              DFR[j,2] * fTRef[iDoF,2,iQ,EdgeTypeR])
          end
        end
        @inbounds for iD in 1:FeHDiv.DoF
          @views uFFLocL[iD] = uFFRef[iD,:,iQ,EdgeTypeL]' * nBarLocL
          @views uFFLocR[iD] = uFFRef[iD,:,iQ,EdgeTypeR]' * nBarLocR
        end
        @. uL = 0.0
        @. uR = 0.0
      
        @inbounds for iD in 1:FeVecDG.DoF 
          uL[1] +=  uVecDGRef[iD,1,iQ,EdgeTypeL] * uVecDGLocLeft[iD]          
          uL[2] +=  uVecDGRef[iD,2,iQ,EdgeTypeL] * uVecDGLocLeft[iD]          
          uL[3] +=  uVecDGRef[iD,3,iQ,EdgeTypeL] * uVecDGLocLeft[iD]          
          uR[1] +=  uVecDGRef[iD,1,iQ,EdgeTypeR] * uVecDGLocRight[iD]
          uR[2] +=  uVecDGRef[iD,2,iQ,EdgeTypeR] * uVecDGLocRight[iD]
          uR[3] +=  uVecDGRef[iD,3,iQ,EdgeTypeR] * uVecDGLocRight[iD]
        end
        @inbounds for iDoF  = 1 : FeTHDiv.DoF 
          @inbounds for jDoF  = 1 : FeHDiv.DoF 
            @inbounds for k = 1 : 3  
              MLoc11[iDoF,jDoF] += WeightsL[iQ] * (fTLocL[k,iDoF]' * uL[k]) * uFFLocL[jDoF]
              MLoc12[iDoF,jDoF] += WeightsL[iQ] * (fTLocL[k,iDoF]' * uR[k]) * uFFLocR[jDoF]
              MLoc21[iDoF,jDoF] += WeightsL[iQ] * (fTLocR[k,iDoF]' * uL[k]) * uFFLocL[jDoF]
              MLoc22[iDoF,jDoF] += WeightsL[iQ] * (fTLocR[k,iDoF]' * uR[k]) * uFFLocR[jDoF]
            end
          end
        end
      end
      #all together over edges
      @inbounds for iD = 1 : FeTHDiv.DoF
        indL = FeTHDiv.Glob[iD,iFL] 
        @views Rhs[indL] += (+0.5 - gammaLoc) * MLoc11[iD,:]' * uHDivLocLeft
        @views Rhs[indL] += (-0.5 + gammaLoc) * MLoc12[iD,:]' * uHDivLocRight
        indR = FeTHDiv.Glob[iD,iFR] 
        @views Rhs[indR] += (+0.5 + gammaLoc) * MLoc21[iD,:]' * uHDivLocLeft
        @views Rhs[indR] += (-0.5 - gammaLoc) * MLoc22[iD,:]' * uHDivLocRight
      end
    end
  end
end
