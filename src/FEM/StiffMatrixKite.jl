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
        fTRefX[iComp,iD,i] = FeT.phi[iD,iComp](1.0,PointsL[i])
        fTRefY[iComp,iD,i] = FeT.phi[iD,iComp](PointsL[i],1.0)
      end
    end
    @inbounds for iD = 1 : FeF.DoF
      fFRefX[1,iD,i] = FeF.phi[iD,1](1.0,PointsL[i])
      fFRefY[1,iD,i] = FeF.phi[iD,2](PointsL[i],1.0)
    end
  end

  RowInd = Int64[]
  ColInd = Int64[]
  Val = Float64[]
  DivLoc = zeros(FeT.DoF,FeF.DoF)

  @inbounds for iF = 1 : Grid.NumFaces
    DivLoc .= 0
    @inbounds for i = 1 : NumQuad
      DivLoc += Grid.Faces[iF].Orientation * Weights[i] * (fTRef[:,:,i]' * fFRef[:,:,i])
    end
    @inbounds for i = 1 : NumQuadL
       DivLoc -= WeightsL[i] * (fTRefX[:,:,i]' * fFRefX[:,:,i] +  
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

function GradRhs1!(backend,FTB,Grad,FeT::HDivKiteDElement,hF,hFeF::ScalarElement,Grid,ElemType::Grids.ElementType,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  hFFRef  = zeros(FeT.Comp,hFeF.DoF,NumQuad)
  fTRef  = zeros(FeT.Comp,FeT.DoF,NumQuad)

  @inbounds for i = 1 : NumQuad
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : FeT.DoF
        fTRef[iComp,iD,i] = FeT.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : hFeF.DoF
        hFFRef[iComp,iD,i] = hFeF.Gradphi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  fTRefX  = zeros(FeT.Comp,FeT.DoF,NumQuadL)
  fTRefY  = zeros(FeT.Comp,FeT.DoF,NumQuadL)
  hFFRefX  = zeros(hFeF.Comp,hFeF.DoF,NumQuadL)
  hFFRefY  = zeros(hFeF.Comp,hFeF.DoF,NumQuadL)
  @inbounds for i = 1 : NumQuadL
    @inbounds for iComp = 1 : hFeF.Comp
      @inbounds for iD = 1 : hFeF.DoF
        hFFRefX[iComp,iD,i] = hFeF.phi[iD,iComp](1.0,PointsL[i])
        hFFRefY[iComp,iD,i] = hFeF.phi[iD,iComp](PointsL[i],1.0)
      end
    end
    @inbounds for iD = 1 : FeT.DoF
      fTRefX[1,iD,i] = FeT.phi[iD,1](1.0,PointsL[i])
      fTRefY[1,iD,i] = FeT.phi[iD,2](PointsL[i],1.0)
    end
  end
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  GradLoc = zeros(FeT.DoF)
  hFLoc = zeros(2,NumQuad)
  hFLocX = zeros(NumQuadL)
  hFLocY = zeros(NumQuadL)

  @inbounds for iF = 1 : Grid.NumFaces
    @. GradLoc = 0
    @. hFLoc = 0
    @. hFLocX = 0
    @. hFLocY = 0
    @inbounds for iDoFhFeF = 1 : hFeF.DoF
      ind = hFeF.Glob[iDoFhFeF,iF]
      @inbounds for iQ = 1 : NumQuad
        hFLoc[1,iQ] += hFFRef[1,iDoFhFeF,iQ] * hF[ind]
        hFLoc[2,iQ] += hFFRef[2,iDoFhFeF,iQ] * hF[ind]
      end
      @inbounds for iQ = 1 : NumQuadL
        hFLocX[iQ] += hFFRefX[1,iDoFhFeF,iQ] * hF[ind]
        hFLocY[iQ] += hFFRefY[1,iDoFhFeF,iQ] * hF[ind]
      end
    end
    @inbounds for iQ = 1 : NumQuad
      @inbounds for iDoFFeT = 1 : FeT.DoF
        GradLoc[iDoFFeT] -=  Grid.Faces[iF].Orientation * Weights[iQ] * (fTRef[1,iDoFFeT,iQ] * hFLoc[1,iQ] +
          fTRef[2,iDoFFeT,iQ] * hFLoc[2,iQ])
      end
    end
    @inbounds for iQ = 1 : NumQuadL
      @inbounds for iDoFFeT = 1 : FeT.DoF
        GradLoc[iDoFFeT] +=  WeightsL[iQ] * (fTRefX[1,iDoFFeT,iQ] * hFLocX[iQ] +
          fTRefY[1,iDoFFeT,iQ] * hFLocY[iQ])
      end
    end   
    @inbounds for iDoFFeT = 1 : FeT.DoF
      ind = FeT.Glob[iDoFFeT,iF]
      Grad[ind] += GradLoc[iDoFFeT]
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
function GradKinHeight!(backend,FTB,Rhs,FeT::HDivElement,h,hFeF::ScalarElement,u,uFeF::HDivElement,
  Grid,ElemType::Grids.ElementType,QuadOrd,Jacobi)

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
  fTRefXM  = zeros(hFeF.Comp,FeT.DoF,NumQuadL)
  fTRefYM  = zeros(hFeF.Comp,FeT.DoF,NumQuadL)
  hFRefXM  = zeros(hFeF.Comp,hFeF.DoF,NumQuadL)
  hFRefYM  = zeros(hFeF.Comp,hFeF.DoF,NumQuadL)
  uFRefXM  = zeros(uFeF.Comp,uFeF.DoF,NumQuadL)
  uFRefYM  = zeros(uFeF.Comp,uFeF.DoF,NumQuadL)
  fTRefXP  = zeros(hFeF.Comp,FeT.DoF,NumQuadL)
  fTRefYP  = zeros(hFeF.Comp,FeT.DoF,NumQuadL)
  hFRefXP  = zeros(hFeF.Comp,hFeF.DoF,NumQuadL)
  hFRefYP  = zeros(hFeF.Comp,hFeF.DoF,NumQuadL)
  uFRefXP  = zeros(uFeF.Comp,uFeF.DoF,NumQuadL)
  uFRefYP  = zeros(uFeF.Comp,uFeF.DoF,NumQuadL)
  @inbounds for iQ = 1 : NumQuadL
    @inbounds for iComp = 1 : hFeF.Comp
      @inbounds for iDoF = 1 : hFeF.DoF
        hFRefXM[iComp,iDoF,iQ] = hFeF.phi[iDoF,iComp](-1.0,PointsL[iQ])
        hFRefYM[iComp,iDoF,iQ] = hFeF.phi[iDoF,iComp](PointsL[iQ],-1.0)
        hFRefXP[iComp,iDoF,iQ] = hFeF.phi[iDoF,iComp](1.0,PointsL[iQ])
        hFRefYP[iComp,iDoF,iQ] = hFeF.phi[iDoF,iComp](PointsL[iQ],1.0)
      end
    end
    @inbounds for iDoF = 1 : uFeF.DoF
      uFRefXM[1,iDoF,iQ] = uFeF.phi[iDoF,1](-1.0,PointsL[iQ])
      uFRefXM[2,iDoF,iQ] = uFeF.phi[iDoF,2](-1.0,PointsL[iQ])
      uFRefYM[1,iDoF,iQ] = uFeF.phi[iDoF,1](PointsL[iQ],-1.0)
      uFRefYM[2,iDoF,iQ] = uFeF.phi[iDoF,2](PointsL[iQ],-1.0)
      uFRefXP[1,iDoF,iQ] = uFeF.phi[iDoF,1](1.0,PointsL[iQ])
      uFRefXP[2,iDoF,iQ] = uFeF.phi[iDoF,2](1.0,PointsL[iQ])
      uFRefYP[1,iDoF,iQ] = uFeF.phi[iDoF,1](PointsL[iQ],1.0)
      uFRefYP[2,iDoF,iQ] = uFeF.phi[iDoF,2](PointsL[iQ],1.0)
    end
    @inbounds for iDoF = 1 : FeT.DoF
      fTRefXM[1,iDoF,iQ] = FeT.phi[iDoF,1](-1.0,PointsL[iQ])
      fTRefYM[1,iDoF,iQ] = FeT.phi[iDoF,2](PointsL[iQ],-1.0)
      fTRefXP[1,iDoF,iQ] = FeT.phi[iDoF,1](1.0,PointsL[iQ])
      fTRefYP[1,iDoF,iQ] = FeT.phi[iDoF,2](PointsL[iQ],1.0)
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
    @inbounds for iDoF = 1 : uFeF.DoF
      ind = uFeF.Glob[iDoF,iF]  
      uF[iDoF] = u[ind]
    end  
    @inbounds for iDoF = 1 : hFeF.DoF
      ind = hFeF.Glob[iDoF,iF]  
      hF[iDoF] = Grav * h[ind]
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
        hLoc += hFRef[1,iDoF,iQ] * hF[iDoF]
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
        u1X += uFRefXM[1,iDoF,iQ] * uF[iDoF]
        u2X += uFRefXM[2,iDoF,iQ] * uF[iDoF]
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
        u1Y += uFRefYM[1,iDoF,iQ] * uF[iDoF]
        u2Y += uFRefYM[2,iDoF,iQ] * uF[iDoF]
      end
      uLoc1Y = 1 / detDFLoc * (DF[1,1] * u1Y + DF[1,2] * u2Y)
      uLoc2Y = 1 / detDFLoc * (DF[2,1] * u1Y + DF[2,2] * u2Y)
      uLoc3Y = 1 / detDFLoc * (DF[3,1] * u1Y + DF[3,2] * u2Y)
      hLocY = 0.5 * (uLoc1Y * uLoc1Y + uLoc2Y * uLoc2Y + uLoc3Y * uLoc3Y)
      @inbounds for iDoF = 1 : hFeF.DoF
        hLocX += hFRefXM[1,iDoF,iQ] * hF[iDoF]
        hLocY += hFRefYM[1,iDoF,iQ] * hF[iDoF]
      end
      @inbounds for iDoF = 1 : FeT.DoF
#       GradLoc[iDoF] -= WeightsL[iQ] * (Grid.Faces[iF].OrientE[4] * fTRefXM[1,iDoF,iQ] * hLocX +
#       Grid.Faces[iF].OrientE[1] * fTRefYM[1,iDoF,iQ] * hLocY)
        GradLoc[iDoF] -= WeightsL[iQ] * (fTRefXM[1,iDoF,iQ] * hLocX +
        fTRefYM[1,iDoF,iQ] * hLocY)
      end
    end   
    @inbounds for iQ = 1 : NumQuadL
      hLocX = 0.0 
      hLocY = 0.0 
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,1.0,PointsL[iQ],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      u1X = 0.0
      u2X = 0.0
      @inbounds for iDoF = 1 : uFeF.DoF
        u1X += uFRefXP[1,iDoF,iQ] * uF[iDoF]
        u2X += uFRefXP[2,iDoF,iQ] * uF[iDoF]
      end
      uLoc1X = 1 / detDFLoc * (DF[1,1] * u1X + DF[1,2] * u2X)
      uLoc2X = 1 / detDFLoc * (DF[2,1] * u1X + DF[2,2] * u2X)
      uLoc3X = 1 / detDFLoc * (DF[3,1] * u1X + DF[3,2] * u2X)
      hLocX = 0.5 * (uLoc1X * uLoc1X + uLoc2X * uLoc2X + uLoc3X * uLoc3X)
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,PointsL[iQ],1.0,Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      u1Y = 0.0
      u2Y = 0.0
      @inbounds for iDoF = 1 : uFeF.DoF
        u1Y += uFRefYP[1,iDoF,iQ] * uF[iDoF]
        u2Y += uFRefYP[2,iDoF,iQ] * uF[iDoF]
      end
      uLoc1Y = 1 / detDFLoc * (DF[1,1] * u1Y + DF[1,2] * u2Y)
      uLoc2Y = 1 / detDFLoc * (DF[2,1] * u1Y + DF[2,2] * u2Y)
      uLoc3Y = 1 / detDFLoc * (DF[3,1] * u1Y + DF[3,2] * u2Y)
      hLocY = 0.5 * (uLoc1Y * uLoc1Y + uLoc2Y * uLoc2Y + uLoc3Y * uLoc3Y)
      @inbounds for iDoF = 1 : hFeF.DoF
        hLocX += hFRefXP[1,iDoF,iQ] * hF[iDoF]
        hLocY += hFRefYP[1,iDoF,iQ] * hF[iDoF]
      end
      @inbounds for iDoF = 1 : FeT.DoF
#       GradLoc[iDoF] += WeightsL[iQ] * (Grid.Faces[iF].OrientE[2] * fTRefXP[1,iDoF,iQ] * hLocX +
#       Grid.Faces[iF].OrientE[3] * fTRefYP[1,iDoF,iQ] * hLocY)
        GradLoc[iDoF] -= WeightsL[iQ] * (fTRefXP[1,iDoF,iQ] * hLocX +
        fTRefYP[1,iDoF,iQ] * hLocY)
      end
    end   
    @inbounds for iDoF = 1 : FeT.DoF
      ind = FeT.Glob[iDoF,iF]  
      Rhs[ind] += GradLoc[iDoF]
    end  
  end
end

function GradRhs!(backend,FTB,Rhs,FeT::HDivKiteDElement,h,hFeF::ScalarElement,
  Grid,ElemType::Grids.ElementType,QuadOrd,Jacobi)

  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  hFRef  = zeros(hFeF.Comp,hFeF.DoF,NumQuad)
  fTRef  = zeros(hFeF.Comp,FeT.DoF,NumQuad)
  @show "GradRhs! LinKite"

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
  end
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  fTRefXP  = zeros(FeT.Comp,FeT.DoF,NumQuadL)
  fTRefYP  = zeros(FeT.Comp,FeT.DoF,NumQuadL)
  hFRefXP  = zeros(hFeF.Comp,hFeF.DoF,NumQuadL)
  hFRefYP  = zeros(hFeF.Comp,hFeF.DoF,NumQuadL)
  @inbounds for iQ = 1 : NumQuadL
    @inbounds for iComp = 1 : hFeF.Comp
      @inbounds for iDoF = 1 : hFeF.DoF
        hFRefXP[iComp,iDoF,iQ] = hFeF.phi[iDoF,iComp](1.0,PointsL[iQ])
        hFRefYP[iComp,iDoF,iQ] = hFeF.phi[iDoF,iComp](PointsL[iQ],1.0)
      end
    end
    @inbounds for iDoF = 1 : FeT.DoF
      fTRefXP[1,iDoF,iQ] = FeT.phi[iDoF,1](1.0,PointsL[iQ])
      fTRefYP[1,iDoF,iQ] = FeT.phi[iDoF,2](PointsL[iQ],1.0)
    end
  end
  GradLoc = zeros(FeT.DoF)
  hF = zeros(hFeF.DoF)

  @inbounds for iF = 1 : Grid.NumFaces
    GradLoc .= 0
    @inbounds for iDoF = 1 : hFeF.DoF
      ind = hFeF.Glob[iDoF,iF]  
      hF[iDoF] = h[ind]
    end  
    @inbounds for iQ = 1 : NumQuad
      hLoc = 0.0
      @inbounds for iDoF = 1 : hFeF.DoF 
        hLoc += hFRef[1,iDoF,iQ] * hF[iDoF]
      end  
      @inbounds for iDoF = 1 : FeT.DoF
        GradLoc[iDoF] += Weights[iQ] * fTRef[1,iDoF,iQ] * hLoc
      end
    end
    @inbounds for iQ = 1 : NumQuadL
      hLocX = 0.0 
      hLocY = 0.0 
      @inbounds for iDoF = 1 : hFeF.DoF
        hLocX += hFRefXP[1,iDoF,iQ] * hF[iDoF]
        hLocY += hFRefYP[1,iDoF,iQ] * hF[iDoF]
      end
      @inbounds for iDoF = 1 : FeT.DoF
        GradLoc[iDoF] -= WeightsL[iQ] * (fTRefXP[1,iDoF,iQ] * hLocX +
        fTRefYP[1,iDoF,iQ] * hLocY)
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
function DivRhs!(backend,FTB,Div,FeT::ScalarElement,u,uFeF::HDivElement,h,hFeF::ScalarElement,Grid,
  ElemType::Grids.ElementType,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  uFFRef  = zeros(uFeF.Comp,uFeF.DoF,NumQuad)
  hFFRef  = zeros(FeT.Comp,hFeF.DoF,NumQuad)
  GradfTRef  = zeros(uFeF.Comp,hFeF.DoF,NumQuad)


  @show "DivRhs! Kite"
  @inbounds for iQ = 1 : NumQuad
    @inbounds for iComp = 1 : uFeF.Comp
      @inbounds for iD = 1 : uFeF.DoF
        uFFRef[iComp,iD,iQ] = uFeF.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    @inbounds for iComp = 1 : uFeF.Comp
      @inbounds for iD = 1 : hFeF.DoF
        GradfTRef[iComp,iD,iQ] = FeT.Gradphi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : hFeF.DoF
        hFFRef[iComp,iD,iQ] = hFeF.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
  end
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  uFFRefXM  = zeros(FeT.Comp,uFeF.DoF,NumQuadL)
  uFFRefYM  = zeros(FeT.Comp,uFeF.DoF,NumQuadL)
  hFFRefXM  = zeros(FeT.Comp,hFeF.DoF,NumQuadL)
  hFFRefYM  = zeros(FeT.Comp,hFeF.DoF,NumQuadL)
  fTRefXM  = zeros(FeT.Comp,FeT.DoF,NumQuadL)
  fTRefYM  = zeros(FeT.Comp,FeT.DoF,NumQuadL)
  uFFRefXP  = zeros(FeT.Comp,uFeF.DoF,NumQuadL)
  uFFRefYP  = zeros(FeT.Comp,uFeF.DoF,NumQuadL)
  hFFRefXP  = zeros(FeT.Comp,hFeF.DoF,NumQuadL)
  hFFRefYP  = zeros(FeT.Comp,hFeF.DoF,NumQuadL)
  fTRefXP  = zeros(FeT.Comp,FeT.DoF,NumQuadL)
  fTRefYP  = zeros(FeT.Comp,FeT.DoF,NumQuadL)
  @inbounds for iQ = 1 : NumQuadL
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : FeT.DoF
        fTRefXM[iComp,iD,iQ] = FeT.phi[iD,iComp](-1.0,PointsL[iQ])
        fTRefYM[iComp,iD,iQ] = FeT.phi[iD,iComp](PointsL[iQ],-1.0)
        fTRefXP[iComp,iD,iQ] = FeT.phi[iD,iComp](1.0,PointsL[iQ])
        fTRefYP[iComp,iD,iQ] = FeT.phi[iD,iComp](PointsL[iQ],1.0)
      end
    end
    @inbounds for iD = 1 : uFeF.DoF
      uFFRefXM[1,iD,iQ] = uFeF.phi[iD,1](-1.0,PointsL[iQ])
      uFFRefYM[1,iD,iQ] = uFeF.phi[iD,2](PointsL[iQ],-1.0)
      uFFRefXP[1,iD,iQ] = uFeF.phi[iD,1](1.0,PointsL[iQ])
      uFFRefYP[1,iD,iQ] = uFeF.phi[iD,2](PointsL[iQ],1.0)
    end
    @inbounds for iD = 1 : hFeF.DoF
      hFFRefXM[1,iD,iQ] = hFeF.phi[iD,1](-1.0,PointsL[iQ])
      hFFRefYM[1,iD,iQ] = hFeF.phi[iD,1](PointsL[iQ],-1.0)
      hFFRefXP[1,iD,iQ] = hFeF.phi[iD,1](1.0,PointsL[iQ])
      hFFRefYP[1,iD,iQ] = hFeF.phi[iD,1](PointsL[iQ],1.0)
    end
  end

  DivLoc = zeros(FeT.DoF)
  uF = zeros(uFeF.DoF)
  hF = zeros(hFeF.DoF)

  @inbounds for iF = 1 : Grid.NumFaces
    @. DivLoc = 0
    @inbounds for iDoFuFeF = 1 : uFeF.DoF
      ind = uFeF.Glob[iDoFuFeF,iF]  
      uF[iDoFuFeF] = u[ind]
    end  
    @inbounds for iDoFhFeF = 1 : hFeF.DoF
      ind = hFeF.Glob[iDoFhFeF,iF]  
      hF[iDoFhFeF] = h[ind]
    end
    @inbounds for iQ = 1 : NumQuad
      uFLoc1 = 0.0
      uFLoc2 = 0.0
      @inbounds for iDoFuFeF = 1 : uFeF.DoF
        uFLoc1 += uFFRef[1,iDoFuFeF,iQ] * uF[iDoFuFeF]
        uFLoc2 += uFFRef[2,iDoFuFeF,iQ] * uF[iDoFuFeF]
      end  
      hFLoc = 0.0
      @inbounds for iDoFhFeF = 1 : hFeF.DoF
        hFLoc += hFFRef[1,iDoFhFeF,iQ] * hF[iDoFhFeF]
      end  
      @inbounds for iDoFFeT = 1 : FeT.DoF
        DivLoc[iDoFFeT] +=  Weights[iQ] * (GradfTRef[1,iDoFFeT,iQ] * uFLoc1 + GradfTRef[2,iDoFFeT,iQ] *uFLoc2) * hFLoc
      end
    end  
    @inbounds for iQ = 1 : NumQuadL
      uFLocX = 0.0
      uFLocY = 0.0
      @inbounds for iDoFuFeF = 1 : uFeF.DoF
        uFLocX += uFFRefXM[1,iDoFuFeF,iQ] * uF[iDoFuFeF]
        uFLocY += uFFRefYM[1,iDoFuFeF,iQ] * uF[iDoFuFeF]
      end
      hFLocX = 0.0
      hFLocY = 0.0
      @inbounds for iDoFhFeF = 1 : hFeF.DoF
        hFLocX += hFFRefXM[1,iDoFhFeF,iQ] * hF[iDoFhFeF]
        hFLocY += hFFRefYM[1,iDoFhFeF,iQ] * hF[iDoFhFeF]
      end
      @inbounds for iDoFFeT = 1 : FeT.DoF
        DivLoc[iDoFFeT] -= WeightsL[iQ] * (Grid.Faces[iF].OrientE[4] * fTRefXM[1,iDoFFeT,iQ] * hFLocX +
        Grid.Faces[iF].OrientE[1] * fTRefYM[1,iDoFFeT,iQ] * hFLocY)
#       DivLoc[iDoFFeT] -=  WeightsL[iQ] * (fTRefXM[1,iDoFFeT,iQ] * uFLocX * hFLocX +
#         fTRefYM[1,iDoFFeT,iQ] * uFLocY * hFLocY)
      end
    end
    @inbounds for iQ = 1 : NumQuadL
      uFLocX = 0.0
      uFLocY = 0.0
      @inbounds for iDoFuFeF = 1 : uFeF.DoF
        uFLocX += uFFRefXP[1,iDoFuFeF,iQ] * uF[iDoFuFeF]
        uFLocY += uFFRefYP[1,iDoFuFeF,iQ] * uF[iDoFuFeF]
      end
      hFLocX = 0.0
      hFLocY = 0.0
      @inbounds for iDoFhFeF = 1 : hFeF.DoF
        hFLocX += hFFRefXP[1,iDoFhFeF,iQ] * hF[iDoFhFeF]
        hFLocY += hFFRefYP[1,iDoFhFeF,iQ] * hF[iDoFhFeF]
      end
      @inbounds for iDoFFeT = 1 : FeT.DoF
#       DivLoc[iDoFFeT] -= WeightsL[iQ] * (Grid.Faces[iF].OrientE[2] * fTRefXP[1,iDoFFeT,iQ] * hFLocX +
#       Grid.Faces[iF].OrientE[3] * fTRefYP[1,iDoFFeT,iQ] * hFLocY)
#       DivLoc[iDoFFeT] -=  WeightsL[iQ] * (fTRefXP[1,iDoFFeT,iQ] * uFLocX * hFLocX +
#         fTRefYP[1,iDoFFeT,iQ] * uFLocY * hFLocY)
      end
    end
    @inbounds for iDoFFeT = 1 : FeT.DoF
      ind = FeT.Glob[iDoFFeT,iF]
      Div[ind] += DivLoc[iDoFFeT]
    end
  end
end

function DivRhs!(backend,FTB,Div,FeT::ScalarElement,u,uFeF::HDivKiteDElement,Grid,
  ElemType::Grids.ElementType,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  uFFRef  = zeros(uFeF.Comp,uFeF.DoF,NumQuad)
  GradfTRef  = zeros(uFeF.Comp,FeT.DoF,NumQuad)
  @show "DivRhs! LinKite"

  @inbounds for iQ = 1 : NumQuad
    @inbounds for iComp = 1 : uFeF.Comp
      @inbounds for iD = 1 : uFeF.DoF
        uFFRef[iComp,iD,iQ] = uFeF.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    @inbounds for iComp = 1 : uFeF.Comp
      @inbounds for iD = 1 : FeT.DoF
        GradfTRef[iComp,iD,iQ] = FeT.Gradphi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
  end
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  uFFRefXM  = zeros(FeT.Comp,uFeF.DoF,NumQuadL)
  uFFRefYM  = zeros(FeT.Comp,uFeF.DoF,NumQuadL)
  fTRefXM  = zeros(FeT.Comp,FeT.DoF,NumQuadL)
  fTRefYM  = zeros(FeT.Comp,FeT.DoF,NumQuadL)
  @inbounds for iQ = 1 : NumQuadL
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : FeT.DoF
        fTRefXM[iComp,iD,iQ] = FeT.phi[iD,iComp](-1.0,PointsL[iQ])
        fTRefYM[iComp,iD,iQ] = FeT.phi[iD,iComp](PointsL[iQ],-1.0)
      end
    end
    @inbounds for iD = 1 : uFeF.DoF
      uFFRefXM[1,iD,iQ] = uFeF.phi[iD,1](-1.0,PointsL[iQ])
      uFFRefYM[1,iD,iQ] = uFeF.phi[iD,2](PointsL[iQ],-1.0)
    end
  end

  DivLoc = zeros(FeT.DoF)
  uF = zeros(uFeF.DoF)

  @inbounds for iF = 1 : Grid.NumFaces
    @. DivLoc = 0
    @inbounds for iDoF = 1 : uFeF.DoF
      ind = uFeF.Glob[iDoF,iF]  
      uF[iDoF] = u[ind]
    end  
    @inbounds for iQ = 1 : NumQuad
      uFLoc1 = 0.0
      uFLoc2 = 0.0
      @inbounds for iDoF = 1 : uFeF.DoF
        uFLoc1 += uFFRef[1,iDoF,iQ] * uF[iDoF]
        uFLoc2 += uFFRef[2,iDoF,iQ] * uF[iDoF]
      end  
      @inbounds for iDoF = 1 : FeT.DoF
        DivLoc[iDoF] +=  Weights[iQ] * (GradfTRef[1,iDoF,iQ] * uFLoc1 + GradfTRef[2,iDoF,iQ] *uFLoc2)
      end
    end  
    @inbounds for iQ = 1 : NumQuadL
      uFLocX = 0.0
      uFLocY = 0.0
      @inbounds for iDoF = 1 : uFeF.DoF
        uFLocX += uFFRefXM[1,iDoF,iQ] * uF[iDoF]
        uFLocY += uFFRefYM[1,iDoF,iQ] * uF[iDoF]
      end
      @inbounds for iDoF = 1 : FeT.DoF
        DivLoc[iDoF] +=  WeightsL[iQ] * (fTRefXM[1,iDoF,iQ] * uFLocX +
          fTRefYM[1,iDoF,iQ] * uFLocY)
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
  hFLoc = zeros(2,NumQuad)
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
        hFLoc[1,iQ] += hFFRef[1,iDoFhFeF,iQ] * h[ind]
        hFLoc[2,iQ] += hFFRef[2,iDoFhFeF,iQ] * h[ind]
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

function CurlVel!(q,FeT::CGKitePrimalStruct,u,uFe::HDivKiteDElement,QuadOrd,ElemType,Grid,Jacobi)
#
#
# int q*v dx = int Curl u * v dx = - int u * rot v dx 
#
#

  @show "CurlKite Primal"
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
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  uFRefXM  = zeros(uFe.Comp,uFe.DoF,NumQuadL)
  uFRefYM  = zeros(uFe.Comp,uFe.DoF,NumQuadL)
  fTRefXM  = zeros(FeT.Comp,FeT.DoF,NumQuadL)
  fTRefYM  = zeros(FeT.Comp,FeT.DoF,NumQuadL)
  @inbounds for iQ = 1 : NumQuadL
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : FeT.DoF
        fTRefXM[iComp,iD,iQ] = FeT.phi[iD,iComp](-1.0,PointsL[iQ])
        fTRefYM[iComp,iD,iQ] = FeT.phi[iD,iComp](PointsL[iQ],-1.0)
      end
    end
    @inbounds for iD = 1 : uFe.DoF
      uFRefXM[1,iD,iQ] = uFe.phi[iD,1](-1.0,PointsL[iQ])
      uFRefXM[2,iD,iQ] = uFe.phi[iD,2](-1.0,PointsL[iQ])
      uFRefYM[1,iD,iQ] = uFe.phi[iD,1](PointsL[iQ],-1.0)
      uFRefYM[2,iD,iQ] = uFe.phi[iD,2](PointsL[iQ],-1.0)
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
        qLoc[iDoF] -= Weights[iQ] * (uFLoc[1] * RotqRef[1,iDoF,iQ] +
          uFLoc[2] * RotqRef[2,iDoF,iQ])  
      end  
    end  
    @inbounds for iQ = 1 : NumQuadL
      @. uFLoc = 0  
      @inbounds for iDoF = 1 : uFe.DoF
        uFLoc[1] += uFRefXM[1,iDoF,iQ] * uLoc[iDoF]
        uFLoc[2] += uFRefXM[2,iDoF,iQ] * uLoc[iDoF]
      end
      Jacobi(DF,detDF,pinvDF,X,ElemType,-1.0,PointsL[iQ],Grid.Faces[iF], Grid)
      uuFLoc[1] = (DF[1,1] * uFLoc[1] + DF[1,2] * uFLoc[2]) / detDF[1]
      uuFLoc[2] = (DF[2,1] * uFLoc[1] + DF[2,2] * uFLoc[2]) / detDF[1]
      uuFLoc[3] = (DF[3,1] * uFLoc[1] + DF[3,2] * uFLoc[2]) / detDF[1]
      uLocX = DF[1,2] * uuFLoc[1] + DF[2,2] * uuFLoc[2] + DF[3,2] * uuFLoc[3]

      @. uFLoc = 0  
      @inbounds for iDoF = 1 : uFe.DoF
        uFLoc[1] += uFRefYM[1,iDoF,iQ] * uLoc[iDoF]
        uFLoc[2] += uFRefYM[2,iDoF,iQ] * uLoc[iDoF]
      end
      Jacobi(DF,detDF,pinvDF,X,ElemType,PointsL[iQ],-1.0,Grid.Faces[iF], Grid)
      uuFLoc[1] = (DF[1,1] * uFLoc[1] + DF[1,2] * uFLoc[2]) / detDF[1]
      uuFLoc[2] = (DF[2,1] * uFLoc[1] + DF[2,2] * uFLoc[2]) / detDF[1]
      uuFLoc[3] = (DF[3,1] * uFLoc[1] + DF[3,2] * uFLoc[2]) / detDF[1]
      uLocY = DF[1,1] * uuFLoc[1] + DF[2,1] * uuFLoc[2] + DF[3,1] * uuFLoc[3]
      @inbounds for iDoF = 1 : FeT.DoF
        qLoc[iDoF] += WeightsL[iQ] * (-fTRefXM[1,iDoF,iQ] * uLocX +
         fTRefYM[1,iDoF,iQ] * uLocY)
      end
    end
    @inbounds for iDoF = 1 : FeT.DoF
      ind = FeT.Glob[iDoF,iF]
      q[ind] += qLoc[iDoF]
    end
  end    
  ldiv!(FeT.LUM,q)
end

function CurlVel!(q,FeT::CGKiteDualStruct,u,uFe::HDivKiteDElement,QuadOrd,ElemType,Grid,Jacobi)
#
#
# int q*v dx = int Curl u * v dx = - int u * rot v dx 
#
#
  @show "CurlKite Dual"

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
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  uFRefXM  = zeros(uFe.Comp,uFe.DoF,NumQuadL)
  uFRefYM  = zeros(uFe.Comp,uFe.DoF,NumQuadL)
  fTRefXM  = zeros(FeT.Comp,FeT.DoF,NumQuadL)
  fTRefYM  = zeros(FeT.Comp,FeT.DoF,NumQuadL)
  uFRefXP  = zeros(uFe.Comp,uFe.DoF,NumQuadL)
  uFRefYP  = zeros(uFe.Comp,uFe.DoF,NumQuadL)
  fTRefXP  = zeros(FeT.Comp,FeT.DoF,NumQuadL)
  fTRefYP  = zeros(FeT.Comp,FeT.DoF,NumQuadL)
  @inbounds for iQ = 1 : NumQuadL
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : FeT.DoF
        fTRefXM[iComp,iD,iQ] = FeT.phi[iD,iComp](-1.0,PointsL[iQ])
        fTRefXM[iComp,iD,iQ] = FeT.phi[iD,iComp](-1.0,PointsL[iQ])
        fTRefYP[iComp,iD,iQ] = FeT.phi[iD,iComp](PointsL[iQ],1.0)
        fTRefYP[iComp,iD,iQ] = FeT.phi[iD,iComp](PointsL[iQ],1.0)
      end
    end
    @inbounds for iD = 1 : uFe.DoF
      uFRefXM[1,iD,iQ] = uFe.phi[iD,1](-1.0,PointsL[iQ])
      uFRefXM[2,iD,iQ] = uFe.phi[iD,2](-1.0,PointsL[iQ])
      uFRefYM[1,iD,iQ] = uFe.phi[iD,1](PointsL[iQ],-1.0)
      uFRefYM[2,iD,iQ] = uFe.phi[iD,2](PointsL[iQ],-1.0)
      uFRefXP[1,iD,iQ] = uFe.phi[iD,1](1.0,PointsL[iQ])
      uFRefXP[2,iD,iQ] = uFe.phi[iD,2](1.0,PointsL[iQ])
      uFRefYP[1,iD,iQ] = uFe.phi[iD,1](PointsL[iQ],1.0)
      uFRefYP[2,iD,iQ] = uFe.phi[iD,2](PointsL[iQ],1.0)
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
    @inbounds for iQ = 1 : NumQuadL
      @. uFLoc = 0
      @inbounds for iDoF = 1 : uFe.DoF
        uFLoc[1] += uFRefXP[1,iDoF,iQ] * uLoc[iDoF]
        uFLoc[2] += uFRefXP[2,iDoF,iQ] * uLoc[iDoF]
      end
      Jacobi(DF,detDF,pinvDF,X,ElemType,1.0,PointsL[iQ],Grid.Faces[iF], Grid)
      uuFLoc[1] = (DF[1,1] * uFLoc[1] + DF[1,2] * uFLoc[2]) / detDF[1]
      uuFLoc[2] = (DF[2,1] * uFLoc[1] + DF[2,2] * uFLoc[2]) / detDF[1]
      uuFLoc[3] = (DF[3,1] * uFLoc[1] + DF[3,2] * uFLoc[2]) / detDF[1]
      uLocX = DF[1,2] * uuFLoc[1] + DF[2,2] * uuFLoc[2] + DF[3,2] * uuFLoc[3]
      @. uFLoc = 0
      @inbounds for iDoF = 1 : uFe.DoF
        uFLoc[1] += uFRefYP[1,iDoF,iQ] * uLoc[iDoF]
        uFLoc[2] += uFRefYP[2,iDoF,iQ] * uLoc[iDoF]
      end
      Jacobi(DF,detDF,pinvDF,X,ElemType,PointsL[iQ],1.0,Grid.Faces[iF], Grid)
      uuFLoc[1] = (DF[1,1] * uFLoc[1] + DF[1,2] * uFLoc[2]) / detDF[1]
      uuFLoc[2] = (DF[2,1] * uFLoc[1] + DF[2,2] * uFLoc[2]) / detDF[1]
      uuFLoc[3] = (DF[3,1] * uFLoc[1] + DF[3,2] * uFLoc[2]) / detDF[1]
      uLocY = DF[1,1] * uuFLoc[1] + DF[2,1] * uuFLoc[2] + DF[3,1] * uuFLoc[3]
      @inbounds for iDoF = 1 : FeT.DoF
#        qLoc[iDoF] -= WeightsL[iQ] * (Grid.Faces[iF].OrientE[2]*fTRefXP[1,iDoF,iQ] * uLocX -
#        Grid.Faces[iF].OrientE[3]*fTRefYP[1,iDoF,iQ] * uLocY)
      end
    end
    @inbounds for iQ = 1 : NumQuadL
      @. uFLoc = 0
      @inbounds for iDoF = 1 : uFe.DoF
        uFLoc[1] += uFRefXM[1,iDoF,iQ] * uLoc[iDoF]
        uFLoc[2] += uFRefXM[2,iDoF,iQ] * uLoc[iDoF]
      end
      Jacobi(DF,detDF,pinvDF,X,ElemType,-1.0,PointsL[iQ],Grid.Faces[iF], Grid)
      uuFLoc[1] = (DF[1,1] * uFLoc[1] + DF[1,2] * uFLoc[2]) / detDF[1]
      uuFLoc[2] = (DF[2,1] * uFLoc[1] + DF[2,2] * uFLoc[2]) / detDF[1]
      uuFLoc[3] = (DF[3,1] * uFLoc[1] + DF[3,2] * uFLoc[2]) / detDF[1]
      uLocX = DF[1,2] * uuFLoc[1] + DF[2,2] * uuFLoc[2] + DF[3,2] * uuFLoc[3]
      @. uFLoc = 0
      @inbounds for iDoF = 1 : uFe.DoF
        uFLoc[1] += uFRefYM[1,iDoF,iQ] * uLoc[iDoF]
        uFLoc[2] += uFRefYM[2,iDoF,iQ] * uLoc[iDoF]
      end
      Jacobi(DF,detDF,pinvDF,X,ElemType,PointsL[iQ],-1.0,Grid.Faces[iF], Grid)
      uuFLoc[1] = (DF[1,1] * uFLoc[1] + DF[1,2] * uFLoc[2]) / detDF[1]
      uuFLoc[2] = (DF[2,1] * uFLoc[1] + DF[2,2] * uFLoc[2]) / detDF[1]
      uuFLoc[3] = (DF[3,1] * uFLoc[1] + DF[3,2] * uFLoc[2]) / detDF[1]
      uLocY = DF[1,1] * uuFLoc[1] + DF[2,1] * uuFLoc[2] + DF[3,1] * uuFLoc[3]
      @inbounds for iDoF = 1 : FeT.DoF
#       qLoc[iDoF] += WeightsL[iQ] * (-Grid.Faces[iF].OrientE[4]*fTRefXM[1,iDoF,iQ] * uLocX +
#        Grid.Faces[iF].OrientE[1]*fTRefYM[1,iDoF,iQ] * uLocY)
      end
    end
    @inbounds for iDoF = 1 : FeT.DoF
      ind = FeT.Glob[iDoF,iF]
      q[ind] += qLoc[iDoF]
    end
  end    
  ldiv!(FeT.LUM,q)
end

