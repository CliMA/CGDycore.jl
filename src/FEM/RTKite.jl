mutable struct RTKiteDualHDiv{FT<:AbstractFloat,
                        IT2<:AbstractArray} <: HDivKiteDElement
  Order::Int                      
  Glob::IT2
  DoF::Int
  DoFN::Int
  DoFE::Int
  Comp::Int                      
  phi::Array{Polynomial,2}                       
  Divphi::Array{Polynomial,2}                       
  NumG::Int
  NumI::Int
  Type::Grids.ElementType
  M::AbstractSparseMatrix
  LUM
# LUM::SparseArrays.UMFPACK.UmfpackLU{Float64, Int64}
end

function RTKiteDualHDiv{FT}(k,::Grids.Quad,backend,Grid) where FT<:AbstractFloat
  DoF, DoFE, DoFF, phi, Divphi = FEM.ConstructRT(k,Grids.Quad())
  DoFN = 0

  Order = k

  Type = Grids.QuadDual()
  Type = Grids.Quad()
  Comp = 2

  NumFaces = Grid.NumFaces
  NumEdges = Grid.NumEdges
  NumNodes = Grid.NumNodes
  Faces = Grid.Faces
  Edges = Grid.Edges
  Nodes = Grid.Nodes

  Glob = KernelAbstractions.zeros(backend,Int,DoF,NumFaces)
  GlobCPU = zeros(Int,DoF,NumFaces)

  EdgesD = zeros(Int,NumEdges)

  NumEdgesD = 0
  @inbounds for iE = 1 : NumEdges
    if Edges[iE].Type == "E1" || Edges[iE].Type == "E2"
      NumEdgesD += 1
      EdgesD[iE] = NumEdgesD
    end
  end

  @inbounds for iF = 1 : NumFaces
    iG = 1
    iGF = 1
    iE = EdgesD[Faces[iF].E[1]]
    @inbounds for i = 1 : DoFE
      GlobCPU[iG,iF] = DoFE * (iE - 1) + i 
      iG += 1
    end
    @inbounds for i = 1 : DoFE
      GlobCPU[iG,iF] = DoFE * NumEdgesD + (2 * DoFE + DoFF) * (iF - 1) + iGF
      iG += 1
      iGF += 1
    end
    @inbounds for i = 1 : DoFE
      GlobCPU[iG,iF] = DoFE * NumEdgesD + (2 * DoFE + DoFF) * (iF - 1) + iGF
      iG += 1
      iGF += 1
    end
    iE = EdgesD[Faces[iF].E[4]]
    @inbounds for i = 1 : DoFE
      GlobCPU[iG,iF] = DoFE * (iE - 1) + i 
      iG += 1
    end
    @inbounds for i = 1 : DoFF
      GlobCPU[iG,iF] = DoFE * NumEdgesD + (2 * DoFE + DoFF) * (iF - 1) + iGF
      iG += 1
      iGF += 1
    end
  end

  NumG = DoFE * NumEdgesD + (2 * DoFE + DoFF) * NumFaces
  NumI = NumG


  copyto!(Glob,GlobCPU)
  M = sparse([1],[1],[1.0])
  return RTKiteDualHDiv{FT,
                  typeof(Glob)}( 
    Order,              
    Glob,
    DoF,
    DoFN,
    DoFE,
    Comp,
    phi,                      
    Divphi,                      
    NumG,
    NumI,
    Type,
    M,
    nothing,
      )
end
