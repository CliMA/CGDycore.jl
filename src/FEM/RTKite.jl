mutable struct RT1KiteDualHDiv{FT<:AbstractFloat,
                        IT2<:AbstractArray} <: HDivKiteDElement
  Order::Int                      
  Glob::IT2
  DoF::Int
  Comp::Int                      
  phi::Array{Polynomial,2}                       
  Divphi::Array{Polynomial,2}                       
  NumG::Int
  NumI::Int
  Type::Grids.ElementType
  M::AbstractSparseMatrix
  LUM::SparseArrays.UMFPACK.UmfpackLU{Float64, Int64}
end

function RT1KiteDualHDiv{FT}(k,::Grids.Quad,backend,Grid) where FT<:AbstractFloat
  DoF, DoFE, DoFF, phi, Divphi = FEM.ConstructRT(k,Grids.Quad())

  @show DoF, DoFE, DoFF
  for iDoF = 1 : DoF
    @show phi[iDoF,:]
  end

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

  NumEdgesD = 0
  NumEdgesP = 0
  EdgesD = zeros(Int,NumEdges)
  for iE = 1 : NumEdges
    if Edges[iE].Type == "E1" || Edges[iE].Type == "E2"
      NumEdgesD += 1
      EdgesD[iE] = NumEdgesD
    else
      NumEdgesP += 1
    end
  end

  Glob = KernelAbstractions.zeros(backend,Int,DoF,NumFaces)
  GlobCPU = zeros(Int,DoF,NumFaces)
  for iF = 1 : NumFaces
    iG = 1
    iGF = 1
    iE = EdgesD[Faces[iF].E[1]]
    for i = 1 : DoFE
      GlobCPU[iG,iF] = DoFE * (iE - 1) + i 
      iG += 1
    end
    for i = 1 : DoFE
      GlobCPU[iG,iF] = DoFE * NumEdgesD + (2 * DoFE + DoFF) * (iF - 1) + iGF
      iG += 1
      iGF += 1
    end
    for i = 1 : DoFE
      GlobCPU[iG,iF] = DoFE * NumEdgesD + (2 * DoFE + DoFF) * (iF - 1) + iGF
      iG += 1
      iGF += 1
    end
    iE = EdgesD[Faces[iF].E[4]]
    for i = 1 : DoFE
      GlobCPU[iG,iF] = DoFE * (iE - 1) + i 
      iG += 1
    end
    for i = 1 : DoFF
      GlobCPU[iG,iF] = DoFE * NumEdgesD + (2 * DoFE + DoFF) * (iF - 1) + iGF
      iG += 1
      iGF += 1
    end
  end



  NumG = DoFE * NumEdgesD + (2 * DoFE + DoFF) * NumFaces
  NumI = NumG
  copyto!(Glob,GlobCPU)
  M = sparse([1],[1],[1.0])
  LUM = lu(M)
  return RT1KiteDualHDiv{FT,
                  typeof(Glob)}( 
    Order,              
    Glob,
    DoF,
    Comp,
    phi,                      
    Divphi,                      
    NumG,
    NumI,
    Type,
    M,
    LUM,
      )
end
