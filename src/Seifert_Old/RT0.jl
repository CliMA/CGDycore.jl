mutable struct RT0Struct{FT<:AbstractFloat,
                        IT2<:AbstractArray} <: HDivElement
  Glob::IT2
  DoF::Int
  Comp::Int                      
  phi::Array{Polynomial,2}                       
  NumG::Int
  NumI::Int
  Type::Grids.ElementType
end

function RT0Struct{FT}(::Grids.Quad,backend,Grid) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  Type = Grids.Quad()
  DoF = 4
  Comp = 2
  @polyvar x1 x2
  phi = Array{Polynomial,2}(undef,DoF,Comp)
  phi[1,1] = 0.0*x1 + 0.0*x2
  phi[1,2] = 0.5 - 0.5*x2  + 0.0*x1

  phi[2,1] = 0.5 + 0.5*x1 + 0.0*x2
  phi[2,2] = 0*x1 + 0.0*x2

  phi[3,1] = 0*x1 + 0.0*x2
  phi[3,2] = 0.5 + 0.5*x2 + 0.0*x1

  phi[4,1] = 0.5 - 0.5*x1 + 0.0*x2
  phi[4,2] = 0*x1 + 0.0*x2

  Glob = KernelAbstractions.zeros(backend,Int,DoF,Grid.NumFaces)
  GlobCPU = zeros(Int,DoF,Grid.NumFaces)
  NumG = Grid.NumEdgesI + Grid.NumEdgesB
  NumI = Grid.NumEdgesI
  for iF = 1 : Grid.NumFaces
    for i = 1 : length(Grid.Faces[iF].E)
      iE = Grid.Faces[iF].E[i]
      GlobCPU[i,iF] = Grid.Edges[iE].E
    end
  end
  copyto!(Glob,GlobCPU)
  return RT0Struct{FT,
                  typeof(Glob)}( 
  Glob,
  DoF,
  Comp,
  phi,                      
  NumG,
  NumI,
  Type,
    )
  end

function RT0Struct{FT}(type::Grids.Tri,backend,Grid) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  Type = Grids.Quad()
  DoF = 3
  Comp = 2
  phi = Array{Polynomial,2}(undef,DoF,Comp)
  @polyvar x1 x2
  phi[1,1] = 0.5 + 0.5*x1 + 0.0*x2
  phi[1,2] = 0.5 - 0.5*x2 + 0.0*x1

  phi[2,1] = 0.5 + 0.5*x1 + 0.0*x2
  phi[2,2] = 0.5 + 0.5*x2 + 0.0*x1

  phi[3,1] = -0.5 + 0.5*x1 + 0.0*x2
  phi[3,2] = -0.5 - 0.5*x2 + 0.0*x1

  Glob = KernelAbstractions.zeros(backend,Int,DoF,Grid.NumFaces)
  GlobCPU = zeros(Int,DoF,Grid.NumFaces)
  NumG = Grid.NumEdgesI + Grid.NumEdgesB
  NumI = Grid.NumEdgesI
  for iF = 1 : Grid.NumFaces
    for i = 1 : length(Grid.Faces[iF].E)
      iE = Grid.Faces[iF].E[i]
      GlobCPU[i,iF] = Grid.Edges[iE].E
    end
  end
  copyto!(Glob,GlobCPU)
  return RT0Struct{FT,
                  typeof(Glob)}( 
    Glob,
    DoF,
    Comp,
    phi,                      
    NumG,
    NumI,
    Type,
      )
end
