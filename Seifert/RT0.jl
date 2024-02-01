abstract type FiniteElement end

#abstract type ElementType end
#struct Tri <: ElementType end
#struct Quad <: ElementType end
mutable struct RT0Struct{FT<:AbstractFloat,
                        IT2<:AbstractArray} <: FiniteElement
  Glob::IT2
  DoF::Int
  Comp::Int                      
  phi::Array{Polynomial,3}                       
  NumG::Int
  NumI::Int
  Type::Grids.ElementType
end

function RT0Struct{FT}(::Grids.Quad,backend,Grid) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  Type = Grids.Quad()
  DoF = 4
  Comp = 2
  phi = Array{Polynomial,3}(undef,DoF,Comp,2)
  phi[1,1,1] = Polynomial([0])
  phi[1,1,2] = Polynomial([0])
  phi[1,2,1] = Polynomial([1])
  phi[1,2,2] = Lagrange([-1,1], [1,0])

  phi[2,1,1] = Lagrange([-1,1], [0,1])
  phi[2,1,2] = Polynomial([1])
  phi[2,2,1] = Polynomial([0])
  phi[2,2,2] = Polynomial([0])

  phi[3,1,1] = Polynomial([0])
  phi[3,1,2] = Polynomial([0])
  phi[3,2,1] = Polynomial([1])
  phi[3,2,2] = Lagrange([-1,1], [0,1])

  phi[4,1,1] = Lagrange([-1,1], [1,0])
  phi[4,1,2] = Polynomial([1])
  phi[4,2,1] = Polynomial([0])
  phi[4,2,2] = Polynomial([0])

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
    phi = Array{Polynomial,3}(undef,DoF,Comp,2)
    phi[1,1,1] = Polynomial([0.5,0.5]) 
    phi[1,1,2] = Polynomial([1])
    phi[1,2,1] = Polynomial([1])
    phi[1,2,2] = Polynomial([0.5,-0.5])

    phi[2,1,1] = Polynomial([0.5,0.5])
    phi[2,1,2] = Polynomial([1])
    phi[2,2,1] = Polynomial([1])
    phi[2,2,2] = Polynomial([0.5,0.5])

    phi[3,1,1] = Polynomial([-0.5,0.5])
    phi[3,1,2] = Polynomial([1])
    phi[3,2,1] = Polynomial([1])
    phi[3,2,2] = Polynomial([-0.5,-0.5])
    
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