mutable struct CG1Struct{FT<:AbstractFloat,
                        IT2<:AbstractArray} <: ScalarElement
  Glob::IT2
  DoF::Int
  Comp::Int                      
  phi::Array{Polynomial,2} 
  NumG::Int
  NumI::Int
  Type::Grids.ElementType
  M::AbstractSparseMatrix
end

#CG1 Quad
#=
4__3
|  |
1__2
=#

function CG1Struct{FT}(::Grids.Quad,backend,Grid) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  Type = Grids.Quad()
  DoF = 4
  Comp = 1
  nu = Array{Polynomial,2}(undef,DoF,Comp)
  phi = Array{Polynomial,2}(undef,DoF,Comp)
  @polyvar x1 x2 ksi1 ksi2
  nu[1,1] = 1.0*ksi1*ksi2 - 1.0*ksi1 - 1.0*ksi2 + 1.0
  nu[2,1] = 1.0*ksi1 * (1.0 - 1.0*ksi2)
  nu[3,1] = 1.0*ksi2 * (1.0 - 1.0*ksi1)
  nu[4,1] = 1.0*ksi1*ksi2 + 0.0
  
  for s = 1 : DoF
    for t = 1 : 1
      phi[s,t] = subs(nu[s,t], ksi1 => (x1+1)/2, ksi2 => (x2+1)/2)
    end
  end
  

  Glob = KernelAbstractions.zeros(backend,Int,DoF,Grid.NumNodes)
  GlobCPU = zeros(Int,DoF,Grid.NumNodes)
  NumG = 4 * Grid.NumNodes
  NumI = 4 * Grid.NumNodes
  for iF = 1 : Grid.NumNodes
    GlobCPU[1,iF] = 4 * Grid.Nodes[iF].N - 3
    GlobCPU[2,iF] = 4 * Grid.Nodes[iF].N - 2
    GlobCPU[3,iF] = 4 * Grid.Nodes[iF].N - 1
    GlobCPU[4,iF] = 4 * Grid.Nodes[iF].N
  end
  copyto!(Glob,GlobCPU)
  M = spzeros(0,0)
  return CG1Struct{FT,
                  typeof(Glob)}( 
    Glob,
    DoF,
    Comp,
    phi,
    NumG,
    NumI,
    Type,
    M,
      )
end

#CG1 Tri
#=
3
| \
1__2
=#

function CG1Struct{FT}(Type::Grids.Tri,backend,Grid) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  DoF = 3
  Comp = 1
  nu = Array{Polynomial,2}(undef,DoF,Comp)
  phi = Array{Polynomial,2}(undef,DoF,Comp)
  @polyvar x1 x2 ksi1 ksi2
  nu[1,1] = -1.0*ksi1 - 1.0*ksi2 + 1.0
  nu[2,1] = 1.0*ksi1 + 0.0*ksi2 + 0.0
  nu[3,1] = 0.0*ksi1 + 1.0*ksi2 + 0.0
  
  for s = 1 : DoF
    for t = 1 : 1
      phi[s,t] = subs(nu[s,t], ksi1 => (x1+1)/2, ksi2 => (x2+1)/2)
    end
  end


  Glob = KernelAbstractions.zeros(backend,Int,DoF,Grid.NumNodes)
  GlobCPU = zeros(Int,DoF,Grid.NumNodes)
  NumG = 3 * Grid.NumNodes
  NumI = 3 * Grid.NumNodes
  for iF = 1 : Grid.NumNodes
    GlobCPU[1,iF] = 3 * Nodes[iF].N - 2
    GlobCPU[2,iF] = 3 * Nodes[iF].N - 1
    GlobCPU[3,iF] = 3 * Nodes[iF].N
  end
  copyto!(Glob,GlobCPU)
  M = spzeros(0,0)
  return CG1Struct{FT,
                  typeof(Glob)}( 
    Glob,
    DoF,
    Comp,
    phi,
    NumG,
    NumI,
    Type,
    M,
      )
end