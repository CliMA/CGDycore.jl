mutable struct CG1Struct{FT<:AbstractFloat,
                        IT2<:AbstractArray} <: ScalarElement
  Glob::IT2
  DoF::Int
  Comp::Int                      
  phi::Array{Polynomial,2} 
  Gradphi::Array{Polynomial,2} 
  NumG::Int
  NumI::Int
  Type::Grids.ElementType
  M::AbstractSparseMatrix
  LUM::SparseArrays.UMFPACK.UmfpackLU{Float64, Int64}
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
  Gradphi = Array{Polynomial,2}(undef,DoF,2)
  for i = 1 : DoF
    Gradphi[i,1] = differentiate(phi[i,1],x1)
    Gradphi[i,2] = differentiate(phi[i,1],x2)
  end
  

  Glob = KernelAbstractions.zeros(backend,Int,DoF,Grid.NumNodes)
  GlobCPU = zeros(Int,DoF,Grid.NumNodes)
  NumG = Grid.NumNodes
  NumI = Grid.NumNodes
  for iF = 1 : Grid.NumFaces
    for i = 1 : length(Grid.Faces[iF].N)
      iN = Grid.Faces[iF].N[i]
      GlobCPU[i,iF] = Grid.Nodes[iN].N
    end  
  end  
  copyto!(Glob,GlobCPU)
  M = sparse([1],[1],[1.0])
  LUM = lu(M)
  return CG1Struct{FT,
                  typeof(Glob)}( 
    Glob,
    DoF,
    Comp,
    phi,
    Gradphi,
    NumG,
    NumI,
    Type,
    M,
    LUM,
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
  Gradphi = Array{Polynomial,2}(undef,DoF,2)
  for i = 1 : DoF
    Gradphi[i,1] = differentiate(phi[i,1],x1)
    Gradphi[i,2] = differentiate(phi[i,1],x2)
  end


  Glob = KernelAbstractions.zeros(backend,Int,DoF,Grid.NumNodes)
  GlobCPU = zeros(Int,DoF,Grid.NumNodes)
  NumG = Grid.NumNodes
  NumI = Grid.NumNodes
  for iF = 1 : Grid.NumFaces
    for i = 1 : length(Grid.Faces[iF].N)
      iN = Grid.Faces[iF].N[i]
      GlobCPU[i,iF] = Grid.Nodes[iE].N
    end  
  end  
  copyto!(Glob,GlobCPU)
  M = sparse([1],[1],[1.0])
  LUM = lu(M)
  return CG1Struct{FT,
                  typeof(Glob)}( 
    Glob,
    DoF,
    Comp,
    phi,
    Gradphi,
    NumG,
    NumI,
    Type,
    M,
    LUM,
      )
end
