mutable struct DG0Struct{FT<:AbstractFloat,
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

#DG0 Quad

function DG0Struct{FT}(::Grids.Quad,backend,Grid) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  Type = Grids.Quad()
  DoF = 1
  Comp = 1
  @polyvar x1 x2
  phi = Array{Polynomial,2}(undef,DoF,Comp)
  phi[1,1] = 1.0 + 0.0*x1 + 0.0*x2
  
    
  Glob = KernelAbstractions.zeros(backend,Int,DoF,Grid.NumFaces)
  GlobCPU = zeros(Int,DoF,Grid.NumFaces)
  NumG = Grid.NumFaces
  NumI = Grid.NumFaces
  for iF = 1 : Grid.NumFaces
    GlobCPU[1,iF] = Grid.Faces[iF].F
  end
  copyto!(Glob,GlobCPU)
  M = spzeros(0,0)
  return DG0Struct{FT,
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

#DG0 Tri

function DG0Struct{FT}(Type::Grids.Tri,backend,Grid) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  DoF = 1
  Comp = 1
  phi = Array{Polynomial,2}(undef,DoF,Comp)
  @polyvar x1 x2
  phi[1,1] = 1.0 + 0.0*x1 + 0.0*x2
        
  Glob = KernelAbstractions.zeros(backend,Int,DoF,Grid.NumFaces)
  GlobCPU = zeros(Int,DoF,Grid.NumFaces)
  NumG = Grid.NumFaces
  NumI = Grid.NumFaces
  for iF = 1 : Grid.NumFaces
    GlobCPU[1,iF] = Grid.Faces[iF].F
  end
  copyto!(Glob,GlobCPU)
  M = spzeros(0,0)
  return DG0Struct{FT,
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

mutable struct DG1Struct{FT<:AbstractFloat,
                        IT2<:AbstractArray} <: ScalarElement
  Glob::IT2
  DoF::Int
  Comp::Int                      
  phi::Array{Polynomial,2} 
  Divphi::Array{Polynomial,2}                      
  NumG::Int
  NumI::Int
  Type::Grids.ElementType
  M::AbstractSparseMatrix
end

#DG1 Quad

function DG1Struct{FT}(::Grids.Quad,backend,Grid) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  Type = Grids.Quad()
  DoF = 4
  Comp = 2
  @polyvar x1 x2
  phi = Array{Polynomial,2}(undef,DoF,Comp)
  Divphi = Array{Polynomial,2}(undef,DoF,1)
  
  phi[1,1] = (1.0-1.0*x1)+0.0*x2
  phi[1,2] = (1.0-1.0*x2)+0.0*x1
  
  phi[2,1] = 1.0*x1+0.0*x2
  phi[2,2] = (1-1.0*x2)+0.0*x1

  phi[3,1] = 1.0*x1+0.0*x2
  phi[3,2] = 1.0*x2+0.0*x1

  phi[4,1] = (1-1.0*x1)+0.0*x2
  phi[4,2] = 1.0*x2+0.0*x1

  for i = 1 : DoF
    Divphi[i,1] = differentiate(phi[i,1],x1) + differentiate(phi[i,2],x2)
  end
   
  Glob = KernelAbstractions.zeros(backend,Int,DoF,Grid.NumFaces)
  GlobCPU = zeros(Int,DoF,Grid.NumFaces)
  NumG = Grid.NumFaces
  NumI = Grid.NumFaces
  for iF = 1 : Grid.NumFaces
    GlobCPU[1,iF] = Grid.Faces[iF].F
  end
  copyto!(Glob,GlobCPU)
  M = spzeros(0,0)
  return DG1Struct{FT,
                  typeof(Glob)}( 
    Glob,
    DoF,
    Comp,
    phi,
    Divphi,                      
    NumG,
    NumI,
    Type,
    M,
      )
end

#DG1 TRi

function DG1Struct{FT}(Type::Grids.Tri,backend,Grid) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  DoF = 3
  Comp = 2
  phi = Array{Polynomial,2}(undef,DoF,Comp)
  @polyvar x1 x2
  phi[1,1] = (1.0-1.0*x1)+0.0*x2
  phi[1,2] = -1.0*x2+0.0*x1
  
  phi[2,1] = 1.0*x1+0.0*x2
  phi[2,2] = 0.0*x2+0.0*x1

  phi[3,1] = 0.0*x1+0.0*x2
  phi[3,2] = 1.0*x2+0.0*x1

  for i = 1 : DoF
    Divphi[i,1] = differentiate(phi[i,1],x1) + differentiate(phi[i,2],x2)
  end
        
  Glob = KernelAbstractions.zeros(backend,Int,DoF,Grid.NumFaces)
  GlobCPU = zeros(Int,DoF,Grid.NumFaces)
  NumG = Grid.NumFaces
  NumI = Grid.NumFaces
  for iF = 1 : Grid.NumFaces
    GlobCPU[1,iF] = Grid.Faces[iF].F
  end
  copyto!(Glob,GlobCPU)
  M = spzeros(0,0)
  return DG1Struct{FT,
                  typeof(Glob)}( 
    Glob,
    DoF,
    Comp,
    phi,
    Divphi,                      
    NumG,
    NumI,
    Type,
    M,
      )
end


