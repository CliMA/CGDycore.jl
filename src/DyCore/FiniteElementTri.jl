abstract type FiniteElement end

#abstract type ElementType end
#struct Tri <: ElementType end
#struct Quad <: ElementType end

mutable struct DG0Struct{FT<:AbstractFloat,
                        IT2<:AbstractArray} <: FiniteElement
  Glob::IT2
  Stencil::IT2
  NumG::Int
  NumI::Int
  Type::Grids.ElementType
end

function DG0Struct{FT}(backend,::Grids.Tri,Grid) where FT<:AbstractFloat
  Stencil = KernelAbstractions.zeros(backend,Int,0,0)
  Type = Grids.Tri()
  if Type == Grid.Type
    Glob = KernelAbstractions.zeros(backend,Int,1,Grid.NumFaces)
    GlobCPU = zeros(Int,1,Grid.NumFaces)
    for iF = 1 : Grid.NumFaces
      GlobCPU[1,iF] = Grid.Faces[iF].F  
    end  
    copyto!(Glob,GlobCPU)
    NumG = Grid.NumFaces
    NumI = Grid.NumFaces

  end  
 return DG0Struct{FT,
                 typeof(Glob)}(
    Glob,
    Stencil,
    NumG,
    NumI,
    Type,
  )
end

function DG0Struct{FT}(backend,::Grids.Quad) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  Stencil = KernelAbstractions.zeros(backend,Int,0,0)
  Type = Grids.Quad()
  NumG = 0
  NumI = 0
 return DG0Struct{FT,
                 typeof(Glob)}(
    Glob,
    Stencil,
    NumG,
    NumI,
    Type,
  )
end

mutable struct RT0Struct{FT<:AbstractFloat,
                        IT2<:AbstractArray} <: FiniteElement
  DoF::Int                      
  Fun::Array{Polynomial,2}
  Glob::IT2                        
  Stencil::IT2
  NumG::Int
  NumI::Int
  Type::Grids.ElementType
end

function RT0Struct{FT}(backend,::Grids.Tri,Grid) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  Stencil = KernelAbstractions.zeros(backend,Int,0,0)
  Type = Grids.Tri()
  if Type == Grid.Type
    DoF = 3  
    Fun = Array{Polynomial,2}(undef,3,2)
    Fun[1,1] = Polynomial([sqrt(2)/2,sqrt(2)/2])
    Fun[1,2] = Polynomial([sqrt(2)/2,sqrt(2)/2])
    Fun[2,1] = Polynomial([1/2,-1/2])
    Fun[2,2] = Polynomial([1/2,+1/2])
    Fun[3,1] = Polynomial([1/2,+1/2])
    Fun[3,2] = Polynomial([1/2,-1/2])
    Glob = KernelAbstractions.zeros(backend,Int,3,Grid.NumFaces)
    GlobCPU = zeros(Int,3,Grid.NumFaces)
    NumG = Grid.NumEdgesI + Grid.NumEdgesB
    NumI = Grid.NumEdgesI
    for iF = 1 : Grid.NumFaces
      for i = 1 : 3 
         iE = Grid.Faces[iF].E[i] 
         Glob[i,iF] = Grid.Edges[iE].E
      end
    end  
  end    
 return RT0Struct{FT,
                 typeof(Glob)}( 
    DoF,             
    Fun,
    Glob,
    Stencil,
    NumG,
    NumI,
    Type,
  )
end

