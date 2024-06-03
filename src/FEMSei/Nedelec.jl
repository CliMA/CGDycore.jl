mutable struct Nedelec0Struct{FT<:AbstractFloat,
                      IT2<:AbstractArray} <: HCurlConfElement
  Glob::IT2
  DoF::Int
  Comp::Int                      
  phi::Array{Polynomial,2}  
  Curlphi::Array{Polynomial,2}                       
  NumG::Int
  NumI::Int
  Type::Grids.ElementType
  M::AbstractSparseMatrix
  LUM::SparseArrays.UMFPACK.UmfpackLU{Float64, Int64}
end

#Nedelec0 Quad

function Nedelec0Struct{FT}(::Grids.Quad,backend,Grid) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  Type = Grids.Quad()
  DoF = 4
  Comp = 2
  @polyvar x1 x2 ksi1 ksi2
  phi = Array{Polynomial,2}(undef,DoF,Comp)
  nu = Array{Polynomial,2}(undef,DoF,Comp)
  Curlphi = Array{Polynomial,2}(undef,DoF,1)

  nu[1,1] = 0.0*ksi1 + 1.0 - 1.0*ksi2
  nu[1,2] = 0.0*ksi1 + 0.0*ksi2

  nu[2,1] = 0.0*ksi1 + 0.0*ksi2
  nu[2,2] = 1.0*ksi1 + 0.0*ksi2

  nu[3,1] = 0.0*ksi1 + 1.0*ksi2
  nu[3,2] = 0.0*ksi1 + 0.0*ksi2

  nu[4,1] = 0.0*ksi1 + 0.0*ksi2
  nu[4,2] = -1.0*ksi1 +1.0 + 0.0*ksi2

  for s = 1 : DoF
    for t = 1 : 2
      phi[s,t] = subs(nu[s,t], ksi1 => (x1+1)/2, ksi2 => (x2+1)/2)
    end
  end

  for i = 1 : DoF
    Curlphi[i,1] = -differentiate(phi[i,1],x2) + differentiate(phi[i,2],x1)
  end

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
  M = sparse([1],[1],[1.0])
  LUM = lu(M)
  return Nedelec0Struct{FT,
                  typeof(Glob)}( 
  Glob,
  DoF,
  Comp,
  phi,                      
  Curlphi,
  NumG,
  NumI,
  Type,
  M,
  LUM,
    )
  end

#Nedelec0 Tri

function Nedelec0Struct{FT}(type::Grids.Tri,backend,Grid) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  Type = Grids.Tri()
  DoF = 3
  Comp = 2
  phi = Array{Polynomial,2}(undef,DoF,Comp) #base function of our used reference element
  nu = Array{Polynomial,2}(undef,DoF,Comp) #base function of standard reference element see defelement.com
  Curlphi = Array{Polynomial,2}(undef,DoF,1)
  @polyvar x1 x2 ksi1 ksi2


  nu[1,1] = 0.0*ksi1 - 1.0*ksi2 + 1.0
  nu[1,2] = 1.0*ksi1 + 0.0*ksi2

  nu[2,1] = 0.0*ksi1 - 1.0*ksi2
  nu[2,2] = 1.0*ksi1 + 0.0*ksi2

  nu[3,1] = 0.0*ksi1 + 1.0*ksi2
  nu[3,2] = -1.0*ksi1 + 0.0*ksi2 + 1.0

  for s = 1 : DoF
    for t = 1 : 2
      phi[s,t] = subs(nu[s,t], ksi1 => (x1+1)/2, ksi2 => (x2+1)/2)
    end
  end

  for i = 1 : DoF
    Curlphi[i,1] = -differentiate(phi[i,1],x2) + differentiate(phi[i,2],x1)
  end


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
  M = sparse([1],[1],[1.0])
  LUM = lu(M)
  return Nedelec0Struct{FT,
                  typeof(Glob)}( 
    Glob,
    DoF,
    Comp,
    phi,
    Curlphi,                      
    NumG,
    NumI,
    Type,
    M,
    LUM,
      )
end

mutable struct Nedelec1Struct{FT<:AbstractFloat,
                      IT2<:AbstractArray} <: HCurlConfElement
  Glob::IT2
  DoF::Int
  Comp::Int                      
  phi::Array{Polynomial,2}  
  Curlphi::Array{Polynomial,2}                       
  NumG::Int
  NumI::Int
  Type::Grids.ElementType
  M::AbstractSparseMatrix
  LUM::SparseArrays.UMFPACK.UmfpackLU{Float64, Int64}
end

#Nedelec1 Quad
#=
__6__7_
|  11  |
3 9 10 5
|  8   |
2      4
|_0__1_|

New
__6__5_
|  12  |
7 10 11 4
|  9   |
8      3
|_1__2_|

=#

function Nedelec1Struct{FT}(::Grids.Quad,backend,Grid) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  Type = Grids.Quad()
  DoF = 12
  Comp = 2
  phi = Array{Polynomial,2}(undef,DoF,Comp) 
  nu = Array{Polynomial,2}(undef,DoF,Comp)
  Curlphi = Array{Polynomial,2}(undef,DoF,1)
  @polyvar x1 x2 ksi1 ksi2

#Nedelec1 DEFelement

  nu[1,1] = -18.0*ksi1*ksi2^2 + 24.0*ksi1*ksi2 - 6.0*ksi1 + 12.0*ksi2^2 - 16.0*ksi2 + 4.0
  nu[1,2] = 0.0*ksi1 + 0.0*ksi2 

  nu[2,1] = 18.0*ksi1*ksi2^2 - 24.0*ksi1*ksi2 + 6.0*ksi1 - 6.0*ksi2^2 + 8.0*ksi2 - 2.0
  nu[2,2] = 0.0*ksi1 + 0.0*ksi2

  nu[3,1] = 0.0*ksi1 + 0.0*ksi2
  nu[3,2] = 2.0*ksi1 * (-9.0*ksi1*ksi2 + 6.0*ksi1 + 6.0*ksi2 - 4.0)

  nu[4,1] = 0.0*ksi1 + 0.0*ksi2
  nu[4,2] = 2.0*ksi1 * (9.0*ksi1*ksi2 - 3.0*ksi1 -6.0*ksi2 + 2.0)

  nu[5,1] = 2.0*ksi2 * (-9.0*ksi1*ksi2 +6.0*ksi1 + 6.0*ksi2 - 4.0)
  nu[5,2] = 0.0*ksi1 + 0.0*ksi2 

  nu[6,1] = 2.0*ksi2 * (9.0*ksi1*ksi2 -6.0*ksi1 - 3.0*ksi2 + 2.0)
  nu[6,2] = 0.0*ksi1 + 0.0*ksi2 

  nu[7,1] = 0.0*ksi1 + 0.0*ksi2
  nu[7,2] = -18.0*ksi1^2*ksi2 + 12.0*ksi1^2 + 24.0*ksi1*ksi2 - 16.0*ksi1 -6.0*ksi2 + 4.0

  nu[8,1] = 0.0*ksi1 + 0.0*ksi2
  nu[8,2] = 18.0*ksi1^2*ksi2 - 6.0*ksi1^2 - 24.0*ksi1*ksi2 + 8.0*ksi1 + 6.0*ksi2 - 2.0

  nu[9,1] = 0.0*ksi1 + 0.0*ksi2
  nu[9,2] = 12.0*ksi1 * (3.0*ksi1*ksi2 - 2.0*ksi1 - 3.0*ksi2 + 2.0)

  nu[10,1] = 12.0*ksi2 * (-3.0*ksi1*ksi2 + 3.0*ksi1 + 2.0*ksi2 - 2.0)
  nu[10,2] = 0.0*ksi1 + 0.0*ksi2

  nu[11,1] = 12.0*ksi2 * (3.0*ksi1*ksi2 - 3.0*ksi1 - 1.0*ksi2 + 1.0) 
  nu[11,2] = 0.0*ksi1 + 0.0*ksi2

  nu[12,1] =  0.0*ksi1 + 0.0*ksi2
  nu[12,2] =  12.0*ksi1 * (-3.0*ksi1*ksi2 + ksi1 + 3.0*ksi2 - 1.0)


  for s = 1 : DoF
    for t = 1 : 2
      phi[s,t] = subs(nu[s,t], ksi1 => (x1+1)/2, ksi2 => (x2+1)/2)
    end
  end
  
  for i = 1 : DoF
    Curlphi[i,1] = -differentiate(phi[i,1],x2) + differentiate(phi[i,2],x1)
  end

  Glob = KernelAbstractions.zeros(backend,Int,DoF,Grid.NumFaces)
  GlobCPU = zeros(Int,DoF,Grid.NumFaces)
  NumG = 2*Grid.NumEdgesI + 2*Grid.NumEdgesB + 4*Grid.NumFaces
  NumI = Grid.NumEdgesI
  for iF = 1 : Grid.NumFaces
    for i = 1 : length(Grid.Faces[iF].E)
      iE = Grid.Faces[iF].E[i]
      GlobCPU[2*i-1,iF] = 2*Grid.Edges[iE].E - 1
      GlobCPU[2*i,iF] = 2*Grid.Edges[iE].E 
    end
    GlobCPU[9,iF] = 2*Grid.NumEdges + 4*Grid.Faces[iF].F - 3
    GlobCPU[10,iF] = 2*Grid.NumEdges + 4*Grid.Faces[iF].F - 2
    GlobCPU[11,iF] = 2*Grid.NumEdges + 4*Grid.Faces[iF].F - 1
    GlobCPU[12,iF] = 2*Grid.NumEdges + 4*Grid.Faces[iF].F
  end
  copyto!(Glob,GlobCPU)
  M = sparse([1],[1],[1.0])
  LUM = lu(M)
  return Nedelec1Struct{FT,
                  typeof(Glob)}( 
    Glob,
    DoF,
    Comp,
    phi,                      
    Curlphi,
    NumG,
    NumI,
    Type,
    M,
    LUM,
      )
  end

#Nedelec1 Tri
#=Numbering

| \
4  2
| 8 \
3 7  1
|_5_6_\

New
| \
5  4
| 8 \
6 7  3
|_1_2_\


=#

function Nedelec1Struct{FT}(type::Grids.Tri,backend,Grid) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  Type = Grids.Tri()
  DoF = 8
  Comp = 2
  phi = Array{Polynomial,2}(undef,DoF,Comp) #base function of our reference triangle
  nu = Array{Polynomial,2}(undef,DoF,Comp) #base function of standard reference element
  Curlphi = Array{Polynomial,2}(undef,DoF,1)
  @polyvar x1 x2 ksi1 ksi2

  nu[1,1] = 8.0*ksi1*ksi2 - 6.0*ksi1 + 8.0*ksi2*ksi2 - 12.0*ksi2 + 4.0
  nu[1,2] = -2.0*ksi1 * (4.0*ksi1 + 4.0*ksi2 - 3.0)

  nu[2,1] = -8.0*ksi1*ksi2 + 6.0*ksi1 + 2.0*ksi2 - 2.0
  nu[2,2] = -4.0*ksi1 * (1.0 - 2.0*ksi1) + 0.0*ksi2

  nu[3,1] = 2.0*ksi2 * (1.0 - 4.0*ksi1)
  nu[3,2] = -4.0*ksi1 * (1.0 - 2.0*ksi1) + 0.0*ksi2 

  nu[4,1] = 4.0*ksi2 * (1.0 - 2.0*ksi2) + 0.0*ksi1
  nu[4,2] = -2.0*ksi1 * (1.0 - 4.0*ksi2) 

  nu[5,1] = 4.0*ksi2 * (2.0*ksi2 - 1.0) + 0.0*ksi1
  nu[5,2] = -8.0*ksi1*ksi2 + 2.0*ksi1 + 6.0*ksi2 - 2.0

  nu[6,1] = 2.0*ksi2 * (-4.0*ksi1 - 4.0*ksi2 + 3.0)
  nu[6,2] = 8.0*ksi1*ksi1 + 8.0*ksi1*ksi2 - 12.0*ksi1 - 6.0*ksi2 + 4.0

  #non-normal

  nu[7,1] = 8.0*ksi2 * (-2.0*ksi1 - 1.0*ksi2 + 1.0)
  nu[7,2] = -8.0*ksi1 * (-2.0*ksi1 - 1.0*ksi2 + 2.0)

  nu[8,1] = 8.0*ksi2 * (-1.0*ksi1 - 2.0*ksi2 + 2.0)
  nu[8,2] = -8.0*ksi1 * (-1.0*ksi1 - 2.0*ksi2 + 1.0)

  for s = 1 : DoF
    for t = 1 : 2
      phi[s,t] = subs(nu[s,t], ksi1 => (x1+1)/2, ksi2 => (x2+1)/2)
    end
  end

  for i = 1 : DoF
    Curlphi[i,1] = -differentiate(phi[i,1],x2) + differentiate(phi[i,2],x1)
  end

  Glob = KernelAbstractions.zeros(backend,Int,DoF,Grid.NumFaces)
  GlobCPU = zeros(Int,DoF,Grid.NumFaces)
  NumG = 2*Grid.NumEdgesI + 2*Grid.NumEdgesB + 2*Grid.NumFaces
  NumI = Grid.NumEdgesI
  for iF = 1 : Grid.NumFaces
      iE1 = Grid.Faces[iF].E[1]
      GlobCPU[1,iF] = 2*Grid.Edges[iE1].E - 1
      GlobCPU[2,iF] = 2*Grid.Edges[iE1].E
      iE2 = Grid.Faces[iF].E[2]
      GlobCPU[3,iF] = 2*Grid.Edges[iE2].E - 1 
      GlobCPU[4,iF] = 2*Grid.Edges[iE2].E 
      iE3 = Grid.Faces[iF].E[3]
      GlobCPU[5,iF] = 2*Grid.Edges[iE3].E 
      GlobCPU[6,iF] = 2*Grid.Edges[iE3].E - 1
    GlobCPU[7,iF] = 2*Grid.NumEdges + 2*Grid.Faces[iF].F - 1
    GlobCPU[8,iF] = 2*Grid.NumEdges + 2*Grid.Faces[iF].F
  end
  copyto!(Glob,GlobCPU)
  M = sparse([1],[1],[1.0])
  LUM = lu(M)
  return Nedelec1Struct{FT,
                  typeof(Glob)}( 
    Glob,
    DoF,
    Comp,
    phi,
    Curlphi,                      
    NumG,
    NumI,
    Type, 
    M,
    LUM,
      )
end

