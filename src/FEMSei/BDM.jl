mutable struct BDM0Struct{FT<:AbstractFloat,
                      IT2<:AbstractArray} <: HDivElement
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

#BDM0 Quad Brezzi-Douglas-Marini

#=Numbering
__6__7_
|      |
3      5
|      |
2      4
|_0__1_|

New
__5__6_
|      |
8      4
|      |
7      3
|_1__2_|

=#

function BDM0Struct{FT}(::Grids.Quad,backend,Grid) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  Type = Grids.Quad()
  DoF = 8
  Comp = 2
  phi = Array{Polynomial,2}(undef,DoF,Comp) 
  nu = Array{Polynomial,2}(undef,DoF,Comp)
  Divphi = Array{Polynomial,2}(undef,DoF,1)
  @polyvar x1 x2 ksi1 ksi2

  nu[1,1] = 3.0*ksi1 * (ksi1 - 1.0) + 0.0*ksi1 + 0.0*ksi2  
  nu[1,2] = 6.0*ksi1*ksi2 - 6.0*ksi1 - 4.0*ksi2 + 4.0

  nu[2,1] = 3.0*ksi1 * (1.0 - ksi1) + 0.0*ksi1 + 0.0*ksi2  
  nu[2,2] = -6.0*ksi1*ksi2 + 6.0*ksi1 + 2.0*ksi2 - 2.0

  nu[3,1] = 2.0*ksi1 * (3.0*ksi2 - 2.0)
  nu[3,2] = 3.0*ksi2 * (1.0*ksi2 - 1.0) + 0.0*ksi1 + 0.0*ksi2 

  nu[4,1] = 2.0*ksi1 * (1.0 - 3.0*ksi2)
  nu[4,2] = 3.0*ksi2 * (1.0 - ksi2) + 0.0*ksi1 + 0.0*ksi2 

  nu[5,1] = 3.0*ksi1 * (1.0 - 1.0*ksi1) + 0.0*ksi1 + 0.0*ksi2 
  nu[5,2] = 2.0*ksi2 * (2.0 - 3.0*ksi1) + 0.0*ksi1 + 0.0*ksi2 

  nu[6,1] = 3.0*ksi1 * (1.0*ksi1 - 1.0) + 0.0*ksi1 + 0.0*ksi2 
  nu[6,2] = 2.0*ksi2 * (3.0*ksi1 - 1.0)

  nu[7,1] = -6.0*ksi1*ksi2 + 4.0*ksi1 + 6.0*ksi2 - 4.0
  nu[7,2] = 3.0*ksi2 * (1.0 - 1.0*ksi2) + 0.0*ksi1 + 0.0*ksi2 

  nu[8,1] = 6.0*ksi1*ksi2 - 2.0*ksi1 - 6.0*ksi2 + 2.0
  nu[8,2] = 3.0*ksi2 * (1.0*ksi2 - 1.0) + 0.0*ksi1 + 0.0*ksi2 

  for s = 1 : DoF
    for t = 1 : 2
      phi[s,t] = subs(nu[s,t], ksi1 => (x1+1)/2, ksi2 => (x2+1)/2)
    end
  end
  
  for i = 1 : DoF
    Divphi[i,1] = differentiate(phi[i,1],x1) + differentiate(phi[i,2],x2)
  end

  Glob = KernelAbstractions.zeros(backend,Int,DoF,Grid.NumFaces)
  GlobCPU = zeros(Int,DoF,Grid.NumFaces)
  NumG = 2*Grid.NumEdgesI + 2*Grid.NumEdgesB
  NumI = Grid.NumEdgesI
  for iF = 1 : Grid.NumFaces
    for i = 1 : length(Grid.Faces[iF].E)
      iE = Grid.Faces[iF].E[i]
      GlobCPU[2*i-1,iF] = 2*Grid.Edges[iE].E - 1
      GlobCPU[2*i,iF] = 2*Grid.Edges[iE].E 
    end
    @show Grid.Edges[iF].F
  end
  copyto!(Glob,GlobCPU)
  M = spzeros(0,0)
  return BDM0Struct{FT,
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
#BDM0 Tri Brezzi-Douglas-Marini
#=Numbering

| \
4  2
|   \
3    1
|_5_6_\

New
| \
5  4
|   \
6    3
|_1_2_\

see RT1
=#

function BDM0Struct{FT}(type::Grids.Tri,backend,Grid) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  Type = Grids.Tri()
  DoF = 6
  Comp = 2
  phi = Array{Polynomial,2}(undef,DoF,Comp) #base function of our reference triangle
  nu = Array{Polynomial,2}(undef,DoF,Comp) #base function of standard reference element
  Divphi = Array{Polynomial,2}(undef,DoF,1)
  @polyvar x1 x2 ksi1 ksi2

  nu[1,1] = 2.0*ksi1 + 0.0*ksi2
  nu[1,2] = -6.0*ksi1 - 4.0*ksi2 + 4.0

  nu[2,1] = -4.0*ksi1 + 0.0*ksi2
  nu[2,2] = 6.0*ksi1 + 2.0*ksi2 - 2.0

  nu[3,1] = -4.0*ksi1 + 0.0*ksi2 
  nu[3,2] = 2.0*ksi2 + 0.0*ksi1

  nu[4,1] = 2.0*ksi1 + 0.0*ksi2 
  nu[4,2] = -4.0*ksi2 + 0.0*ksi1

  nu[5,1] = -2.0*ksi1 - 6.0*ksi2 + 2.0
  nu[5,2] = 4.0*ksi2 + 0.0*ksi1

  nu[6,1] = 4.0*ksi1 + 6.0*ksi2 - 4.0
  nu[6,2] = -2.0*ksi2 + 0.0*ksi1
  for s = 1 : DoF
    for t = 1 : 2
      phi[s,t] = subs(nu[s,t], ksi1 => (x1+1)/2, ksi2 => (x2+1)/2)
    end
  end

  for i = 1 : DoF
    Divphi[i,1] = differentiate(phi[i,1],x1) + differentiate(phi[i,2],x2)
  end

  Glob = KernelAbstractions.zeros(backend,Int,DoF,Grid.NumFaces)
  GlobCPU = zeros(Int,DoF,Grid.NumFaces)
  NumG = 2*Grid.NumEdgesI + 2*Grid.NumEdgesB
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
  end
  copyto!(Glob,GlobCPU)
  M = spzeros(0,0)
  return BDM0Struct{FT,
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

#BDM1 Brezzi-Douglas-Marini

mutable struct BDM1Struct{FT<:AbstractFloat,
                      IT2<:AbstractArray} <: HDivElement
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

#BDM1 Brezzi-Douglas-Marini Quad
#=
 _7_8_9_ 
|       |
12      6
|  13   |
11   14 5
|       |
10      4
|_1_2_3_|

=#

function BDM1Struct{FT}(::Grids.Quad,backend,Grid) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  Type = Grids.Quad()
  DoF = 14
  Comp = 2
  phi = Array{Polynomial,2}(undef,DoF,Comp) 
  nu = Array{Polynomial,2}(undef,DoF,Comp)
  Divphi = Array{Polynomial,2}(undef,DoF,1)
  @polyvar x1 x2 ksi1 ksi2

  nu[13,1] = 6.0*ksi1 * (1.0 - ksi1) + 0.0*ksi2
  nu[13,2] = 0.0 + 0.0*ksi1 + 0.0*ksi2

  nu[14,1] = 0.0 + 0.0*ksi1 + 0.0*ksi2
  nu[14,2] = 6.0*ksi2 * (1.0 - ksi2) + 0.0*ksi1


  nu[1,1] = 5.0*ksi1 * (-2.0*ksi1^2 + 3.0*ksi1 - 1.0) + 0.0*ksi2
  nu[1,2] = -30.0*ksi1^2*ksi2 + 30.0*ksi1^2 + 36.0*ksi1*ksi2 - 36.0*ksi1 + 3.0*ksi2^2 - 12.0*ksi2 + 9

  nu[2,1] = 0.5 * (5.0*ksi1 * (2.0*ksi1^2 - 3.0*ksi1 + 1.0)) + 0.0*ksi2
  nu[2,2] = 15.0*ksi1^2*ksi2 - 15.0*ksi1^2 - 15.0*ksi1*ksi2 + 15.0*ksi1 + 3.0*ksi2^2 - 1.5*ksi2 - 1.5

  nu[3,1] = 5.0*ksi1 * (-2.0*ksi1^2 + 3.0*ksi1 - 1.0) + 0.0*ksi2
  nu[3,2] = -30.0*ksi1^2*ksi2 + 30.0*ksi1^2 + 24.0*ksi1*ksi2 - 24.0*ksi1 + 3.0*ksi2^2 - 6.0*ksi2 + 3.0

  nu[4,1] = 3.0*ksi1 * (-1.0*ksi1 -10.0*ksi2^2 + 12.0*ksi2 - 2.0)
  nu[4,2] = 5.0*ksi2 * (-2.0*ksi2^2 + 3.0*ksi2 - 1.0) + 0.0*ksi1

  nu[5,1] = (3.0*ksi1 * (-2.0*ksi1 + 10.0*ksi2^2 - 10.0*ksi2 + 3.0)) * 0.5 
  nu[5,2] = (5.0*ksi2 * (2.0*ksi2^2 - 3.0*ksi2 + 1.0)) * 0.5  + 0.0*ksi1

  nu[6,1] = 3.0*ksi1 * (-1.0*ksi1 - 10.0*ksi2^2 + 8.0*ksi2)
  nu[6,2] = 5.0*ksi2 * (-2.0*ksi2^2 + 3.0*ksi2 - 1.0) + 0.0*ksi1

  nu[7,1] = 5.0*ksi1 * (2.0*ksi1^2 - 3.0*ksi1 + 1.0) + 0.0*ksi2 
  nu[7,2] = 3.0*ksi2 * (10.0*ksi1^2 - 12.0*ksi1 + 1.0*ksi2 + 2.0)

  nu[8,1] = (5.0*ksi1 * (-2.0*ksi1^2 + 3.0*ksi1 - 1.0)) * 0.5  + 0.0*ksi2
  nu[8,2] = (3.0*ksi2 * (-10.0*ksi1^2 + 10.0*ksi1 + 2.0*ksi2 - 3.0)) * 0.5 + 0.0*ksi1

  nu[9,1] = 5.0*ksi1 * (2.0*ksi1^2 - 3.0*ksi1 + 1.0) + 0.0*ksi2 
  nu[9,2] = 3.0*ksi2 * (10.0*ksi1^2 - 8.0*ksi1 + 1.0*ksi2)

  nu[10,1] = -3.0*ksi1^2 + 30.0*ksi1*ksi2^2 - 36.0*ksi1*ksi2 + 12.0*ksi1 - 30.0*ksi2^2 + 36.0*ksi2 - 9.0
  nu[10,2] = 5.0*ksi2 * (2.0*ksi2^2 - 3.0*ksi2 + 1.0) + 0.0*ksi1

  nu[11,1] = -3.0*ksi1^2 - 15.0*ksi1*ksi2^2 + 15.0*ksi1*ksi2 + 1.5*ksi1 + 15.0*ksi2^2 - 15.0*ksi2 + 1.5
  nu[11,2] = (5.0*ksi2 * (-2.0*ksi2^2 + 3.0*ksi2 - 1.0) + 0.0*ksi1) * 0.5 

  nu[12,1] = -3.0*ksi1^2 + 30.0*ksi1*ksi2^2 - 24.0*ksi1*ksi2 + 6.0*ksi1 - 30.0*ksi2^2 + 24.0*ksi2 - 3.0
  nu[12,2] = 5.0*ksi2 * (2.0*ksi2^2 - 3.0*ksi2 + 1.0) + 0.0*ksi1

  for s = 1 : DoF
    for t = 1 : 2
      phi[s,t] = subs(nu[s,t], ksi1 => (x1+1)/2, ksi2 => (x2+1)/2)
    end
  end
  
  for i = 1 : DoF
    Divphi[i,1] = differentiate(phi[i,1],x1) + differentiate(phi[i,2],x2)
  end

  Glob = KernelAbstractions.zeros(backend,Int,DoF,Grid.NumFaces)
  GlobCPU = zeros(Int,DoF,Grid.NumFaces)
  NumG = 3*Grid.NumEdgesI + 3*Grid.NumEdgesB + 2*Grid.NumFaces
  NumI = Grid.NumEdgesI
  for iF = 1 : Grid.NumFaces
    for i = 1 : length(Grid.Faces[iF].E)
      iE = Grid.Faces[iF].E[i]
      GlobCPU[3*i-1,iF] = 3*Grid.Edges[iE].E - 2
      GlobCPU[3*i-2,iF] = 3*Grid.Edges[iE].E - 1
      GlobCPU[3*i,iF] = 3*Grid.Edges[iE].E 
    end
    GlobCPU[13,iF] = 3*Grid.NumEdges + 2*Grid.Faces[iF].F - 1
    GlobCPU[14,iF] = 3*Grid.NumEdges + 2*Grid.Faces[iF].F
  end
  copyto!(Glob,GlobCPU)
  M = spzeros(0,0)
  return BDM1Struct{FT,
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
#BDM1 Brezzi-Douglas-Marini Tri
#=Numbering

| \
7  6
|   \
8 12 5
| 11  \
9 10   4
|_1_2_3_\

=#

function BDM1Struct{FT}(type::Grids.Tri,backend,Grid) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  Type = Grids.Tri()
  DoF = 12
  Comp = 2
  phi = Array{Polynomial,2}(undef,DoF,Comp) #base function of our reference triangle
  nu = Array{Polynomial,2}(undef,DoF,Comp) #base function of standard reference element
  Divphi = Array{Polynomial,2}(undef,DoF,1)
  @polyvar x1 x2 ksi1 ksi2

  nu[10,1] = 12.0*ksi1 * (-1.0*ksi1 - 4.0*ksi2 + 1.0)
  nu[10,2] = 12.0*ksi2 * (4.0*ksi1 + 1.0*ksi2 - 1.0)

  nu[11,1] = 12.0*ksi1 * (1.0*ksi1 + 2.0*ksi2 - 1.0)
  nu[11,2] = 12.0*ksi2 * (-4.0*ksi1 - 3.0*ksi2 + 3.0)

  nu[12,1] = 12.0*ksi1 * (-3.0*ksi1 - 4.0*ksi2 + 3.0)
  nu[12,2] = 12.0*ksi2 * (2.0*ksi1 + 1.0*ksi2 - 1.0)


  nu[1,1] = 3.0*ksi1 * (-4.0*ksi1 - 4.0*ksi2 + 3.0)
  nu[1,2] = 30.0*ksi1^2 + 48.0*ksi1*ksi2 - 36.0*ksi1 + 18.0*ksi2^2 - 27.0*ksi2 + 9.0

  nu[2,1] = (3.0*ksi1 * (5.0*ksi1 + 6.0*ksi2 - 4.0)) * 0.5 
  nu[2,2] = -15.0*ksi1^2 - 15.0*ksi1*ksi2 + 15.0*ksi1 + 1.5*ksi2^2 - 1.5

  nu[3,1] = 9.0*ksi1 * (1.0 - 2.0*ksi1) + 0.0*ksi2
  nu[3,2] = 30.0*ksi1^2 + 12.0*ksi1*ksi2 - 24.0*ksi1 - 3.0*ksi2 + 3.0
#
  nu[4,1] = 9.0*ksi1 * (1.0 - 2.0*ksi1) + 0.0*ksi2
  nu[4,2] = 3.0*ksi2 * (4.0*ksi1 - 1.0)

  nu[5,1] = (3.0*ksi1 * (-1.0*ksi1 - 6.0*ksi2 + 2.0)) * 0.5
  nu[5,2] = (3.0*ksi2 * (-6.0*ksi1 - 1.0*ksi2 + 2.0)) * 0.5

  nu[6,1] = 3.0*ksi1 * (4.0*ksi2 - 1.0)
  nu[6,2] = 9.0*ksi2 * (1.0 - 2.0*ksi2) + 0.0*ksi1
#
  nu[7,1] = -12.0*ksi1*ksi2 + 3.0*ksi1 - 30.0*ksi2^2 + 24.0*ksi2 - 3.0
  nu[7,2] = 9.0*ksi2 * (2.0*ksi2 - 1.0) + 0.0*ksi1

  nu[8,1] = -1.5*ksi1^2 + 15.0*ksi1*ksi2 + 15.0*ksi2^2 - 15.0*ksi2 + 1.5
  nu[8,2] = (3.0*ksi2 * (-6.0*ksi1 - 5.0*ksi2 + 4.0)) * 0.5

  nu[9,1] = -18.0*ksi1^2 - 48.0*ksi1*ksi2 + 27.0*ksi1 - 30.0*ksi2^2 + 36.0*ksi2 - 9.0
  nu[9,2] = 3.0*ksi2 * (4.0*ksi1 + 4.0*ksi2 - 3.0)

  for s = 1 : DoF
    for t = 1 : 2
      phi[s,t] = subs(nu[s,t], ksi1 => (x1+1)/2, ksi2 => (x2+1)/2)
    end
  end

  for i = 1 : DoF
    Divphi[i,1] = differentiate(phi[i,1],x1) + differentiate(phi[i,2],x2)
  end

  Glob = KernelAbstractions.zeros(backend,Int,DoF,Grid.NumFaces)
  GlobCPU = zeros(Int,DoF,Grid.NumFaces)
  NumG = 3*Grid.NumEdgesI + 3*Grid.NumEdgesB + 3*Grid.NumFaces
  NumI = Grid.NumEdgesI
  for iF = 1 : Grid.NumFaces
      iE1 = Grid.Faces[iF].E[1]
      GlobCPU[1,iF] = 3*Grid.Edges[iE1].E - 1
      GlobCPU[2,iF] = 3*Grid.Edges[iE1].E - 2
      GlobCPU[3,iF] = 3*Grid.Edges[iE1].E
      iE2 = Grid.Faces[iF].E[2]
      GlobCPU[4,iF] = 3*Grid.Edges[iE2].E - 1 
      GlobCPU[5,iF] = 3*Grid.Edges[iE2].E - 2
      GlobCPU[6,iF] = 3*Grid.Edges[iE2].E  
      iE3 = Grid.Faces[iF].E[3]
      GlobCPU[7,iF] = 3*Grid.Edges[iE3].E 
      GlobCPU[8,iF] = 3*Grid.Edges[iE3].E - 2
      GlobCPU[9,iF] = 3*Grid.Edges[iE3].E - 1
    GlobCPU[10,iF] = 3*Grid.NumEdges + 3*Grid.Faces[iF].F - 2
    GlobCPU[11,iF] = 3*Grid.NumEdges + 3*Grid.Faces[iF].F - 1
    GlobCPU[12,iF] = 3*Grid.NumEdges + 3*Grid.Faces[iF].F
  end
  copyto!(Glob,GlobCPU)
  M = spzeros(0,0)
  return BDM1Struct{FT,
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