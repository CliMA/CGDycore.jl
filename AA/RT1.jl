mutable struct RT1Struct{FT<:AbstractFloat,
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

function RT1Struct{FT}(::Grids.Quad,backend,Grid) where FT<:AbstractFloat
Glob = KernelAbstractions.zeros(backend,Int,0,0)
Type = Grids.Quad()
DoF = 12
Comp = 2
@polyvar x1 x2
phi = Array{Polynomial,2}(undef,DoF,Comp)
Divphi = Array{Polynomial,2}(undef,DoF,1)
phi[1,1] = 0.0*x1 + 0.0*x2
phi[1,2] = -0.5 + 0.5*x2  + 0.0*x1

phi[2,1] = 0.5 + 0.5*x1 + 0.0*x2
phi[2,2] = 0.0*x1 + 0.0*x2

phi[3,1] = 0.0*x1 + 0.0*x2
phi[3,2] = -0.5 - 0.5*x2 + 0.0*x1

phi[4,1] = 0.5 - 0.5*x1 + 0.0*x2
phi[4,2] = 0.0*x1 + 0.0*x2
#ToDO 
phi[5,1] = 0.0*x1 + 0.0*x2
phi[5,2] = -0.5 + 0.5*x2  + 0.0*x1

phi[6,1] = 0.5 + 0.5*x1 + 0.0*x2
phi[6,2] = 0*x1 + 0.0*x2

phi[7,1] = 0*x1 + 0.0*x2
phi[7,2] = -0.5 - 0.5*x2 + 0.0*x1

phi[8,1] = 0.5 - 0.5*x1 + 0.0*x2
phi[8,2] = 0*x1 + 0.0*x2

phi[9,1] = 0.0*x1 + 0.0*x2
phi[9,2] = -0.5 + 0.5*x2  + 0.0*x1

phi[10,1] = 0.5 + 0.5*x1 + 0.0*x2
phi[10,2] = 0*x1 + 0.0*x2

phi[11,1] = 0*x1 + 0.0*x2
phi[11,2] = -0.5 - 0.5*x2 + 0.0*x1

phi[12,1] = 0.5 - 0.5*x1 + 0.0*x2
phi[12,2] = 0*x1 + 0.0*x2

for i = 1 : DoF
Divphi[i,1] = differentiate(phi[i,1],x1) + differentiate(phi[i,2],x2)
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
M = spzeros(0,0)
return RT1Struct{FT,
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

function RT1Struct{FT}(type::Grids.Tri,backend,Grid) where FT<:AbstractFloat
Glob = KernelAbstractions.zeros(backend,Int,0,0)
Type = Grids.Tri()
DoF = 8
Comp = 2
phi = Array{Polynomial,2}(undef,DoF,Comp)
Divphi = Array{Polynomial,2}(undef,DoF,1)
@polyvar x1 x2 l1x1 l2x2
g1 = -(1/sqrt(3))
g2 = (1/sqrt(3))
l1x1 = (x1-g2)/(g1-g2)
l2x1 = (x1-g1)/(g2-g1)
l1x2 = (x2-g2)/(g1-g2)
l2x2 = (x2-g1)/(g2-g1)

phi[1,1] = l1x1 * (0.5 + 0.5*x1 + 0.0*x2)
phi[1,2] = -0.5 + 0.5*x2 + 0.0*x1

phi[2,1] = l2x1 * (0.5 + 0.5*x1 + 0.0*x2)
phi[2,2] = -0.5 + 0.5*x2 + 0.0*x1

phi[3,1] = l2x1 * (0.5 + 0.5*x1 + 0.0*x2)
phi[3,2] = 0.5 + 0.5*x2 + 0.0*x1

phi[4,1] = l1x1 * (0.5 + 0.5*x1 + 0.0*x2)
phi[4,2] = 0.5 + 0.5*x2 + 0.0*x1

phi[5,1] = 0.5 - 0.5*x1 + 0.0*x2
phi[5,2] = l1x2 * (-0.5 - 0.5*x2 + 0.0*x1)

phi[6,1] = 0.5 - 0.5*x1 + 0.0*x2
phi[6,2] = l2x2 * (-0.5 - 0.5*x2 + 0.0*x1)

Divphi[1,1] = differentiate(phi[1,1],x1) + differentiate(phi[1,2],x2)
Divphi[2,1] = differentiate(phi[2,1],x1) + differentiate(phi[2,2],x2)
Divphi[3,1] = differentiate(phi[3,1],x1) + differentiate(phi[3,2],x2)
Divphi[4,1] = differentiate(phi[4,1],x1) + differentiate(phi[4,2],x2)
Divphi[5,1] = differentiate(phi[5,1],x1) + differentiate(phi[5,2],x2)
Divphi[6,1] = differentiate(phi[6,1],x1) + differentiate(phi[6,2],x2)

#non-normal

phi[7,1] = x1 * (0.5 + 0.5*x1)
phi[7,2] = x1 * (-0.5 + 0.5*x2)

phi[8,1] = x2 * (0.5 + 0.5*x1)
phi[8,2] = x2 * (0.5 + 0.5*x2)

Divphi[7,1] = differentiate(phi[7,1],x1) + differentiate(phi[7,2],x2)
Divphi[8,1] = differentiate(phi[8,1],x1) + differentiate(phi[8,2],x2)


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
M = spzeros(0,0)
return RT1Struct{FT,
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
