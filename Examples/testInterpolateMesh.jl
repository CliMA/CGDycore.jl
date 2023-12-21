import CGDycore:
  Examples, Parallels, Models, Grids, Outputs, Integration,  GPU, DyCore
using MPI
using Base
using CUDA
using AMDGPU
using Metal
using KernelAbstractions
using StaticArrays
using ArgParse
using MPI
using NearestNeighbors
using Distances
using LinearAlgebra
using SparseArrays

EPS = eps(Float64)

struct SphericalAngleCart <: Metric  end 

#spherical_angular_cart(a, b) = 1.0 - cos(atan(norm(cross(a,b))/dot(a,b)))
spherical_angular_cart(a, b) = acos(max(min(dot(a,b)/(norm(a)*norm(b)),1),-1))
(::SphericalAngleCart)(a, b) = spherical_angular_cart(a, b)

FTB = Float64
Problem = "AdvectionSphereSlottedCylinder"
Param = Examples.Parameters(FTB,Problem)
Phys=DyCore.PhysParameters{FTB}()
Profile = Examples.DivergentSphereExample()(Param,Phys)

function main()
  nz = 1
  nPanel = 30
  RefineLevel = 5
  RadEarth = 1.0
  FT = Float64
  backend = CPU()

  SrcGrid = Grids.CubedGrid(backend,FT,nPanel,Grids.OrientFaceSphere,RadEarth,nz)

  IcosahedronGrid = Grids.CreateIcosahedronGrid()
  for iRef = 1 : RefineLevel
    Grids.RefineEdgeTriangularGrid!(IcosahedronGrid)
    Grids.RefineFaceTriangularGrid!(IcosahedronGrid)
  end
  Grids.NumberingTriangularGrid!(IcosahedronGrid)
  DestGrid = Grids.TriangularGridToGrid(backend,FT,IcosahedronGrid,RadEarth,nz)

  Points1 = zeros(3,SrcGrid.NumNodes)
  for i = 1 : SrcGrid.NumNodes
    Points1[1,i] = SrcGrid.Nodes[i].P.x   
    Points1[2,i] = SrcGrid.Nodes[i].P.y   
    Points1[3,i] = SrcGrid.Nodes[i].P.z   
  end  

  tree1 = BallTree(Points1, SphericalAngleCart(), leafsize=10)
  #tree1 = KDTree(Points1, Euclidean(), leafsize=10)
    

  Val = Array{Float64,1}(undef,0)
  RowInd = Array{Int,1}(undef,0)
  ColInd = Array{Int,1}(undef,0)
  for iFD = 1 : DestGrid.NumFaces
    Mid = SVector{3}(DestGrid.Faces[iFD].Mid.x,DestGrid.Faces[iFD].Mid.y,DestGrid.Faces[iFD].Mid.z)
    DestPolygon  = Grids.Polygon(DestGrid.Faces[iFD],DestGrid.Nodes) 
    DestTriangles = Grids.triangulate(DestPolygon)
    areaDest = Grids.area(DestTriangles)
    r = 0.0
    for N in DestGrid.Faces[iFD].N
      P = SVector{3}(DestGrid.Nodes[N].P.x,DestGrid.Nodes[N].P.y,DestGrid.Nodes[N].P.z)
      r = max(r, SphericalAngleCart()(Mid,P))
    end  
    idxsLnn, dist = nn(tree1, Mid)
    r = max(r, 1.5*dist + 10.0 * EPS)
    idxsL = inrange(tree1, Mid, r, true)
    FaceSource = zeros(Int,0)
    for N in idxsL
      for F in SrcGrid.Nodes[N].F
        push!(FaceSource,F)
      end
    end
    FaceSource = unique(FaceSource)
    if iFD == 107
      @show FaceSource  
    end  
    if length(FaceSource) == 0
      @show iFD
      @show idxsL
      @show idxsLnn
      @show r
      stop
    end  
    for iFS in FaceSource
      areaLoc = Grids.intersect(DestGrid.Faces[iFD],DestGrid,SrcGrid.Faces[iFS],SrcGrid,2.0*EPS)
      push!(Val,areaLoc/areaDest)
      push!(RowInd,iFD)
      push!(ColInd,iFS)
    if iFD == 107
      @show areaLoc,areaDest
    end  
    end  
  end
  Inter = sparse(RowInd,ColInd,Val)
  #Test 
  c1 = ones(SrcGrid.NumFaces)
  c2 = Inter * c1
  @show maximum(c2)
  @show minimum(c2)
  for iFD = 1 : DestGrid.NumFaces
    if c2[iFD] < 0.9
      @show iFD,c2[iFD]  
    end
  end  

  for iFS = 1 : SrcGrid.NumFaces
    Mid = SVector{3}(SrcGrid.Faces[iFS].Mid.x,SrcGrid.Faces[iFS].Mid.y,SrcGrid.Faces[iFS].Mid.z)
    _,_,_,_,Tr = Profile(Mid,0.0)
    c1[iFS] = Tr
  end    
  c2 = Inter * c1
  @show size(Inter),DestGrid.NumFaces

  vtkSkeletonMeshSrc = Outputs.vtkStruct{Float64}(backend,SrcGrid)
  Outputs.vtkSkeleton(vtkSkeletonMeshSrc,"Source", 1, 1, c1)

  vtkSkeletonMeshDest = Outputs.vtkStruct{Float64}(backend,DestGrid)
  Outputs.vtkSkeleton(vtkSkeletonMeshDest,"Destination", 1, 1, c2)
end
