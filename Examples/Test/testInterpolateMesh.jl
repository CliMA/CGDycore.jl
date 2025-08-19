import CGDycore:
  Examples, Parallels, Models, Grids, Outputs, Integration, GPU, DyCore
using MPI
using Base
using CUDA
using AMDGPU
using Metal
using KernelAbstractions
using StaticArrays
using ArgParse
using MPI


function Interpolate(SrcGrid,DestGrid,Profile,Example)

  Inter = Grids.interpolate(SrcGrid,DestGrid)
  #Test 
  c1 = ones(SrcGrid.NumFaces)
  c2 = Inter * c1
  @show maximum(c2)
  @show minimum(c2)
  for iFD = 1 : DestGrid.NumFaces
    if c2[iFD] < 0.9999999
      @show iFD,c2[iFD]  
    end
  end  
  for iFS = 1 : SrcGrid.NumFaces
    Mid = SVector{3}(SrcGrid.Faces[iFS].Mid.x,SrcGrid.Faces[iFS].Mid.y,SrcGrid.Faces[iFS].Mid.z)
    _,_,_,_,Tr = Profile(Mid,0.0)
    c1[iFS] = Tr
  end
  c2 = Inter * c1
  vtkSkeletonMeshSrc = Outputs.vtkStruct{Float64}(backend,SrcGrid)
  Outputs.vtkSkeleton(vtkSkeletonMeshSrc,Example*"Source", 1, 1, c1)
  vtkSkeletonMeshDest = Outputs.vtkStruct{Float64}(backend,DestGrid)
  Outputs.vtkSkeleton(vtkSkeletonMeshDest,Example*"Destin", 1, 1, c2)
end  

FTB = Float64
backend = CPU()

Problem = "AdvectionSphereSlottedCylinder"
Param = Examples.Parameters(FTB,Problem)
Phys=DyCore.PhysParameters{FTB}()
Profile = Examples.DivergentSphereExample()(Param,Phys)

nz = 1
nPanel = 30
RefineLevel = 5
RadEarth = 1.0

SrcGrid = Grids.CubedGrid(backend,FTB,nPanel,Grids.OrientFaceSphere,RadEarth,nz)
DestGrid = Grids.DelaunayGrid(backend,FTB,RefineLevel,RadEarth,nz)
Interpolate(SrcGrid,DestGrid,Profile,"CubeDelau")

nPanelSrc = 40
nPanelDest = 47
SrcGrid = Grids.CubedGrid(backend,FTB,nPanelSrc,Grids.OrientFaceSphere,RadEarth,nz)
DestGrid = Grids.CubedGrid(backend,FTB,nPanelDest,Grids.OrientFaceSphere,RadEarth,nz)
Interpolate(SrcGrid,DestGrid,Profile,"CubeCube")

nPanelSrc = 40
nLon = 200
nLat = 100
LatB = (1.0 - 0.001) * pi / 2
SrcGrid = Grids.CubedGrid(backend,FTB,nPanelSrc,Grids.OrientFaceSphere,RadEarth,nz)
DestGrid = Grids.SphericalGrid(backend,FTB,nLon,nLat,LatB,Grids.OrientFaceSphere,RadEarth,nz)
Interpolate(SrcGrid,DestGrid,Profile,"CubeSphere")
