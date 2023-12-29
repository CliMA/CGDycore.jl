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


function main()
  FTB = Float64
  Problem = "AdvectionSphereSlottedCylinder"
  Param = Examples.Parameters(FTB,Problem)
  Phys=DyCore.PhysParameters{FTB}()
  Profile = Examples.DivergentSphereExample()(Param,Phys)

  nz = 1
  nPanel = 30
  RefineLevel = 5
  RadEarth = 1.0
  FT = Float64
  backend = CPU()

  SrcGrid = Grids.CubedGrid(backend,FT,nPanel,Grids.OrientFaceSphere,RadEarth,nz)
# DestGrid = Grids.TriangularGrid(backend,FT,RefineLevel,RadEarth,nz)
  DestGrid = Grids.DelaunayGrid(backend,FT,RefineLevel,RadEarth,nz)

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
  Outputs.vtkSkeleton(vtkSkeletonMeshSrc,"Source1CuSp", 1, 1, c1)
  vtkSkeletonMeshDest = Outputs.vtkStruct{Float64}(backend,DestGrid)
  Outputs.vtkSkeleton(vtkSkeletonMeshDest,"Destination1Del", 1, 1, c2)

  nPanelSrc = 40
  nPanelDest = 47
  SrcGrid = Grids.CubedGrid(backend,FT,nPanelSrc,Grids.OrientFaceSphere,RadEarth,nz)
  DestGrid = Grids.CubedGrid(backend,FT,nPanelDest,Grids.OrientFaceSphere,RadEarth,nz)

  Inter = Grids.interpolate(SrcGrid,DestGrid)
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
  vtkSkeletonMeshSrc = Outputs.vtkStruct{Float64}(backend,SrcGrid)
  Outputs.vtkSkeleton(vtkSkeletonMeshSrc,"Source2CuSp", 1, 1, c1)
  vtkSkeletonMeshDest = Outputs.vtkStruct{Float64}(backend,DestGrid)
  Outputs.vtkSkeleton(vtkSkeletonMeshDest,"Destination2CuSp", 1, 1, c2)
  return nothing
end
