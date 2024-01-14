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

backend = CPU()
FT = Float32
nz = 1
Rad = 1.0
RefineLevel = 3

GridTri = Grids.TriangularGrid(backend,FT,RefineLevel,Rad,nz)
vtkSkeletonMeshTri = Outputs.vtkStruct{Float64}(backend,GridTri)
c = ones(FT,GridTri.NumFaces)
Outputs.vtkSkeleton(vtkSkeletonMeshTri,"IcosahedronTri", 1, 1,c)

GridDel = Grids.DelaunayGrid(backend,FT,RefineLevel,Rad,nz)
vtkSkeletonMeshDel = Outputs.vtkStruct{Float64}(backend,GridDel)
c = ones(FT,GridDel.NumFaces)
Outputs.vtkSkeleton(vtkSkeletonMeshDel,"IcosahedronDel", 1, 1,c)


