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

MPI.Init()
comm = MPI.COMM_WORLD
Proc = MPI.Comm_rank(comm) + 1
ProcNumber = MPI.Comm_size(comm)
ParallelCom = DyCore.ParallelComStruct()
ParallelCom.Proc = Proc
ParallelCom.ProcNumber  = ProcNumber

FTB = Float32

# Physical parameters
Phys = DyCore.PhysParameters{FTB}()

#ModelParameters
Model = DyCore.ModelStruct{FTB}()

backend = CPU()
nz = 1
Rad = 1.0
RefineLevel = 6
OrdPoly = 3
nPanel = 10
RadEarth = 1.0
Decomp = "EqualArea"

FileNumber = 1
GridType = "MPAS"
Grid, Exchange = Grids.InitGridSphere(backend,FTB,OrdPoly,nz,nPanel,RefineLevel,GridType,Decomp,RadEarth,Model,ParallelCom)
vtkSkeletonMesh = Outputs.vtkStruct{Float64}(backend,Grid)
c = ones(FTB,Grid.NumFaces) * Proc
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType, Proc, ProcNumber , c, FileNumber)

@show "MPASO"
GridType = "MPASO"
Grid, Exchange = Grids.InitGridSphere(backend,FTB,OrdPoly,nz,nPanel,RefineLevel,GridType,Decomp,RadEarth,Model,ParallelCom)
vtkSkeletonMesh = Outputs.vtkStruct{Float64}(backend,Grid)
c = ones(FTB,Grid.NumFaces) * Proc
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType, Proc, ProcNumber, c, FileNumber)

GridType = "Msh"
Grid, Exchange = Grids.InitGridSphere(backend,FTB,OrdPoly,nz,nPanel,RefineLevel,GridType,Decomp,RadEarth,Model,ParallelCom)
vtkSkeletonMesh = Outputs.vtkStruct{Float64}(backend,Grid)
c = ones(FTB,Grid.NumFaces) * Proc
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType, Proc, ProcNumber , c, FileNumber)

GridType = "DelaunaySphere"
Grid, Exchange = Grids.InitGridSphere(backend,FTB,OrdPoly,nz,nPanel,RefineLevel,GridType,Decomp,RadEarth,Model,ParallelCom)
vtkSkeletonMesh = Outputs.vtkStruct{Float64}(backend,Grid)
c = ones(FTB,Grid.NumFaces) * Proc
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType, Proc, ProcNumber , c, FileNumber)

GridType = "TriangularSphere"
Grid, Exchange = Grids.InitGridSphere(backend,FTB,OrdPoly,nz,nPanel,RefineLevel,GridType,Decomp,RadEarth,Model,ParallelCom)
vtkSkeletonMesh = Outputs.vtkStruct{Float64}(backend,Grid)
c = ones(FTB,Grid.NumFaces) * Proc
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType, Proc, ProcNumber , c, FileNumber)

GridType = "SQuadGen"
Grid, Exchange = Grids.InitGridSphere(backend,FTB,OrdPoly,nz,nPanel,RefineLevel,GridType,Decomp,RadEarth,Model,ParallelCom)
vtkSkeletonMesh = Outputs.vtkStruct{Float64}(backend,Grid)
c = ones(FTB,Grid.NumFaces) * Proc
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType, Proc, ProcNumber , c, FileNumber)

GridType = "HealPix"
Grid, Exchange = Grids.InitGridSphere(backend,FTB,OrdPoly,nz,nPanel,RefineLevel,GridType,Decomp,RadEarth,Model,ParallelCom)
vtkSkeletonMesh = Outputs.vtkStruct{Float64}(backend,Grid)
c = ones(FTB,Grid.NumFaces) * Proc
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType, Proc, ProcNumber , c, FileNumber)



