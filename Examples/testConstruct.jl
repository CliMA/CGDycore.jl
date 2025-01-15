import CGDycore:
  Examples, Parallels, Models, Grids, Outputs, Integration,  GPU, DyCore, FEMSei, FiniteVolumes
using MPI
using Base
using CUDA
using AMDGPU
using Metal
using KernelAbstractions
using StaticArrays
using ArgParse
using DynamicPolynomials


backend = CPU()
FTB = Float64
MPI.Init()
comm = MPI.COMM_WORLD
Proc = MPI.Comm_rank(comm) + 1
ProcNumber = MPI.Comm_size(comm)
ParallelCom = DyCore.ParallelComStruct()
ParallelCom.Proc = Proc
ParallelCom.ProcNumber  = ProcNumber
#ModelParameters
Phys = DyCore.PhysParameters{FTB}()
Model = DyCore.ModelStruct{FTB}()

#Grid construction
RefineLevel = 3
nz = 1
nQuad = 3
nQuadM = 3 #2
nQuadS = 3 #3
Decomp = "EqualArea"
nLat = 0
nLon = 0
LatB = 0.0

#Quad
GridType = "TriangularSphere"
nPanel =  1
OrdPoly = 1
RadEarth = 1.0
#GridType = "HealPix"
ns = 10
Grid, Exchange = Grids.InitGridSphere(backend,FTB,OrdPoly,nz,nPanel,RefineLevel,ns,nLat,nLon,LatB,GridType,Decomp,RadEarth,
  Model,ParallelCom)
RT = FEMSei.RT1Struct{FTB}(Grids.Tri(),backend,Grid)
@show RT.phi
print(RT.phi)
s = @polyvar x[1:2]
phi,PRT = FEMSei.ConstructRT_k(1)
u = zeros(RT.NumG)
uN = zeros(RT.NumG)
uN1 = zeros(RT.NumG)

Problem = "AdvectionSphereSpherical"
Param = Examples.Parameters(FTB,Problem)
Examples.InitialProfile!(Model,Problem,Param,Phys)
QuadOrd = 3
#FEMSei.Interpolate!(backend,FTB,uN,RT,Grid.Type,Grid,QuadOrd,FEMSei.Jacobi!,Model.InitialProfile)

#FEMSei.InterpolateCons!(backend,FTB,uN1,RT,Grid,QuadOrd,FEMSei.Jacobi!,Model.InitialProfile)

FEMSei.InterpolateRT!(u,RT,FEMSei.Jacobi!,Grid,Model.InitialProfile)
@show minimum(u),maximum(u)

