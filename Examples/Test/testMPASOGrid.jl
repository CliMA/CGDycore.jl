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

MPI.Init()
comm = MPI.COMM_WORLD
Proc = MPI.Comm_rank(comm) + 1
ProcNumber = MPI.Comm_size(comm)
ParallelCom = DyCore.ParallelComStruct()
ParallelCom.Proc = Proc
ParallelCom.ProcNumber  = ProcNumber

FTB = Float64

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
GridType = "MPASO"
Grid, Exchange = Grids.InitGridSphere(backend,FTB,OrdPoly,nz,nPanel,RefineLevel,GridType,Decomp,RadEarth,Model,ParallelCom;order=false)
vtkSkeletonMesh = Outputs.vtkStruct{Float64}(backend,Grid)
c = ones(FTB,Grid.NumFaces) * Proc
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType, Proc, ProcNumber, c, FileNumber)

KiteGrid = Grids.Grid2KiteGrid(backend,FTB,Grid,Grids.OrientFaceSphere)
vtkSkeletonKite = Outputs.vtkStruct{Float64}(backend,KiteGrid)
pKite = ones(KiteGrid.NumFaces,1)
FileNumber += 1
Outputs.vtkSkeleton!(vtkSkeletonKite, "KiteGrid", Proc, ProcNumber, pKite, FileNumber)

CG1KiteP = FEMSei.CG1KitePrimalStruct{FTB}(Grids.Quad(),backend,KiteGrid)
CG1KiteP.M = FEMSei.MassMatrix(backend,FTB,CG1KiteP,KiteGrid,1,FEMSei.Jacobi!)
FpM = lu(CG1KiteP.M)
CG1KiteDHDiv = FEMSei.CG1KiteDualHDiv{FTB}(Grids.Quad(),backend,KiteGrid)
CG1KiteDHDiv.M = FEMSei.MassMatrix(backend,FTB,CG1KiteDHDiv,KiteGrid,1,FEMSei.Jacobi!)
FuM = lu(CG1KiteDHDiv.M)
Div = FEMSei.DivMatrix(backend,FTB,CG1KiteDHDiv,CG1KiteP,KiteGrid,3,FEMSei.Jacobi!)
Grad = FEMSei.GradMatrix(backend,FTB,CG1KiteP,CG1KiteDHDiv,KiteGrid,3,FEMSei.Jacobi!)

u = zeros(FTB,CG1KiteDHDiv.NumG)
uNeu = zeros(FTB,CG1KiteDHDiv.NumG)

p = FEMSei.Project(backend,FTB,CG1KiteP,KiteGrid,3,FEMSei.Jacobi!,FEMSei.fp)
#u = FEMSei.Project(backend,FTB,CG1KiteD,KiteGrid,3,FEMSei.Jacobi,FEMSei.up)
p0 = deepcopy(p)
FileNumber += 1
pM = FEMSei.ComputeScalar(backend,FTB,CG1KiteP,KiteGrid,p)
VelCa = zeros(KiteGrid.NumFaces,Grid.Dim)
VelSp = zeros(KiteGrid.NumFaces,2)
pM = FEMSei.ComputeScalar(backend,FTB,CG1KiteP,KiteGrid,p)
#FEMSei.ConvertVelocityCart!(backend,FTB,VelCa,u,CG1KiteD,KiteGrid,FEMSei.Jacobi)
Outputs.vtkSkeleton!(vtkSkeletonKite, "KiteGridMPAS", Proc, ProcNumber, [pM VelCa VelSp], FileNumber)
pNeu = similar(p)

nAdveVel = 1000
dtau = 0.001
time = 0.0

rp = similar(p)
ru = similar(u)
for i = 1 : nAdveVel
  @show size(Div),size(u),size(rp)  
  mul!(rp,Div,u)
  ldiv!(FpM,rp)  
  mul!(ru,Div',p)
  @. ru = -ru
  @. ru[CG1KiteDHDiv.ListB] .= 0
  ldiv!(FuM,ru)  

  @. uNeu = u + 0.5 * dtau * ru
  @. pNeu = p + 0.5 * dtau * rp

  mul!(rp,Div,uNeu)
  ldiv!(FpM,rp)  
  mul!(ru,Div',pNeu)
  @. ru = -ru
  @. ru[CG1KiteDHDiv.ListB] .= 0
  ldiv!(FuM,ru)  

  @. u = u + dtau * ru
  @. p = p + dtau * rp

end
FileNumber += 1
pM = FEMSei.ComputeScalar(backend,FTB,CG1KiteP,KiteGrid,p)
FEMSei.ConvertVelocityCart!(backend,FTB,VelCa,u,CG1KiteDHDiv,KiteGrid,FEMSei.Jacobi!)
FEMSei.ConvertVelocitySp!(backend,FTB,VelSp,u,CG1KiteDHDiv,KiteGrid,FEMSei.Jacobi!)
Outputs.vtkSkeleton!(vtkSkeletonKite, "KiteGridMPAS", Proc, ProcNumber, [pM VelCa VelSp], FileNumber)

