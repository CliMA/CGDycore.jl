import CGDycore:
  Parameters, Examples, Parallels, Models, Grids, FiniteElements, Outputs, Integration,  GPU, DyCore, FEM, FiniteVolumes
using MPI
using Base
using CUDA
using AMDGPU
using Metal
using KernelAbstractions
using StaticArrays
using ArgParse
using LinearAlgebra

GridForm = "Spherical"
VelocityForm = "Spherical"
Problem = "GalewskySphere"
GridType = "CubedSphere" 
OrdPoly = 4
nz = 2
nPanel = 10
RefineLevel = 5
ns = 20
nLat = 0
nLon = 0
LatB = 0
Decomp = "EqualArea"
k = 1
Flat = true
H = 1.0
OrdPoly = 5
OrdPolyZ = 5
OrdPrint = 1
OrdPrintZ = 1
NumberThreadGPU = 256

MPI.Init()
Device = "CPU"
FloatTypeBackend = "Float64"

if Device == "CPU" 
  backend = CPU()
elseif Device == "GPU" 
  if GPUType == "CUDA"
    backend = CUDABackend()
    CUDA.allowscalar(true)
#   CUDA.device!(MPI.Comm_rank(MPI.COMM_WORLD))
  elseif GPUType == "AMD"
    backend = ROCBackend()
    AMDGPU.allowscalar(false)
  elseif GPUType == "Metal"
    backend = MetalBackend()
    Metal.allowscalar(true)
  end
else
  backend = CPU()
end

if FloatTypeBackend == "Float64"
  FTB = Float64
elseif FloatTypeBackend == "Float32"
  FTB = Float32
else
  @show "False FloatTypeBackend"
  stop
end

MPI.Init()
comm = MPI.COMM_WORLD
Proc = MPI.Comm_rank(comm) + 1
ProcNumber = MPI.Comm_size(comm)
ParallelCom = DyCore.ParallelComStruct()
ParallelCom.Proc = Proc
ParallelCom.ProcNumber  = ProcNumber

# Physical parameters
Phys = DyCore.PhysParameters{FTB}()

# ModelParameters
Model = DyCore.ModelStruct{FTB}()
Model.MetricType = "DGMetricCurl"

# Grid construction
RadEarth = Phys.RadEarth

# Parameters
Param = Examples.Parameters(FTB,Problem)

if VelocityForm == "Spherical"
   VelForm = Examples.VelocityS()
elseif VelocityForm == "Cartesian"
   VelForm = Examples.VelocityC()
end

# Grid construction
  Grid, CellToProc = Grids.InitGridSphere(backend,FTB,OrdPoly,nz,nPanel,RefineLevel,ns,
  nLat,nLon,LatB,GridType,Decomp,RadEarth,Model,ParallelCom;ChangeOrient=2)
Grid.AdaptGrid = Grids.AdaptGrid(FTB,"Sleve",FTB(H))
Grids.AddVerticalGrid!(Grid,nz,H)



DG = FiniteElements.DGQuad{FTB}(backend,OrdPoly,OrdPolyZ,OrdPrint,OrdPrint,Grid,ParallelCom.Proc)
wH = zeros(DG.DoF)
iDoF = 0
for i = 1 : OrdPoly + 1
  for j = 1 : OrdPoly +1
    global iDoF += 1  
    wH[iDoF] = DG.w[i] * DG.w[j]  
  end
end  

Exchange = Parallels.ExchangeStruct{FTB}(backend,Grid,DG,CellToProc,Proc,ProcNumber,
      false;Discretization="DG")


zS = zeros(DG.DoF,Grid.NumFaces)
MetricFV3 = FiniteVolumes.MetricVolume3(backend,FTB,Grid)
Grids.EdgesInNodes!(Grid.Nodes,Grid.Edges,Grid.Faces)
FiniteVolumes.Grad3D!(MetricFV3,Grid)
stop


