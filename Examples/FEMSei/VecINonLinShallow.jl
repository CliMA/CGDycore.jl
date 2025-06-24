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
using LinearAlgebra

# Model
parsed_args = DyCore.parse_commandline()
Problem = parsed_args["Problem"]
ProfRho = parsed_args["ProfRho"]
ProfTheta = parsed_args["ProfTheta"]
PertTh = parsed_args["PertTh"]
ProfVel = parsed_args["ProfVel"]
ProfVelGeo = parsed_args["ProfVelGeo"]
RhoVPos = parsed_args["RhoVPos"]
RhoCPos = parsed_args["RhoCPos"]
RhoIPos = parsed_args["RhoIPos"]
RhoRPos = parsed_args["RhoRPos"]
HorLimit = parsed_args["HorLimit"]
Upwind = parsed_args["Upwind"]
Damping = parsed_args["Damping"]
Relax = parsed_args["Relax"]
StrideDamp = parsed_args["StrideDamp"]
Geos = parsed_args["Geos"]
Coriolis = parsed_args["Coriolis"]
CoriolisType = parsed_args["CoriolisType"]
Buoyancy = parsed_args["Buoyancy"]
Equation = parsed_args["Equation"]
RefProfile = parsed_args["RefProfile"]
ProfpBGrd = parsed_args["ProfpBGrd"]
ProfRhoBGrd = parsed_args["ProfRhoBGrd"]
Microphysics = parsed_args["Microphysics"]
TypeMicrophysics = parsed_args["TypeMicrophysics"]
RelCloud = parsed_args["RelCloud"]
Rain = parsed_args["Rain"]
Source = parsed_args["Source"]
Forcing = parsed_args["Forcing"]
VerticalDiffusion = parsed_args["VerticalDiffusion"]
JacVerticalDiffusion = parsed_args["JacVerticalDiffusion"]
JacVerticalAdvection = parsed_args["JacVerticalAdvection"]
SurfaceFlux = parsed_args["SurfaceFlux"]
SurfaceFluxMom = parsed_args["SurfaceFluxMom"]
NumV = parsed_args["NumV"]
NumTr = parsed_args["NumTr"]
Curl = parsed_args["Curl"]
ModelType = parsed_args["ModelType"]
Thermo = parsed_args["Thermo"]

# Parallel
Decomp = parsed_args["Decomp"]

# Time integration
SimDays = parsed_args["SimDays"]
SimHours = parsed_args["SimHours"]
SimMinutes = parsed_args["SimMinutes"]
SimSeconds = parsed_args["SimSeconds"]
SimTime = parsed_args["SimTime"]
StartAverageDays = parsed_args["StartAverageDays"]
dtau = parsed_args["dtau"]
IntMethod = parsed_args["IntMethod"]
Table = parsed_args["Table"]

# Grid
nz = parsed_args["nz"]
nPanel = parsed_args["nPanel"]
ns = parsed_args["ns"]
RefineLevel = parsed_args["RefineLevel"]
nLon = parsed_args["nLon"]
nLat = parsed_args["nLat"]
LatB = parsed_args["LatB"]
H = parsed_args["H"]
Stretch = parsed_args["Stretch"]
StretchType = parsed_args["StretchType"]
TopoS = parsed_args["TopoS"]
GridForm = parsed_args["GridForm"]
GridType = parsed_args["GridType"]
RadEarth = parsed_args["RadEarth"]

# flat
nx = parsed_args["nx"]
ny = parsed_args["ny"]
nz = parsed_args["nz"]
H = parsed_args["H"]
Lx = parsed_args["Lx"]
Ly = parsed_args["Ly"]
x0 = parsed_args["x0"]
y0 = parsed_args["y0"]
BoundaryWE = parsed_args["BoundaryWE"]
BoundarySN = parsed_args["BoundarySN"]
BoundaryBT = parsed_args["BoundaryBT"]

# CG Element
OrdPoly = parsed_args["OrdPoly"]

# Viscosity
HyperVisc = parsed_args["HyperVisc"]
HyperDCurl = parsed_args["HyperDCurl"]
HyperDGrad = parsed_args["HyperDGrad"]
HyperDRhoDiv = parsed_args["HyperDRhoDiv"]
HyperDDiv = parsed_args["HyperDDiv"]
HyperDDivW = parsed_args["HyperDDivW"]

# Output
PrintDays = parsed_args["PrintDays"]
PrintHours = parsed_args["PrintHours"]
PrintMinutes = parsed_args["PrintMinutes"]
PrintSeconds = parsed_args["PrintSeconds"]
PrintTime = parsed_args["PrintTime"]
PrintStartTime = parsed_args["PrintStartTime"]
Flat = parsed_args["Flat"]
vtkFileName = parsed_args["vtkFileName"]

# Device
Device = parsed_args["Device"]
GPUType = parsed_args["GPUType"]
FloatTypeBackend = parsed_args["FloatTypeBackend"]
NumberThreadGPU = parsed_args["NumberThreadGPU"]

# Finite elements
k = parsed_args["OrderFEM"]

# Grid Output Refine
ref = parsed_args["RefineOutput"]

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

# Grid construction
RadEarth = Phys.RadEarth

# Parameters
Param = Examples.Parameters(FTB,Problem)

if GridForm == "Spherical"
  # Grid construction
  Grid, CellToProc = Grids.InitGridSphere(backend,FTB,OrdPoly,nz,nPanel,RefineLevel,ns,
  nLat,nLon,LatB,GridType,Decomp,RadEarth,Model,ParallelCom;ChangeOrient=3)

  # Jacobi
  Jacobi = FEMSei.Jacobi!

  # Settings
  GridLengthMin,GridLengthMax = Grids.GridLength(Grid)
  cS = Param.cS
  dtau = GridLengthMin / cS / sqrt(2) * .2 / (k + 1)
  EndTime = SimTime + 3600*24*SimDays + 3600 * SimHours + 60 * SimMinutes + SimSeconds
  nAdveVel::Int = round(EndTime / dtau)
  dtau = EndTime / nAdveVel
  PrintT = PrintTime + 3600*24*PrintDays + 3600 * PrintHours + 60 * PrintMinutes + PrintSeconds
  nPrint::Int = ceil(PrintT/dtau)
  FileNameOutput = vtkFileName
  @show GridLengthMin,GridLengthMax
  @show nAdveVel
  @show dtau
  @show nPrint
  
else
  # Grid construction
  Boundary = Grids.Boundary()
  Boundary.WE = BoundaryWE
  Boundary.SN = BoundarySN
  Boundary.BT = BoundaryBT
  Grid, Exchange = Grids.InitGridCart(backend,FTB,OrdPoly,nx,ny,Lx,Ly,x0,y0,Boundary,nz,Model,ParallelCom;GridType=GridType)
  
  # Jacobi
  Jacobi = FEMSei.JacobiCart!

  # Settings
  GridLengthMin,GridLengthMax = Grids.GridLength(Grid)
  cS = Param.cS
  dtau = 0.05 * 2π / sqrt(nx * ny) / (k + 1)
  EndTime = SimTime + 3600*24*SimDays + 3600 * SimHours + 60 * SimMinutes + SimSeconds
  nAdveVel = round(EndTime / dtau)
  dtau = EndTime / nAdveVel
  PrintT = PrintTime + 3600*24*PrintDays + 3600 * PrintHours + 60 * PrintMinutes + PrintSeconds
  nPrint = ceil(PrintT/dtau)
  nAdveVel = 4000
  nPrint = 400
  FileNameOutput = vtkFileName
  @show GridLengthMin,GridLengthMax
  @show nAdveVel
  @show dtau
  @show nPrint
end

Examples.InitialProfile!(backend,FTB,Model,Problem,Param,Phys)

# Output
vtkSkeletonMesh = Outputs.vtkStruct{Float64}(backend,Grid,Grid.NumFaces,Flat;Refine=ref)

# Quadrature rules
if Grid.Type == Grids.Quad()
  nQuad = 2
  nQuadM = 2
  nQuadS = 2
elseif Grid.Type == Grids.Tri()
  nQuad = 4
  nQuadM = 4
  nQuadS = 4
end

# Finite elements
DG = FEMSei.DGStruct{FTB}(backend,k,Grid.Type,Grid)
CG = FEMSei.CGStruct{FTB}(backend,k+1,Grid.Type,Grid)
RT = FEMSei.RTStruct{FTB}(backend,k,Grid.Type,Grid)
ND = FEMSei.NDStruct{FTB}(backend,k,Grid.Type,Grid)

ModelFEM = FEMSei.ModelFEMVecI(backend,FTB,ND,RT,CG,DG,Grid,nQuadM,nQuadS,Jacobi)

pPosS = ModelFEM.pPosS
pPosE = ModelFEM.pPosE
uPosS = ModelFEM.uPosS
uPosE = ModelFEM.uPosE
U = zeros(FTB,ModelFEM.DG.NumG+ModelFEM.RT.NumG)
@views Up = U[pPosS:pPosE]
@views Uu = U[uPosS:uPosE]

# Interpolation
FEMSei.InterpolateDG!(Up,DG,Jacobi,Grid,Grid.Type,Model.InitialProfile)
FEMSei.InterpolateRT!(Uu,RT,Jacobi,Grid,Grid.Type,nQuad,Model.InitialProfile)

# Time integration
FEMSei.TimeStepperVecI(backend,FTB,U,dtau,FEMSei.FcnVecINonLinShallow!,ModelFEM,Grid,nQuadM,nQuadS,Jacobi,
  nAdveVel,FileNameOutput,Proc,ProcNumber,nPrint,Flat,ref)

