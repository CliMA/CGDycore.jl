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
GridType = parsed_args["GridType"]
RadEarth = parsed_args["RadEarth"]
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

# Grid Output Refine
ref = parsed_args["RefineOutput"]

# Finite elements
k = parsed_args["OrderFEM"]
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

#ModelParameters
Model = DyCore.ModelStruct{FTB}()

#Grid construction
RadEarth = Phys.RadEarth
Grid, Exchange = Grids.InitGridSphere(backend,FTB,OrdPoly,nz,nPanel,RefineLevel,ns,
  nLat,nLon,LatB,GridType,Decomp,RadEarth,Model,ParallelCom;ChangeOrient=3)


Param = Examples.Parameters(FTB,Problem)

if Problem == "GalewskySphere"
  GridLengthMin,GridLengthMax = Grids.GridLength(Grid)
  cS = sqrt(Phys.Grav * Param.H0G)
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
elseif Problem == "HaurwitzSphere"
  GridLengthMin,GridLengthMax = Grids.GridLength(Grid)
  cS = sqrt(Phys.Grav * Param.h0)
  dtau = GridLengthMin / cS / sqrt(2) * .2 / (k + 1)
  EndTime = SimTime + 3600*24*SimDays + 3600 * SimHours + 60 * SimMinutes + SimSeconds
  nAdveVel::Int = round(EndTime / dtau)
  dtau = EndTime / nAdveVel
  PrintT = PrintTime + 3600*24*PrintDays + 3600 * PrintHours + 60 * PrintMinutes + PrintSeconds
  nPrint::Int =ceil(PrintT/dtau)
  FileNameOutput = vtkFileName
  @show GridLengthMin,GridLengthMax
  @show nAdveVel
  @show dtau
  @show nPrint
else
    print("Error")
end

Examples.InitialProfile!(backend,FTB,Model,Problem,Param,Phys)

#Output
vtkSkeletonMesh = Outputs.vtkStruct{Float64}(backend,Grid,Grid.NumFaces,Flat;Refine=ref)

#Quadrature rules
if Grid.Type == Grids.Quad()
  nQuad = 2
  nQuadM = 2
  nQuadS = 2
elseif Grid.Type == Grids.Tri()
  nQuad = 3
  nQuadM = 3
  nQuadS = 3
end
DoF, DoFE, DoFF, phiB = FEMSei.ConstructCG(1,Grids.Tri())
@show phiB
DoF, DoFE, DoFF, phiB = FEMSei.ConstructCG(2,Grids.Tri())
@show phiB
DoF, DoFE, DoFF, phiB = FEMSei.ConstructCG(3,Grids.Tri())
@show phiB
DoF, DoFE, DoFF, phiB = FEMSei.ConstructCG(4,Grids.Tri())
@show phiB
stop

#Finite elements
DG = FEMSei.DGStruct{FTB}(backend,k,Grid.Type,Grid)
RT = FEMSei.RTStruct{FTB}(backend,k,Grid.Type,Grid)
ND = FEMSei.NDStruct{FTB}(backend,k,Grid.Type,Grid)
CG = FEMSei.CGStruct{FTB}(backend,k+1,Grid.Type,Grid)

ModelFEM = FEMSei.ModelFEM(backend,FTB,ND,RT,CG,DG,Grid,nQuadM,nQuadS,FEMSei.Jacobi!)

pPosS = ModelFEM.pPosS
pPosE = ModelFEM.pPosE
uPosS = ModelFEM.uPosS
uPosE = ModelFEM.uPosE
U = zeros(FTB,ModelFEM.DG.NumG+ModelFEM.RT.NumG)
@views Up = U[pPosS:pPosE]
@views Uu = U[uPosS:uPosE]

FEMSei.InterpolateDG!(Up,DG,FEMSei.Jacobi!,Grid,Grid.Type,Model.InitialProfile)
FEMSei.InterpolateRT!(Uu,RT,FEMSei.Jacobi!,Grid,Grid.Type,nQuad,Model.InitialProfile)

cName = ["h";"Vort";"uS";"vS"]

@show  nAdveVel
FEMSei.TimeStepper(backend,FTB,U,dtau,FEMSei.FcnNonLinShallow!,ModelFEM,Grid,nQuadM,nQuadS,FEMSei.Jacobi!,
  nAdveVel,FileNameOutput,Proc,ProcNumber,nPrint,ref)
