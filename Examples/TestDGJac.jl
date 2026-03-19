import CGDycore:
  Parameters, Thermodynamics, Examples, Sources, Parallels, Models, Grids, Surfaces,  Outputs, Integration, FiniteElements, DGSEM, CGSEM, DyCore, IMEXRosenbrock
using MPI
using Base
using CUDA
#using AMDGPU
#using Metal
using KernelAbstractions
#using StaticArrays
using ArgParse


#=
# Model
parsed_args = Parameters.parse_commandline()
Problem = parsed_args["Problem"]
Discretization = parsed_args["Discretization"]
VelocityForm = parsed_args["VelocityForm"]
FluxDG = parsed_args["FluxDG"]
InterfaceFluxDG = parsed_args["InterfaceFluxDG"]
ProfRho = parsed_args["ProfRho"]
ProfTheta = parsed_args["ProfTheta"]
PertTh = parsed_args["PertTh"]
ProfVel = parsed_args["ProfVel"]
ProfpBGrd = parsed_args["ProfpBGrd"]
ProfRhoBGrd = parsed_args["ProfRhoBGrd"]
RhoTPos = parsed_args["RhoTPos"]
RhoVPos = parsed_args["RhoVPos"]
RhoCPos = parsed_args["RhoCPos"]
RhoIPos = parsed_args["RhoIPos"]
RhoRPos = parsed_args["RhoRPos"]
TkePos = parsed_args["TkePos"]
pAuxPos = parsed_args["pAuxPos"]
GPAuxPos = parsed_args["GPAuxPos"]
NumV = parsed_args["NumV"]
NumAux = parsed_args["NumAux"]
NumTr = parsed_args["NumTr"]
HorLimit = parsed_args["HorLimit"]
Upwind = parsed_args["Upwind"]
Damping = parsed_args["Damping"]
Relax = parsed_args["Relax"]
StrideDamp = parsed_args["StrideDamp"]
Forcing = parsed_args["Forcing"]
BoundaryWE = parsed_args["BoundaryWE"]
BoundarySN = parsed_args["BoundarySN"]
BoundaryBT = parsed_args["BoundaryBT"]
Thermo = parsed_args["Thermo"]
State = parsed_args["State"]
RefProfile = parsed_args["RefProfile"]
Profile = parsed_args["Profile"]
Curl = parsed_args["Curl"]
ModelType = parsed_args["ModelType"]
Equation = parsed_args["Equation"]
Microphysics = parsed_args["Microphysics"]
Sedimentation = parsed_args["Sedimentation"]
TypeMicrophysics = parsed_args["TypeMicrophysics"]
RelCloud = parsed_args["RelCloud"]
Rain = parsed_args["Rain"]
#Orography
TopoS = parsed_args["TopoS"]
P1 = parsed_args["P1"]
P2 = parsed_args["P2"]
P3 = parsed_args["P3"]
P4 = parsed_args["P4"]

# Parallel
Decomp = parsed_args["Decomp"]
SimDays = parsed_args["SimDays"]
SimHours = parsed_args["SimHours"]
SimMinutes = parsed_args["SimMinutes"]
SimSeconds = parsed_args["SimSeconds"]
SimTime = parsed_args["SimTime"]
dtau = parsed_args["dtau"]
dtauSmall = parsed_args["dtauSmall"]
IntMethod = parsed_args["IntMethod"]
IntMethodFast = parsed_args["IntMethodFast"]
Table = parsed_args["Table"]
GridForm = parsed_args["GridForm"]
GridType = parsed_args["GridType"]
AdaptGridType = parsed_args["AdaptGridType"]
RadEarth = parsed_args["RadEarth"]
ScaleFactor = parsed_args["ScaleFactor"]
Coriolis = parsed_args["Coriolis"]
CoriolisType = parsed_args["CoriolisType"]
Buoyancy = parsed_args["Buoyancy"]
Turbulence = parsed_args["Turbulence"]
Source = parsed_args["Source"]
VerticalDiffusion = parsed_args["VerticalDiffusion"]
JacVerticalDiffusion = parsed_args["JacVerticalDiffusion"]
JacVerticalAdvection = parsed_args["JacVerticalAdvection"]
VerticalDiffusionMom = parsed_args["VerticalDiffusionMom"]
SurfaceFlux = parsed_args["SurfaceFlux"]
SurfaceFluxMom = parsed_args["SurfaceFluxMom"]
SurfaceScheme = parsed_args["SurfaceScheme"]
# Grid
nx = parsed_args["nx"]
ny = parsed_args["ny"]
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
Lx = parsed_args["Lx"]
Ly = parsed_args["Ly"]
x0 = parsed_args["x0"]
# CG Method
y0 = parsed_args["y0"]
OrdPoly = parsed_args["OrdPoly"]
OrdPolyZ = parsed_args["OrdPolyZ"]
# Viscosity
HyperVisc = parsed_args["HyperVisc"]
HyperDCurl = parsed_args["HyperDCurl"]
HyperDGrad = parsed_args["HyperDGrad"]
HyperDDiv = parsed_args["HyperDDiv"]
# Output
OrdPrint = parsed_args["OrdPrint"]
OrdPrintZ = parsed_args["OrdPrintZ"]
PrintDays = parsed_args["PrintDays"]
PrintHours = parsed_args["PrintHours"]
PrintMinutes = parsed_args["PrintMinutes"]
PrintSeconds = parsed_args["PrintSeconds"]
PrintTime = parsed_args["PrintTime"]
PrintStartTime = parsed_args["PrintStartTime"]
vtkFileName = parsed_args["vtkFileName"]
Flat = parsed_args["Flat"]
# Device
Device = parsed_args["Device"]
GPUType = parsed_args["GPUType"]
FloatTypeBackend = parsed_args["FloatTypeBackend"]
NumberThreadGPU = parsed_args["NumberThreadGPU"]
NumberThreadTriGPU = parsed_args["NumberThreadTriGPU"]
# Examples
aC = parsed_args["LengthOfAgnesiHill"]
=#

MPI.Init()
comm = MPI.COMM_WORLD
Proc = MPI.Comm_rank(comm) + 1
ProcNumber = MPI.Comm_size(comm)
ParallelCom = DyCore.ParallelComStruct()
ParallelCom.Proc = Proc
ParallelCom.ProcNumber  = ProcNumber

JuliaDevice = get(ENV, "JuliaDevice", "CPU")
JuliaGPU = get(ENV, "JuliaGPU", "CUDA")
machine = get(ENV, "machine", "")

if JuliaDevice == "CPU"
  backend = CPU()
elseif JuliaDevice == "GPU"
  if JuliaGPU == "CUDA"
    backend = CUDABackend()
    CUDA.allowscalar(false)
    if machine == "levante" || machine == "derecho"
    else
      CUDA.device!(Proc-1)
    end
#  elseif JuliaGPU == "AMD"
#    backend = ROCBackend()
#    AMDGPUS.allowscalar(false)
#  elseif JuliaGPU == "Metal"
#    backend = MetalBackend()
#    Metal.allowscalar(true)
  end
else
  backend = CPU()
end
FloatTypeBackend = "Float64"
if FloatTypeBackend == "Float64"
  FTB = Float64
elseif FloatTypeBackend == "Float32"
  FTB = Float32
else
  @show "False FloatTypeBackend"
  stop
end

Problem = "BaroWaveDrySphere"

Param = Examples.Parameters(FTB,Problem)

KernelAbstractions.synchronize(backend)

Parallel = true

# Physical parameters
Phys = DyCore.PhysParameters{FTB}()

RefineLevel = 6
nz = 4
nQuad = 3
nQuadM = 3 #2
nQuadS = 3 #3
Decomp = "EqualArea"
nLat = 0
nLon = 0
LatB = 0.0

#Quad
GridType = "CubedSphere"
nPanel = 2
#GridType = "HealPix"
ns = 57
OrdPoly = 3
OrdPolyZ = 7

#Grid construction
RadEarth = Phys.RadEarth
Discretization = "DiscretizationDG"
#ModelParameters
Model = DyCore.ModelStruct{FTB}()

Grid, CellToProc = Grids.InitGridSphere(backend,FTB,OrdPoly,nz,nPanel,RefineLevel,ns,nLon,nLat,LatB,
    GridType,Decomp,RadEarth,Model,ParallelCom;Discretization=Discretization,ChangeOrient=2)

OrdPrint = 1
OrdPrintZ = 1
TimeStepper = DyCore.TimeStepperStruct{FTB}(backend)

Output = DyCore.OutputStruct()
DoF = (OrdPoly + 1) * (OrdPoly + 1)
NumV = 5
NumTr = 0
Global = DyCore.GlobalStruct{FTB}(backend,Grid,Model,TimeStepper,ParallelCom,Output,DoF,nz,
    NumV,NumTr)
DG = FiniteElements.DGQuad{FTB}(backend,OrdPoly,OrdPolyZ,OrdPrint,OrdPrintZ,Global.Grid,ParallelCom.Proc)

dSdS,dSdM,dMdS,dMdM = DGSEM.InitJacDG(DG,nz,Param)

U = ones(OrdPolyZ+1,nz,DG.NumI,NumV)
fac = 1.0
Invfac = 1.0 / fac
dz = ones(nz,DG.NumI)
Jac,JacLU = DGSEM.JacDG(U,DG,fac,dSdS,dSdM,dMdS,dMdM,dz,Phys)

JCache = DGSEM.JacDGVert{FTB}(backend,OrdPolyZ+1,nz,DG.NumI)

DGSEM.FillJacDGVert!(JCache,U,DG,dz,Invfac,Phys)

p = DGSEM.Permutation(OrdPolyZ+1,nz)

J = Jac[1]

JP = J[p,p]

nb = OrdPolyZ + 1 + 2 * (OrdPolyZ - 1)

A = JP[1:nz*nb,1:nz*nb]
B = JP[1:nz*nb,nz*nb+1:end]
C = JP[nz*nb+1:end,1:nz*nb]
D = JP[nz*nb+1:end,nz*nb+1:end]
