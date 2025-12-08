import CGDycore:
  Thermodynamics, Examples, Parallels, Models, Grids, Surfaces,  Outputs, Integration,  GPUS, DyCore
using MPI
using Base
using CUDA
#using AMDGPU
using Metal
#using oneAPI
using KernelAbstractions
using StaticArrays
using ArgParse
using MPI
using NetCDF
using NCDatasets


# Model
parsed_args = DyCore.parse_commandline()
Problem = parsed_args["Problem"]
Discretization = parsed_args["Discretization"]
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
NumV = parsed_args["NumV"]
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
IntMethod = parsed_args["IntMethod"]
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
#    AMDGPU.allowscalar(false)
  elseif JuliaGPU == "Metal"
    backend = MetalBackend()
    Metal.allowscalar(false)
#  elseif JuliaGPU == "oneAPI"
#    backend = oneAPIBackend()
#    oneAPI.allowscalar(false)
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

KernelAbstractions.synchronize(backend)

OrdPolyZ=1
Parallel = true

# Physical parameters
Phys = DyCore.PhysParameters{FTB}(;ScaleFactor)

#ModelParameters
Model = DyCore.ModelStruct{FTB}()

# Initial conditions
Model.NumV=NumV
Model.NumTr=NumTr
Model.NumThermo = 4
if State == "MoistInternalEnergy" || State == "MoistTotalEnergy"
  Model.NumThermo +=2
elseif State == "IceInternalEnergy"
  Model.NumThermo +=3
end  
Model.Problem=Problem
if ProfRho == ""
  Model.ProfRho = Problem
else
  Model.ProfRho = ProfRho  
end  
if ProfTheta == ""
  Model.ProfTheta = Problem
else
  Model.ProfTheta = ProfTheta  
end  
if ProfVel == ""
  Model.ProfVel = Problem
else
  Model.ProfVel = ProfVel  
end  
Model.ProfpBGrd = ProfpBGrd
Model.ProfRhoBGrd = ProfRhoBGrd
Model.RefProfile = RefProfile
Model.Profile = Profile
Model.RhoPos=1
Model.uPos=2
Model.vPos=3
Model.wPos=4
Model.ThPos=5
Model.RhoTPos  = RhoTPos
Model.RhoVPos  = RhoVPos
Model.RhoCPos  = RhoCPos
Model.RhoIPos  = RhoIPos
Model.RhoRPos  = RhoRPos
Model.TkePos  = TkePos
Model.ModelType = ModelType
Model.HorLimit = HorLimit
Model.Upwind = Upwind
Model.Damping = Damping
Model.StrideDamp = StrideDamp
Model.Relax = Relax
Model.Forcing = Forcing
Model.Coriolis = Coriolis
Model.CoriolisType = CoriolisType
Model.Buoyancy = Buoyancy
Model.Turbulence = Turbulence
Model.VerticalDiffusion = VerticalDiffusion
Model.JacVerticalDiffusion = JacVerticalDiffusion
Model.JacVerticalAdvection = JacVerticalAdvection
Model.VerticalDiffusionMom = VerticalDiffusionMom
Model.Source = Source
Model.Microphysics = Microphysics
Model.Sedimentation = Sedimentation
Model.RelCloud = RelCloud
Model.Rain = Rain
Model.Source = Source
Model.SurfaceFlux = SurfaceFlux
Model.SurfaceFluxMom = SurfaceFluxMom
Model.Thermo = Thermo
Model.Curl = Curl
Model.Stretch = Stretch
Model.StretchType = StretchType
Model.State = State
Model.ModelType = ModelType
Model.HyperVisc = HyperVisc
Model.HyperDCurl = HyperDCurl # =7.e15
Model.HyperDGrad = HyperDGrad # =7.e15
Model.HyperDDiv = HyperDDiv # =7.e15


# Grid
if GridForm == "Cartesian"
  Boundary = Grids.Boundary()
  Boundary.WE = BoundaryWE
  Boundary.SN = BoundarySN
  Boundary.BT = BoundaryBT
  Topography=(TopoS=TopoS,
              H=H,
              P1=P1,
              P2=P2,
              P3=P3,
              P4=P4,
              )
  Grid, CellToProc = Grids.InitGridCart(backend,FTB,OrdPoly,nx,ny,Lx,Ly,x0,y0,Boundary,nz,Model,ParallelCom)
else  
  ns=50
  if RadEarth == 0.0
    RadEarth = Phys.RadEarth
    if ScaleFactor != 0.0
      RadEarth = RadEarth / ScaleFactor
    end
  end
  Grid, CellToProc = Grids.InitGridSphere(backend,FTB,OrdPoly,nz,nPanel,RefineLevel,ns,nLon,nLat,LatB,
    GridType,Decomp,RadEarth,Model,ParallelCom;Discretization=Discretization)
  Topography = (TopoS=TopoS,H=H,Rad=RadEarth)
end  


HeightCG = DyCore.ComputeOrography(backend,FTB,OrdPoly,OrdPolyZ,OrdPrint,H,Topography,Model,
  Phys,CellToProc,Grid,ParallelCom)
OrdPolyS = "$OrdPoly"
nPanelS = "$nPanel"
filename = "TopoO" * OrdPolyS * "P" * nPanelS * ".nc"
ds = NCDataset(filename,"c")
v = defVar(ds,"Height",HeightCG,("DoF","NF"))
close(ds)
