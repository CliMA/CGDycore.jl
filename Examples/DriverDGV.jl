import CGDycore:
  Thermodynamics, Examples, Parallels, Models, Grids, Surfaces,  Outputs, Integration, FiniteElements, DGVertical, GPU, DyCore
using MPI
using Base
using CUDA
using AMDGPU
using Metal
using KernelAbstractions
using StaticArrays
using ArgParse


# Model
parsed_args = DyCore.parse_commandline()
Problem = parsed_args["Problem"]
Discretization = parsed_args["Discretization"]
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
  elseif JuliaGPU == "AMD"
    backend = ROCBackend()
    AMDGPU.allowscalar(false)
  elseif JuliaGPU == "Metal"
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
Param = Examples.Parameters(FTB,Problem)

KernelAbstractions.synchronize(backend)

Parallel = true

# Physical parameters
Phys = DyCore.PhysParameters{FTB}()

#ModelParameters
Model = DyCore.ModelStruct{FTB}()

# Initial conditions
Model.NumV = NumV
Model.NumAux = NumAux
Model.NumTr = NumTr
Model.NumThermo = 4
if State == "MoistInternalEnergy"
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

# Equation
if Equation == "CompressibleShallow"
  Model.Equation = Models.CompressibleShallow()
elseif Equation == "CompressibleDeep"
  Model.Equation = Models.CompressibleDeep()
end

z,zP,dzeta = Grids.AddVerticalGrid(nz,H)

DG1 = FiniteElements.DG1{FTB}(backend,OrdPolyZ,OrdPrintZ)

X = KernelAbstractions.zeros(backend,FTB,DG1.OrdPolyZ+1,nz)
dXdxI = KernelAbstractions.zeros(backend,FTB,DG1.OrdPolyZ+1,nz)
J = KernelAbstractions.zeros(backend,FTB,DG1.OrdPolyZ+1,nz)

Grids.JacobiDG1GPU!(X,dXdxI,J,DG1,z)
NumV = 3
NumAux = 2
U = KernelAbstractions.zeros(backend,FTB,DG1.OrdPolyZ+1,nz,NumV)
Profile = Examples.StratifiedExample()(Param,Phys)
time = 0.0
for iz = 1 : nz
  for K = 1 : DG1.OrdPolyZ + 1
    xS = SVector{3}(0.0,0.0,X[K,iz])
    RhoP,_,_,_,ThP= Profile(xS,time)
    U[K,iz,1] = RhoP
    U[K,iz,3] = RhoP * ThP
  end
end  

F = similar(U)
UNew = similar(U)
CacheU = KernelAbstractions.zeros(backend,FTB,DG1.OrdPolyZ+1,nz,NumAux)
Pressure, dPresdRhoTh, dPresdRho = Models.DryDG()(Phys)

RhoPos = 1
wPos = 2
ThPos = 3
pAuxPos = 1
GPAuxPos = 2

FluxAverage = DGVertical.KennedyGruberGravV()(RhoPos,wPos,ThPos,pAuxPos,GPAuxPos)
RiemannSolver = DGVertical.RiemannLMARSV()(Param,Phys,RhoPos,wPos,ThPos,pAuxPos)

Jac = DGVertical.Jacobian(U,DG1,J)
stop

dtau = 0.1
nIter = 1000
for Iter = 1 : nIter
   @show Iter,sum(abs.(U))   
   @show Iter,sum(abs.(U[:,:,wPos]))   
   DGVertical.FcnGPUVert!(F,U,DG1,X,dXdxI,J,CacheU,Pressure,Phys,FluxAverage,RiemannSolver) 
   @. UNew = U + 1/3 * dtau *F
   DGVertical.FcnGPUVert!(F,UNew,DG1,X,dXdxI,J,CacheU,Pressure,Phys,FluxAverage,RiemannSolver) 
   @. UNew = U + 1/2 * dtau *F
   DGVertical.FcnGPUVert!(F,UNew,DG1,X,dXdxI,J,CacheU,Pressure,Phys,FluxAverage,RiemannSolver) 
   @. U = U + dtau *F
 end  
