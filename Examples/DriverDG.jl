import CGDycore:
  Thermodynamics, Examples, Parallels, Models, Grids, Surfaces,  Outputs, Integration, FiniteElements, DGSEM, GPUS, DyCore
using MPI
using Base
using CUDA
using AMDGPU
using Metal
using KernelAbstractions
#using StaticArrays
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
dtauSmall = parsed_args["dtauSmall"]
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
# Examples
aC = parsed_args["LengthOfAgnesiHill"]

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
    AMDGPUS.allowscalar(false)
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
  Grid, CellToProc = Grids.InitGridCart(backend,FTB,OrdPoly,nx,ny,Lx,Ly,x0,y0,Boundary,nz,Model,ParallelCom;
    Discretization=Discretization,GridType=GridType,ChangeOrient=2)
  Trans = Outputs.TransCartX!
else  
  if RadEarth == 0.0
    RadEarth = Phys.RadEarth
    if ScaleFactor != 0.0
      RadEarth = RadEarth / ScaleFactor
    end
  end
  Grid, CellToProc = Grids.InitGridSphere(backend,FTB,OrdPoly,nz,nPanel,RefineLevel,ns,nLon,nLat,LatB,
    GridType,Decomp,RadEarth,Model,ParallelCom;Discretization=Discretization,ChangeOrient=2)
  Topography = (TopoS=TopoS,H=H,Rad=RadEarth)
  Trans = Outputs.TransSphereX!
end  


#Topography
if TopoS == "AgnesiHill"
  TopoProfile = Examples.AgnesiHill{FTB}(;aC=aC)()
elseif TopoS == "SchaerHill"
  TopoProfile = Examples.SchaerHill()()
elseif TopoS == "BaroWaveHill"
  TopoProfile = Examples.BaroWaveHill()()
elseif TopoS == "SchaerSphereCircle"
  TopoProfile = Examples.SchaerSphereCircle()(Param,Phys)
else
  TopoProfile = Examples.Flat()()  
end  

Grid.AdaptGrid = Grids.AdaptGrid(FTB,AdaptGridType,FTB(H))
DGMethod = "Kubatko2LGL"

if GridForm == "Cartesian"
  if ParallelCom.Proc == 1
    @show "InitCart"
  end
  (DG, Metric, Exchange, Global) = DyCore.InitCartDG(backend,FTB,OrdPoly,OrdPolyZ,DGMethod,
    OrdPrint,OrdPrintZ,H,Topography,Model,
    Phys,TopoProfile,CellToProc,Grid,ParallelCom)
else
  (DG, Metric, Exchange, Global) = DyCore.InitSphereDG(backend,FTB,OrdPoly,OrdPolyZ,DGMethod,
    OrdPrint,OrdPrintZ,H,Topography,Model,
    Phys,TopoProfile,CellToProc,Grid,ParallelCom)
end


# Initial values
Examples.InitialProfile!(backend,FTB,Model,Problem,Param,Phys)
U = GPUS.InitialConditions(backend,FTB,DG,Metric,Phys,Global,Model.InitialProfile,Param)

pAuxPos = 1
GPAuxPos = 2

if InterfaceFluxDG == "RiemannLMARS"
  RiemannSolver = DGSEM.RiemannLMARS()(Param,Phys,Model.RhoPos,Model.uPos,Model.vPos,Model.wPos,Model.ThPos,1)
  Model.RiemannSolver = RiemannSolver
elseif InterfaceFluxDG == "RiemannExLMARS"
  RiemannSolver = DGSEM.RiemannExLMARS()(Param,Phys,Model.RhoPos,Model.uPos,Model.vPos,Model.wPos,Model.ThPos,1)
  Model.RiemannSolver = RiemannSolver
elseif InterfaceFluxDG == "RiemannExnerLMARS"
  RiemannSolver = DGSEM.RiemannExnerLMARS()(Param,Phys,Model.RhoPos,Model.uPos,Model.vPos,Model.wPos,Model.ThPos,1)
  Model.RiemannSolver = RiemannSolver
elseif InterfaceFluxDG == "ArtianoEnergyStable"
  RiemannSolver = DGSEM.ArtianoEnergyStable()(Param,Phys,Model.RhoPos,Model.uPos,Model.vPos,Model.wPos,Model.ThPos,1)
  Model.RiemannSolver = RiemannSolver
elseif InterfaceFluxDG == "RiemannExPLMARS"
  RiemannSolver = DGSEM.RiemannExPLMARS()(Param,Phys,Model.RhoPos,Model.uPos,Model.vPos,Model.wPos,Model.ThPos,1)
  Model.RiemannSolver = RiemannSolver
elseif InterfaceFluxDG == "RiemannBoussinesqLMARS"  
  RiemannSolver = DGSEM.RiemannBoussinesqLMARS()(Param,Model.RhoPos,Model.uPos,Model.vPos,Model.wPos,Model.ThPos)
  Model.RiemannSolver = RiemannSolver
end  

NonConservativeFlux = DGSEM.BuoyancyFlux()(Model.RhoPos,GPAuxPos)
Model.NonConservativeFlux = NonConservativeFlux

if FluxDG == "KennedyGruber"
  Model.FluxAverage = DGSEM.KennedyGruber()(Model.RhoPos,Model.uPos,Model.vPos,Model.wPos,Model.ThPos,1)
elseif FluxDG == "KennedyGruberGrav"  
  Model.FluxAverage = DGSEM.KennedyGruberGrav()(Model.RhoPos,Model.uPos,Model.vPos,Model.wPos,
    Model.ThPos,pAuxPos,GPAuxPos)
elseif FluxDG == "ArtianoExner"  
  Model.FluxAverage = DGSEM.ArtianoExner()(Model.RhoPos,Model.uPos,Model.vPos,Model.wPos,
    Model.ThPos,pAuxPos,GPAuxPos,Phys)
elseif FluxDG == "KennedyGruberExPGrav"  
  Model.FluxAverage = DGSEM.ArtianoExPGrav()(Model.RhoPos,Model.uPos,Model.vPos,Model.wPos,
    Model.ThPos,pAuxPos,GPAuxPos,Phys)
elseif FluxDG == "ArtianoGrav"  
  Model.FluxAverage = DGSEM.ArtianoGrav()(Model.RhoPos,Model.uPos,Model.vPos,Model.wPos,
  Model.ThPos,pAuxPos,GPAuxPos,Phys)
elseif FluxDG == "ArtianoExGrav"  
  Model.FluxAverage = DGSEM.ArtianoExGrav()(Model.RhoPos,Model.uPos,Model.vPos,Model.wPos,
  Model.ThPos,pAuxPos,GPAuxPos,Phys)
elseif FluxDG == "LinearBoussinesqFlux"
  Model.Flux = DGSEM.LinearBoussinesqFlux()(Param,Model.RhoPos,Model.uPos,Model.vPos,Model.wPos,Model.ThPos)
end  

if ModelType == "Boussinesq"
  Model.BuoyancyFun = GPUS.BuoyancyBoussinesq()(Param,Model.wPos,Model.ThPos)
end  
    

# Pressure
if State == "Dry"
  Pressure, dPresdRhoTh, dPresdRho = Models.DryDG()(Phys)
  Model.Pressure = Pressure
  Model.dPresdRhoTh = dPresdRhoTh
  Model.dPresdRho = dPresdRho
elseif State  == "ShallowWater"  
  Model.Pressure = Models.ShallowWaterStateDG()(Phys)
end

if Damping
  Damp = GPUS.DampingW()(FTB(H),FTB(StrideDamp),FTB(Relax),Model.wPos)
  Model.Damp = Damp
end
if Grid.Form == "Sphere"
  Model.GeoPotential = GPUS.GeoPotentialDeep()(Phys)
else
  Model.GeoPotential = GPUS.GeoPotentialCart()(Phys)
end  
    



Global.ParallelCom.NumberThreadGPU = NumberThreadGPU

if ModelType == "Conservative"
  Global.Output.cNames = [
    "Rho",
    "Rhou",
    "Rhov",
    "Rhow",
    "RhoTh",
#   "w",
#   "Th",
#   "Vort",
    ]
elseif ModelType == "Boussinesq"
  Global.Output.cNames = [
    "Rho",
    "u",
    "wDG",
    "BDG",
    ]
end

Global.Output.Flat = Flat
Global.Output.PrintDays = PrintDays
Global.Output.PrintHours = PrintHours
Global.Output.PrintMinutes = PrintMinutes
Global.Output.PrintSeconds = PrintSeconds
Global.Output.PrintTime = PrintTime
Global.Output.PrintStartTime = PrintStartTime
if OrdPrint <  0
  Global.Output.OrdPrint = DG.OrdPoly
  Global.Output.OrdPrintZ = DG.OrdPolyZ
else
  Global.Output.OrdPrint = OrdPrint
  Global.Output.OrdPrintZ = OrdPrintZ
end
Global.Output.nPanel = nPanel
Global.Output.dTol = pi/30
Global.Output.vtkFileName = vtkFileName
Global.vtkCache = Outputs.vtkStruct{FTB}(backend,Global.Output.OrdPrint,Global.Output.OrdPrintZ,Trans,DG,Metric,Global)

Parallels.InitExchangeData3D(backend,FTB,nz*(OrdPolyZ+1),NumV+NumAux+1,Exchange)


# Simulation time
GridLengthMin,GridLengthMax = Grids.GridLength(Grid)
if dtau == 0.0
  dtau = GridLengthMin / Param.cS / sqrt(2)  / (OrdPoly + 1)^1.5
end  
if nz > 1 && IntMethod == "RungeKutta"
  dtau = min(Grid.H / nz / Param.cS / (OrdPolyZ + 1)^1.5, dtau)
end  
EndTime = SimTime + 3600*24*SimDays + 3600 * SimHours + 60 * SimMinutes + SimSeconds
IterTime::Int = round(EndTime / dtau)
dtau = EndTime / IterTime
PrintT = PrintTime + 3600*24*PrintDays + 3600 * PrintHours + 60 * PrintMinutes + PrintSeconds
nPrint::Int = ceil(PrintT/dtau)
EndTime = IterTime * dtau / 3600

if Proc == 1
@show GridLengthMin,GridLengthMax
@show dtau
@show EndTime
@show IterTime
@show nPrint
end

if IntMethod == "Rosenbrock"
  Ros = Integration.RosenbrockStruct{FTB}(Table)
  DGSEM.Rosenbrock(Ros,U,DGSEM.FcnGPUSplit!,dtau,IterTime,nPrint,DG,Exchange,Metric,
    Trans,Phys,Param,Grid,Global,Grid.Type)
elseif IntMethod == "MIS"
  Ros = Integration.RosenbrockStruct{FTB}(Table)
  Mis = DGSEM.MISStruct{FTB}("MISRK4")
DGSEM.MIS_Method(Ros,Mis,U,DGSEM.FcnGPUSplitSlow!,DGSEM.FcnGPUSplitFast!,dtauSmall,dtau,IterTime,nPrint,DG,Exchange,Metric,Trans,Phys,Param,Grid,Global)
elseif IntMethod == "RungeKutta"    
  DGSEM.RK3(U,DGSEM.FcnGPUSplit!,dtau,IterTime,nPrint,DG,Exchange,Metric,
    Trans,Phys,Grid,Global)
elseif IntMethod == "RungeKuttaNonConservative"    
  DGSEM.RK3(U,DGSEM.FcnGPUNonConservativeSplit!,dtau,IterTime,nPrint,DG,Exchange,Metric,
    Trans,Phys,Grid,Global)
end  
MPI.Finalize()
