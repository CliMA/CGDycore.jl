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
using MPI
# Model
parsed_args = DyCore.parse_commandline()
Problem = parsed_args["Problem"]
ProfRho = parsed_args["ProfRho"]
ProfTheta = parsed_args["ProfTheta"]
ProfVel = parsed_args["ProfVel"]
ProfVelW = parsed_args["ProfVelW"]
ProfpBGrd = parsed_args["ProfpBGrd"]
ProfRhoBGrd = parsed_args["ProfRhoBGrd"]
ProfTr = parsed_args["ProfTr"]
HorLimit = parsed_args["HorLimit"]
Upwind = parsed_args["Upwind"]
Damping = parsed_args["Damping"]
Relax = parsed_args["Relax"]
StrideDamp = parsed_args["StrideDamp"]
NumV = parsed_args["NumV"]
NumTr = parsed_args["NumTr"]
BoundaryWE = parsed_args["BoundaryWE"]
BoundarySN = parsed_args["BoundarySN"]
BoundaryBT = parsed_args["BoundaryBT"]
Thermo = parsed_args["Thermo"]
RefProfile = parsed_args["RefProfile"]
Profile = parsed_args["Profile"]
Curl = parsed_args["Curl"]
ModelType = parsed_args["ModelType"]
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
GridType = parsed_args["GridType"]
Coriolis = parsed_args["Coriolis"]
CoriolisType = parsed_args["CoriolisType"]
Microphysics = parsed_args["Microphysics"]
Source = parsed_args["Source"]
VerticalDiffusion = parsed_args["VerticalDiffusion"]
VerticalDiffusionMom = parsed_args["VerticalDiffusionMom"]
SurfaceFlux = parsed_args["SurfaceFlux"]
SurfaceFluxMom = parsed_args["SurfaceFluxMom"]
# Grid
RadEarth = parsed_args["RadEarth"]
nPanel = parsed_args["nPanel"]
nz = parsed_args["nz"]
H = parsed_args["H"]
Stretch = parsed_args["Stretch"]
StretchType = parsed_args["StretchType"]
OrdPoly = parsed_args["OrdPoly"]
Lx = parsed_args["Lx"]
Ly = parsed_args["Ly"]
x0 = parsed_args["x0"]
y0 = parsed_args["y0"]
# Viscosity
HyperVisc = parsed_args["HyperVisc"]
HyperDCurl = parsed_args["HyperDCurl"]
HyperDGrad = parsed_args["HyperDGrad"]
HyperDDiv = parsed_args["HyperDDiv"]
# Output
PrintDays = parsed_args["PrintDays"]
PrintHours = parsed_args["PrintHours"]
PrintMinutes = parsed_args["PrintMinutes"]
PrintSeconds = parsed_args["PrintSeconds"]
PrintTime = parsed_args["PrintTime"]
Flat = parsed_args["Flat"]
# Device
Device = parsed_args["Device"]
GPUType = parsed_args["GPUType"]
FloatTypeBackend = parsed_args["FloatTypeBackend"]

if Device == "CPU"
  backend = CPU()
elseif Device == "GPU"
  if GPUType == "CUDA"
    backend = CUDABackend()
    CUDA.allowscalar(false)
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

Param = Examples.Parameters(FTB,Problem)

KernelAbstractions.synchronize(backend)

MPI.Init()

OrdPolyZ=1
Parallel = true

# Physical parameters
Phys=DyCore.PhysParameters{FTB}()

#ModelParameters
Model = DyCore.ModelStruct{FTB}()
# Initial conditions
Model.Equation="Compressible"
Model.NumV=NumV
Model.NumTr=NumTr
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
if ProfVelW == ""
  Model.ProfVelW = Problem
else
  Model.ProfVelW = ProfVelW  
end  
Model.ProfpBGrd = ProfpBGrd
Model.ProfRhoBGrd = ProfRhoBGrd
Model.ProfTr = ProfTr
Model.RefProfile = RefProfile
Model.Profile = Profile
Model.RhoPos=1
Model.uPos=2
Model.vPos=3
Model.wPos=4
Model.ThPos=5
Model.HorLimit = HorLimit
Model.Upwind = Upwind
Model.Damping = Damping
Model.StrideDamp = StrideDamp
Model.Relax = Relax
Model.Coriolis = Coriolis
Model.CoriolisType = CoriolisType
Model.VerticalDiffusion = VerticalDiffusion
Model.VerticalDiffusionMom = VerticalDiffusionMom
Model.Source = Source
Model.Microphysics = Microphysics
Model.Source = Source
Model.SurfaceFlux = SurfaceFlux
Model.SurfaceFluxMom = SurfaceFluxMom
Model.Thermo = Thermo
Model.Curl = Curl
Model.Stretch = Stretch
Model.StretchType = StretchType
Model.ModelType = ModelType
Model.HyperVisc = HyperVisc
Model.HyperDCurl = HyperDCurl # =7.e15
Model.HyperDGrad = HyperDGrad # =7.e15
Model.HyperDDiv = HyperDDiv # =7.e15




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

# Grid
if RadEarth == 0.0
  RadEarth = Phys.RadEarth
end
Topography=(TopoS=TopoS,H=H,Rad=RadEarth)

OrdPolyZ = 1
@show "vor InitSphere"
(CG,Metric,Exchange,Global) = DyCore.InitSphere(backend,FTB,OrdPoly,OrdPolyZ,nz,nPanel,H,GridType,
  Topography,Decomp,Model,Phys,RadEarth)

if Problem == "AdvectionDCMIP"
  Profile = Examples.AdvectionSphereDCMIP()(Param,Phys)
elseif Problem == "AdvectionSphereSlottedCylinder"
  Profile = Examples.DivergentSphereExample()(Param,Phys)
end  

U = GPU.InitialConditionsAdvection(backend,FTB,CG,Metric,Phys,Global,Profile,Param)

# Output
  Global.Output.vtkFileName=string(Problem*"_")
  Global.Output.vtk=0
  Global.Output.Flat = Flat
  Global.Output.nPanel = nPanel
  Global.Output.H=H
  Global.Output.cNames = [
      "Rho",
      "u",
      "v",
      "Tr1",
      ]
  Global.Output.PrintDays = PrintDays
  Global.Output.PrintSeconds = PrintSeconds
  Global.Output.PrintTime = PrintTime
  Global.Output.PrintStartTime = 0
  Global.Output.OrdPrint=CG.OrdPoly
  @show "vor Global.vtkCache"
  Global.vtkCache = Outputs.vtkStruct{FTB}(backend,Global.Output.OrdPrint,Grids.TransSphereX!,CG,Metric,Global)


  # TimeStepper
  time=[0.0]
  Global.TimeStepper.IntMethod = IntMethod
  Global.TimeStepper.Table = Table
  Global.TimeStepper.dtau = dtau
  Global.TimeStepper.SimDays = SimDays
  Global.TimeStepper.SimHours = SimHours
  Global.TimeStepper.SimMinutes = SimMinutes
  Global.TimeStepper.SimSeconds = SimSeconds
  Global.TimeStepper.SimTime = SimTime

  nT = NumV + NumTr
  Parallels.InitExchangeData3D(backend,FTB,nz,nT,Exchange)  
  @show "vor TimeStepperGPUAdvection!"
  if Device == "CPU" || Device == "GPU" 
    Integration.TimeStepperGPUAdvection!(U,GPU.FcnAdvectionGPU!,Grids.TransSphereX,CG,Metric,Phys,Exchange,Global,Param,Profile)
  else
    Integration.TimeStepperAdvection!(U,DyCore.FcnTracer!,Grids.TransSphereX,CG,Metric,Phys,Exchange,Global,Param,Profile)
  end
