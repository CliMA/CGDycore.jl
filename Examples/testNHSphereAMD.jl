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
Geos = parsed_args["Geos"]
Relax = parsed_args["Relax"]
StrideDamp = parsed_args["StrideDamp"]
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
StartAverageDays = parsed_args["StartAverageDays"]
dtau = parsed_args["dtau"]
IntMethod = parsed_args["IntMethod"]
Table = parsed_args["Table"]
# Grid
nz = parsed_args["nz"]
nPanel = parsed_args["nPanel"]
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
# Output
PrintDays = parsed_args["PrintDays"]
PrintHours = parsed_args["PrintHours"]
PrintMinutes = parsed_args["PrintMinutes"]
PrintSeconds = parsed_args["PrintSeconds"]
PrintStartTime = parsed_args["PrintStartTime"]
Flat = parsed_args["Flat"]

# Device
Device = parsed_args["Device"]
GPUType = parsed_args["GPUType"]
FloatTypeBackend = parsed_args["FloatTypeBackend"]
NumberThreadGPU = parsed_args["NumberThreadGPU"]

MPI.Init()

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



# Physical parameters
Phys = DyCore.PhysParameters{FTB}()

#ModelParameters
Model = DyCore.ModelStruct{FTB}()

# Initial conditions
Model.Equation=Equation
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
Model.PertTh = PertTh
if ProfVel == ""
  Model.ProfVel = Problem
else
  Model.ProfVel = ProfVel  
end  
if ProfVelGeo == ""
  Model.ProfVelGeo = Problem
else
  Model.ProfVelGeo = ProfVelGeo  
end  
Model.RefProfile = RefProfile
Model.ProfpBGrd = ProfpBGrd
Model.ProfRhoBGrd = ProfRhoBGrd
Model.RhoPos = 1
Model.uPos = 2
Model.vPos = 3
Model.wPos = 4
Model.ThPos = 5
Model.RhoVPos  = RhoVPos
Model.RhoCPos  = RhoCPos
Model.RhoIPos  = RhoIPos
Model.RhoRPos  = RhoRPos
Model.HorLimit = HorLimit
Model.Upwind = Upwind
Model.Damping = Damping
Model.Geos = Geos
Model.StrideDamp = StrideDamp
Model.Relax = Relax
Model.Coriolis = Coriolis
Model.CoriolisType = CoriolisType
Model.Buoyancy = Buoyancy
Model.VerticalDiffusion = VerticalDiffusion
Model.JacVerticalDiffusion = JacVerticalDiffusion
Model.JacVerticalAdvection = JacVerticalAdvection
Model.Source = Source
Model.Forcing = Forcing
Model.Microphysics = Microphysics
Model.TypeMicrophysics = TypeMicrophysics
Model.RelCloud = RelCloud
Model.Rain = Rain
Model.SurfaceFlux = SurfaceFlux
Model.SurfaceFluxMom = SurfaceFluxMom
Model.Thermo = Thermo
Model.Curl = Curl
Model.ModelType = ModelType
Model.Stretch = Stretch
Model.StretchType = StretchType
Model.HyperVisc = HyperVisc
Model.HyperDCurl = HyperDCurl
Model.HyperDGrad = HyperDGrad
Model.HyperDRhoDiv = HyperDRhoDiv
Model.HyperDDiv = HyperDDiv

OrdPolyZ = 1
if RadEarth == 0.0
  RadEarth = Phys.RadEarth
end

Topography = (TopoS=TopoS,H=H,Rad=RadEarth)

@show "InitSphere"
(CG, Metric, Exchange, Global) = DyCore.InitSphere(backend,FTB,OrdPoly,OrdPolyZ,nz,nPanel,H,
  GridType,Topography,Decomp,Model,Phys,RadEarth)

# Initial values
if Problem == "Galewski"
  Profile = Examples.GalewskiExample()(Param,Phys)
elseif Problem == "BaroWaveDrySphere"
  Profile = Examples.BaroWaveExample()(Param,Phys)
elseif Problem == "HeldSuarezDrySphere"
  Model.InitialProfile, Model.Force = Examples.HeldSuarezDryExample()(Param,Phys)
elseif Problem == "HeldSuarezMoistSphere"
  Profile, Force, Eddy = Examples.HeldSuarezMoistExample()(Param,Phys)
  Model.InitialProfile = Profile
  Model.Force = Force
  Model.Eddy = Eddy
end  

U = GPU.InitialConditions(backend,FTB,CG,Metric,Phys,Global,Profile,Param)

# Pressure
if Equation == "Compressible"
  Pressure = Models.Compressible()(Phys)
  Model.Pressure = Pressure
elseif Equation == "CompressibleMoist"
  Pressure = Models.CompressibleMoist()(Phys,Model.RhoPos,Model.ThPos,
    Model.RhoVPos+NumV,Model.RhoCPos+NumV)
  Model.Pressure = Pressure
end  
# Microphysics
if Microphysics
  if TypeMicrophysics == "SimpleMicrophysics"
    MicrophysicsSource  = Models.SimpleMicrophysics()(Phys,Model.RhoPos,Model.ThPos,
      Model.RhoVPos+NumV,Model.RhoCPos+NumV,Model.RelCloud,Model.Rain)
    Model.MicrophysicsSource = MicrophysicsSource
  else
    @show "False Type Microphysics"  
  end
end  

# Surface flux
if Model.SurfaceFlux || Model.VerticalDiffusion
  if Problem == "HeldSuarezMoistSphere"
    SurfaceValues, SurfaceData = GPU.HeldSuarezMoistSurface()(Phys,Param,Model.uPos,Model.vPos,Model.wPos)
    Model.SurfaceValues = SurfaceValues
    Model.SurfaceData = SurfaceData
  end  
end

# Output
Global.Output.vtkFileName = string(Problem*"_")
Global.Output.vtk = 0
Global.Output.Flat = Flat
Global.Output.nPanel = nPanel
Global.Output.RadPrint = H
Global.Output.H = H
if ModelType == "VectorInvariant" || ModelType == "Advection"
  if Model.Equation == "Compressible"  
    Global.Output.cNames = [
      "Rho",
      "u",
      "v",
      "wB",
      "Th",
#     "Vort",
#     "Pres",
      ]
  elseif Model.Equation == "CompressibleMoist"  
    Global.Output.cNames = [
      "Rho",
      "u",
      "v",
      "wB",
      "Th",
#     "Vort",
#     "Pres",
      "Tr1",
      "Tr2",
      ]
  elseif Model.Equation == "Shallow"  
    Global.Output.cNames = [
      "Rho",
#     "u",
#     "v",
      "Th",
      "Vort",
#     "Pres",
      ]
  end  
elseif ModelType == "Conservative"
  Global.Output.cNames = [
    "Rho",
    "Rhou",
    "Rhov",
    "w",
    "Th",
#   "Vort",
#   "Pres",
    ]
end

Global.Output.PrintDays = PrintDays
Global.Output.PrintHours = PrintHours
Global.Output.PrintMinutes = PrintMinutes
Global.Output.PrintSeconds = PrintSeconds
Global.Output.StartAverageDays = StartAverageDays
Global.Output.PrintStartTime = PrintStartTime
Global.Output.OrdPrint = CG.OrdPoly

Global.vtkCache = Outputs.vtkStruct{FTB}(backend,Global.Output.OrdPrint,Grids.TransSphereX!,CG,Metric,Global)

# TimeStepper
time=[0.0]
Global.TimeStepper.IntMethod = IntMethod
Global.TimeStepper.Table = Table
Global.TimeStepper.dtau = dtau
Global.TimeStepper.SimDays = SimDays
Global.TimeStepper.SimMinutes = SimMinutes
Global.TimeStepper.SimSeconds = SimSeconds
if ModelType == "VectorInvariant" || ModelType == "Advection"
  DiscType = Val(:VectorInvariant)  
elseif ModelType == "Conservative"
  DiscType = Val(:Conservative)  
end  


if Device == "CPU"  || Device == "GPU"
  Global.ParallelCom.NumberThreadGPU = NumberThreadGPU   
  nT = max(7 + NumTr, NumV + NumTr)
  Parallels.InitExchangeData3D(backend,FTB,nz,nT,Exchange)
  @show "vor Timestepper GPU"
  Integration.TimeStepper!(U,GPU.FcnGPUAMD!,GPU.FcnPrepareGPU!,DyCore.JacSchurGPU!,
    Grids.TransSphereX,CG,Metric,Phys,Exchange,Global,Param,DiscType)
else
  nT = max(7 + NumTr, NumV + NumTr)
  Parallels.InitExchangeData3D(backend,FTB,nz,nT,Exchange)
  Integration.TimeStepper!(U,DyCore.Fcn!,DyCore.FcnPrepare!,DyCore.JacSchurGPU!,
    Grids.TransSphereX,CG,Metric,Phys,Exchange,Global,Param,DiscType)
end
