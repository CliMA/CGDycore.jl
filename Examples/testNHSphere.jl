import CGDycore:
  Examples, Parallels, Grids, FiniteElements, Surfaces, Models, Outputs, Integration,  GPU, DyCore
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
ProfRho = parsed_args["ProfRho"]
ProfTheta = parsed_args["ProfTheta"]
PertTh = parsed_args["PertTh"]
ProfVel = parsed_args["ProfVel"]
ProfVelGeo = parsed_args["ProfVelGeo"]
RhoVPos = parsed_args["RhoVPos"]
RhoCPos = parsed_args["RhoCPos"]
RhoIPos = parsed_args["RhoIPos"]
RhoRPos = parsed_args["RhoRPos"]
TkePos = parsed_args["TkePos"]
HorLimit = parsed_args["HorLimit"]
Upwind = parsed_args["Upwind"]
Damping = parsed_args["Damping"]
Relax = parsed_args["Relax"]
StrideDamp = parsed_args["StrideDamp"]
Geos = parsed_args["Geos"]
Coriolis = parsed_args["Coriolis"]
Gravitation = parsed_args["Coriolis"]
CoriolisType = parsed_args["CoriolisType"]
Buoyancy = parsed_args["Buoyancy"]
Turbulence = parsed_args["Turbulence"]
Equation = parsed_args["Equation"]
Thermo = parsed_args["Thermo"]
State = parsed_args["State"]
RefProfile = parsed_args["RefProfile"]
ProfpBGrd = parsed_args["ProfpBGrd"]
ProfRhoBGrd = parsed_args["ProfRhoBGrd"]
EDMF = parsed_args["EDMF"]
NDEDMF = parsed_args["NDEDMF"]
Microphysics = parsed_args["Microphysics"]
TypeMicrophysics = parsed_args["TypeMicrophysics"]
RelCloud = parsed_args["RelCloud"]
Rain = parsed_args["Rain"]
Source = parsed_args["Source"]
Forcing = parsed_args["Forcing"]
VerticalDiffusion = parsed_args["VerticalDiffusion"]
VerticalDiffusionMom = parsed_args["VerticalDiffusionMom"]
JacVerticalDiffusion = parsed_args["JacVerticalDiffusion"]
JacVerticalAdvection = parsed_args["JacVerticalAdvection"]
SurfaceFlux = parsed_args["SurfaceFlux"]
SurfaceFluxMom = parsed_args["SurfaceFluxMom"]
SurfaceScheme = parsed_args["SurfaceScheme"]
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
RefineLevel = parsed_args["RefineLevel"]
nLon = parsed_args["nLon"]
nLat = parsed_args["nLat"]
LatB = parsed_args["LatB"]
H = parsed_args["H"]
Stretch = parsed_args["Stretch"]
StretchType = parsed_args["StretchType"]
TopoS = parsed_args["TopoS"]
GridType = parsed_args["GridType"]
AdaptGridType = parsed_args["AdaptGridType"]
RadEarth = parsed_args["RadEarth"]
ScaleFactor = parsed_args["ScaleFactor"]
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
OrdPrint = parsed_args["OrdPrint"]
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



# Physical parameters
Phys = DyCore.PhysParameters{FTB}()

#ModelParameters
Model = DyCore.ModelStruct{FTB}()

# Initial conditions
Model.NumV=NumV
Model.NumTr=NumTr
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
#Model.PertTh = PertTh
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
Model.TkePos  = TkePos
Model.HorLimit = HorLimit
Model.Upwind = Upwind
Model.Damping = Damping
Model.Geos = Geos
Model.StrideDamp = StrideDamp
Model.Relax = Relax
Model.Coriolis = Coriolis
Model.CoriolisType = CoriolisType
Model.Buoyancy = Buoyancy
Model.Turbulence = Turbulence
Model.VerticalDiffusion = VerticalDiffusion
Model.VerticalDiffusionMom = VerticalDiffusionMom
Model.JacVerticalDiffusion = JacVerticalDiffusion
Model.JacVerticalAdvection = JacVerticalAdvection
Model.Source = Source
Model.Forcing = Forcing
Model.Thermo = Thermo
Model.Microphysics = Microphysics
Model.TypeMicrophysics = TypeMicrophysics
Model.RelCloud = RelCloud
Model.Rain = Rain
Model.SurfaceFlux = SurfaceFlux
Model.SurfaceFluxMom = SurfaceFluxMom
Model.Curl = Curl
Model.ModelType = ModelType
Model.State = State
Model.Stretch = Stretch
Model.StretchType = StretchType
Model.HyperVisc = HyperVisc
Model.HyperDCurl = HyperDCurl 
Model.HyperDGrad = HyperDGrad
Model.HyperDRhoDiv = HyperDRhoDiv
Model.HyperDDiv = HyperDDiv
Model.HyperDDivW = HyperDDivW
Model.EDMF = EDMF
Model.NDEDMF = NDEDMF

OrdPolyZ = 1
if RadEarth == 0.0
  RadEarth = Phys.RadEarth
  if ScaleFactor != 0.0
    RadEarth = RadEarth / ScaleFactor  
  end  
end

# Equation
if Equation == "CompressibleShallow"
  Model.Equation = Models.CompressibleShallow()  
elseif Equation == "CompressibleDeep"
  Model.Equation = Models.CompressibleDeep()  
end  

ns=50
Grid, Exchange = Grids.InitGridSphere(backend,FTB,OrdPoly,nz,nPanel,RefineLevel,ns,nLon,nLat,LatB,GridType,Decomp,RadEarth,Model,ParallelCom)


Topography = (TopoS=TopoS,H=H,Rad=RadEarth)
#Topography
if TopoS == "BaroWaveHill"
  TopoProfile = Examples.BaroWaveHill()()
elseif TopoS == "SchaerSphereCircle"
  TopoProfile = Examples.SchaerSphereCircle()(Param,Phys)
else
  TopoProfile = Examples.Flat()()
end

Grid.AdaptGrid = Grids.AdaptGrid(FTB,AdaptGridType,FTB(H))

if ParallelCom.Proc == 1
  @show "InitSphere"
end  
(CG, Metric, Global) = DyCore.InitSphere(backend,FTB,OrdPoly,OrdPolyZ,H,Topography,Model,
  Phys,TopoProfile,Exchange,Grid,ParallelCom)
#(CG, Metric, Exchange, Global) = DyCore.InitSphere(backend,FTB,OrdPoly,OrdPolyZ,nz,nPanel,H,
# GridType,Topography,Decomp,Model,Phys,RadEarth,TopoProfile)

# Initial values
Examples.InitialProfile!(backend,FTB,Model,Problem,Param,Phys)
U = GPU.InitialConditions(backend,FTB,CG,Metric,Phys,Global,Model.InitialProfile,Param)

#Coriolis
if Coriolis
  if CoriolisType == "Shallow"
    CoriolisFun = GPU.CoriolisShallow()(Phys)
    Model.CoriolisFun = CoriolisFun  
  elseif CoriolisType == "Deep"
    CoriolisFun = GPU.CoriolisDeep()(Phys)
    Model.CoriolisFun = CoriolisFun  
  else  
    CoriolisFun = GPU.CoriolisNo()()
    Model.CoriolisFun = CoriolisFun
  end  
else
  CoriolisFun = GPU.CoriolisNo()()
  Model.CoriolisFun = CoriolisFun
end

#Buoyancy
if Buoyancy
  if Equation == "CompressibleShallow"
    GravitationFun = GPU.GravitationShallow()(Phys)
    Model.GravitationFun = GravitationFun
  elseif Equation == "CompressibleDeep"
    GravitationFun = GPU.GravitationDeep()(Phys)
    Model.GravitationFun = GravitationFun
  else
    GravitationFun = GPU.GravitationNo()()
    Model.GravitationFun = GravitationFun
  end
else
  GravitationFun = GPU.GravitationNo()()
  Model.GravitationFun = GravitationFun
end

# Pressure
if State == "Dry"
  Pressure, dPresdRhoTh = Models.Dry()(Phys)
  Model.Pressure = Pressure
  Model.dPresdRhoTh = dPresdRhoTh
elseif State == "DryInternalEnergy"
  Pressure, dPresdRhoTh = Models.DryInternalEnergy()(Phys)
  Model.Pressure = Pressure
  Model.dPresdRhoTh = dPresdRhoTh
elseif State == "DryTotalEnergy"
  Pressure, dPresdRhoTh = Models.DryTotalEnergy()(Phys)
  Model.Pressure = Pressure
  Model.dPresdRhoTh = dPresdRhoTh
elseif State == "Moist"
  Pressure, dPresdRhoTh = Models.Moist()(Phys,Model.RhoPos,Model.ThPos,
    Model.RhoVPos,Model.RhoCPos)
  Model.Pressure = Pressure
  Model.dPresdRhoTh = dPresdRhoTh
elseif State == "MoistInternalEnergy"
  Pressure, dPresdRhoTh = Models.MoistInternalEnergy()(Phys,Model.RhoPos,Model.ThPos,
    Model.RhoTPos)
  Model.Pressure = Pressure
  Model.dPresdRhoTh = dPresdRhoTh
elseif State == "IceInternalEnergy"
  Pressure, dPresdRhoTh = Models.IceInternalEnergy()(Phys,Model.RhoPos,Model.ThPos,
    Model.RhoTPos)
  Model.Pressure = Pressure
  Model.dPresdRhoTh = dPresdRhoTh
elseif State == "ShallowWater"
  Pressure = Models.ShallowWaterState()(Phys)
  Model.Pressure = Pressure
end

#Microphysics
if Microphysics
  if TypeMicrophysics == "SimpleMicrophysics"
    MicrophysicsSource  = Models.SimpleMicrophysics()(Phys,Model.RhoPos,Model.ThPos,
      Model.RhoVPos+NumV,Model.RhoCPos+NumV,Model.RelCloud,Model.Rain)
    Model.MicrophysicsSource = MicrophysicsSource
  elseif TypeMicrophysics == "OneMomentMicrophysicsMoistEquil"
    T_TPos = 2
    T_RhoVPos = 5
    T_RhoCPos = 6
    MicrophysicsSource, SedimentationSource = Models.OneMomentMicrophysicsMoistEquil()(Phys,
      Model.RhoPos,Model.ThPos,Model.RhoTPos,Model.RhoRPos,T_TPos,T_RhoVPos,T_RhoCPos,Grid.nz)
    Model.MicrophysicsSource = MicrophysicsSource
    @show MicrophysicsSource
    Model.SedimentationSource = SedimentationSource
  else
    @show "False Type Microphysics"
  end
end

# Damping
if Damping
  Damp = GPU.DampingW()(FTB(H),FTB(StrideDamp),FTB(Relax),Model.wPos)
  Model.Damp = Damp
end

# Surface flux
Global.SurfaceData = Surfaces.SurfaceData{FTB}(backend,CG.NumG)
SD = Global.SurfaceData
Global.LandUseData = Surfaces.LandUseData{FTB}(backend,CG.NumG)
@. Global.LandUseData.z0M = 0.01
@. Global.LandUseData.z0H = 0.01
@. Global.LandUseData.LandClass = 5
if Model.SurfaceFlux || Model.VerticalDiffusion || Model.SurfaceFluxMom || Model.VerticalDiffusionMom
  if SurfaceScheme == ""
    if Problem == "HeldSuarezMoistSphere" || Problem == "HeldSuarezMoistSphereOro" ||
      Problem == "HeldSuarezDrySphere" || Problem == "HeldSuarezDrySphereOro"   
      SurfaceValues, SurfaceFluxValues = Surfaces.HeldSuarezMoistSurface()(Phys,Param,Model.uPos,
        Model.vPos,Model.wPos)
      Model.SurfaceValues = SurfaceValues
      Model.SurfaceFluxValues = SurfaceFluxValues
    elseif Problem == "HeldSuarezDrySphere" || Problem == "HeldSuarezDrySphereOro"   
      SurfaceValues, SurfaceFluxValues = Surfaces.HeldSuarezDrySurface()(Phys,Param,Model.uPos,Model.vPos,Model.wPos)
      Model.SurfaceValues = SurfaceValues
      Model.SurfaceFluxValues = SurfaceFluxValues
    elseif Problem == "FriersonSphere" 
      SurfaceValues, SurfaceFluxValues = Surfaces.FriersonSurface()(Phys,Param,Model.RhoPos,Model.uPos,
        Model.vPos,Model.wPos,Model.ThPos)
      Model.SurfaceValues = SurfaceValues
      Model.SurfaceFluxValues = SurfaceFluxValues
    end  
  elseif SurfaceScheme == "MOST"
    SurfaceValues, SurfaceFluxValues = Surfaces.MOSurface()(Surfaces.Businger(),Phys,Model.RhoPos,Model.uPos,
      Model.vPos,Model.wPos,Model.ThPos)
    Model.SurfaceValues = SurfaceValues
    Model.SurfaceFluxValues = SurfaceFluxValues
    if Problem == "HeldSuarezMoistSphere" || Problem == "HeldSuarezMoistSphereOro" 
      SurfaceValues, SurfaceFluxValues = Surfaces.HeldSuarezMoistSurface()(Phys,Param,Model.uPos,Model.vPos,Model.wPos)
      Model.SurfaceValues = SurfaceValues
    end   
  end  
end
if Model.SurfaceFlux
  Model.SurfaceFluxRhs = Surfaces.SurfaceFlux(Phys,Param,Model.ThPos,Model.RhoPos,Model.RhoVPos)  
end  

#Vertical Diffusion
if Model.VerticalDiffusion || Model.VerticalDiffusionMom
  if Model.Turbulence
    Model.Eddy = Models.TkeKoefficient()(Param,Phys,TkePos,Model.RhoPos)
  elseif Problem == "FriersonSphere"   
    Model.Eddy = Surfaces.FriersonKoefficient()(Param,Phys,Model.RhoPos,Model.uPos,Model.vPos,
      Model.ThPos)
  else  
    Model.Eddy = Surfaces.SimpleKoefficient()(Param,Phys)
  end
end  

#Turbulence
if Model.Turbulence
  Model.TurbulenceSource = Examples.TKEModel()(Param,Phys,Model.RhoPos,Model.uPos,
    Model.vPos,Model.ThPos,Model.TkePos)
end  

# HyperViscosity
  MinGridLength, MaxGridLength = Grids.GridLength(Grid)
  if ParallelCom.Proc == 1
    @show MinGridLength
    @show MaxGridLength
  end  

# Output
Global.Output.vtkFileName = string(Problem*"_")
Global.Output.vtk = 0
Global.Output.Flat = Flat
Global.Output.nPanel = nPanel
Global.Output.RadPrint = H
Global.Output.H = H
if ModelType == "VectorInvariant" || ModelType == "Advection"
  if State == "Dry" || State == "DryTotalEnergy"
    Global.Output.cNames = [
      "Rho",
      "u",
      "v",
      "wB",
      "Th",
      "Pres",
      ]
    if TkePos > 0
      push!(Global.Output.cNames,"Tke")
    end  
    if VerticalDiffusion
      push!(Global.Output.cNames,"DiffKoeff")
    end  
  elseif State == "Moist"  
    Global.Output.cNames = [
      "Rho",
      "u",
      "v",
      "wB",
      "Th",
#     "Vort",
      "Pres",
      "Tr1",
      "Tr2",
      ]
    Model.TrPos = NumV + 1
    if TkePos > 0
      push!(Global.Output.cNames,"Tke")
      Model.TrPos += 1
    end  
    if VerticalDiffusion
      push!(Global.Output.cNames,"DiffKoeff")
    end  

  elseif State == "ShallowWater"  
    Global.Output.cNames = [
      "Rho",
      "u",
      "v",
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
if OrdPrint == 0
  Global.Output.OrdPrint = CG.OrdPoly
else  
  Global.Output.OrdPrint = OrdPrint
end  

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
if ModelType == "VectorInvariant" || ModelType == "Advection"
  DiscType = Val(:VectorInvariant)  
elseif ModelType == "Conservative"
  DiscType = Val(:Conservative)  
end  



Global.ParallelCom.NumberThreadGPU = NumberThreadGPU   
nT = max(7 + NumTr, NumV + NumTr) + Model.NDEDMF*(4 + NumTr)
Parallels.InitExchangeData3D(backend,FTB,nz,nT,Exchange)
if ParallelCom.Proc == 1
  @show "vor Timestepper"
end  
Integration.TimeStepper!(U,GPU.FcnGPU!,GPU.FcnPrepareGPU!,DyCore.JacSchurGPU!,
  Grids.TransSphereX,CG,Metric,Phys,Exchange,Global,Param,Model.Equation)
