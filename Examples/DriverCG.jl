import CGDycore:
  Parameters, Thermodynamics, Examples, Sources, Parallels, Models, Grids, Surfaces,  Outputs, Integration,  CGSEM, DyCore
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


# Model
parsed_args = Parameters.parse_commandline()
Problem = parsed_args["Problem"]
Discretization = parsed_args["Discretization"]
VelocityForm = parsed_args["VelocityForm"]
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


Param = Examples.Parameters(FTB,Problem)

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
Model.pAuxPos  = 1
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

if VelocityForm == "Spherical"
   VelForm = Examples.VelocityS()
elseif VelocityForm == "Cartesian"
   VelForm = Examples.VelocityC()
end

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


#Topography
if TopoS == "AgnesiHill"
  TopoProfile = Examples.AgnesiHill()()
elseif TopoS == "SchaerHill"
  TopoProfile = Examples.SchaerHill()()
elseif TopoS == "BaroWaveHill"
  TopoProfile = Examples.BaroWaveHill()()
elseif TopoS == "SchaerSphereCircle"
  TopoProfile = Examples.SchaerSphereCircle()(Param,Phys)
elseif TopoS == "GapHillSphere"
  TopoProfile = Examples.GapHillSphere()(Phys,ScaleFactor)
elseif TopoS == "VortexHillSphere"
  TopoProfile = Examples.VortexHillSphere()(Phys,ScaleFactor)
else
  TopoProfile = Examples.Flat()()  
end  

Grid.AdaptGrid = Grids.AdaptGrid(FTB,AdaptGridType,FTB(H))

if GridForm == "Cartesian"
  if ParallelCom.Proc == 1  
    @show "InitCart"
  end  
  Trans = Grids.TransCartX!
  (CG, Metric, Exchange, Global) = DyCore.InitCart(backend,FTB,OrdPoly,OrdPolyZ,OrdPrint,H,
    Topography,Model,Phys,TopoProfile,CellToProc,Grid,ParallelCom)
else  
  if ParallelCom.Proc == 1  
    @show "InitSphere"
  end  
  Trans = Grids.TransSphereX!
  (CG, Metric, Exchange, Global) = DyCore.InitSphere(backend,FTB,OrdPoly,OrdPolyZ,OrdPrint,H,Topography,Model,
    Phys,TopoProfile,CellToProc,Grid,ParallelCom)
end  

# Initial values
Examples.InitialProfile!(backend,FTB,Model,Problem,Param,Phys,VelForm)
U = Examples.InitialConditions(backend,FTB,CG,Metric,Exchange,Phys,Global,Model.InitialProfile,Param)

#Coriolis
if Coriolis
  if CoriolisType == "Shallow"
    CoriolisFun = Sources.CoriolisShallow()(Phys)
    Model.CoriolisFun = CoriolisFun
  elseif CoriolisType == "Deep"
    CoriolisFun = Sources.CoriolisDeep()(Phys)
    Model.CoriolisFun = CoriolisFun
  elseif CoriolisType == "FPlane"
    CoriolisFun = Sources.FPlane()(Param,Phys)
    Model.CoriolisFun = CoriolisFun
  else
    CoriolisFun = Sources.CoriolisNo()()
    Model.CoriolisFun = CoriolisFun
  end
else
  CoriolisFun = Sources.CoriolisNo()()
  Model.CoriolisFun = CoriolisFun
end

#Buoyancy
if Buoyancy
  if Equation == "CompressibleShallow"
    GravitationFun = Sources.GravitationShallow()(Phys)
    Model.GravitationFun = GravitationFun
  elseif Equation == "CompressibleDeep"
    GravitationFun = Sources.GravitationDeep()(Phys)
    Model.GravitationFun = GravitationFun
  else
    GravitationFun = Sources.GravitationNo()()
    Model.GravitationFun = GravitationFun
  end
else
  GravitationFun = Sources.GravitationNo()()
  Model.GravitationFun = GravitationFun
end

# Pressure
if State == "Dry"
  Pressure, dPresdRhoTh, dPresdRho = Models.Dry()(Phys)
  Model.Pressure = Pressure
  Model.dPresdRhoTh = dPresdRhoTh
  Model.dPresdRho = dPresdRho
elseif State == "DryInternalEnergy"
  Pressure, dPresdRhoTh, dPresdRho = Models.DryInternalEnergy()(Phys)
  Model.Pressure = Pressure
  Model.dPresdRhoTh = dPresdRhoTh
  Model.dPresdRho = dPresdRho
elseif State == "DryTotalEnergy"
  Pressure, dPresdRhoTh, dPresdRho = Models.DryTotalEnergy()(Phys)
  Model.Pressure = Pressure
  Model.dPresdRhoTh = dPresdRhoTh
  Model.dPresdRho = dPresdRho
elseif State == "Moist"
  Pressure, dPresdRhoTh, dPresdRho = Models.Moist()(Phys,Model.RhoPos,Model.ThPos,
    Model.RhoVPos,Model.RhoCPos)
  Model.Pressure = Pressure
  Model.dPresdRhoTh = dPresdRhoTh
  Model.dPresdRho = dPresdRho
elseif State == "MoistInternalEnergy"
  if Model.RhoTPos > 0
    Pressure, dPresdRhoTh, dPresdRho = Models.MoistInternalEnergy()(Phys,Model.RhoPos,Model.ThPos,
      Model.RhoTPos)
  elseif Model.RhoVPos > 0  
    Pressure, dPresdRhoTh, dPresdRho = Models.MoistInternalEnergy()(Phys,Model.RhoPos,Model.ThPos,
      Model.RhoVPos,Model.RhoCPos)
  end  
  Model.Pressure = Pressure
  Model.dPresdRhoTh = dPresdRhoTh
  Model.dPresdRho = dPresdRho
elseif State == "MoistTotalEnergy"
  Pressure, dPresdRhoTh, dPresdRho = Models.MoistTotalEnergy()(Phys,Model.RhoPos,Model.ThPos,
    Model.RhoTPos)
  Model.Pressure = Pressure
  Model.dPresdRhoTh = dPresdRhoTh
  Model.dPresdRho = dPresdRho
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
  if TypeMicrophysics == "SimpleMicrophysicsPot"
    MicrophysicsSource  = Models.SimpleMicrophysicsPot()(Phys,Model.RhoPos,Model.ThPos,
      Model.RhoVPos,Model.RhoCPos,Model.RelCloud,Model.Rain)
    Model.MicrophysicsSource = MicrophysicsSource
  elseif TypeMicrophysics == "SimpleMicrophysicsIE"
    MicrophysicsSource  = Models.SimpleMicrophysicsPot()(Phys,Model.RhoPos,Model.ThPos,
      Model.RhoVPos,Model.RhoCPos,Model.RelCloud,Model.Rain)
    Model.MicrophysicsSource = MicrophysicsSource
  elseif TypeMicrophysics == "OneMomentMicrophysicsMoistEquil"
    T_TPos = 2
    T_RhoVPos = 5
    T_RhoCPos = 6
    MicrophysicsSource, SedimentationSource = Models.OneMomentMicrophysicsMoistEquil()(Phys,
      Model.RhoPos,Model.ThPos,Model.RhoTPos,Model.RhoRPos,T_TPos,T_RhoVPos,T_RhoCPos,Grid.nz)
    Model.MicrophysicsSource = MicrophysicsSource
    Model.SedimentationSource = SedimentationSource
  else
    @show "False Type Microphysics"
  end
end

# Surface flux
if Model.SurfaceFlux
# SurfaceValues
  Global.SurfaceData = Surfaces.SurfaceData{FTB}(backend,CG.NumG)
  Global.LandUseData = Surfaces.LandUseData{FTB}(backend,CG.NumG)
  @. Global.LandUseData.LandClass = 2 # Need more Oswald
  if Problem == "HeldSuarezMoistSphere"  
    SurfaceValues = Surfaces.HeldSuarezMoistSurface()(Phys,Param,Model.uPos,Model.vPos,Model.wPos)  
    Model.SurfaceValues = SurfaceValues
  else
    SurfaceValues = Surfaces.DefaultSurface()(Phys,Param,Model.uPos,Model.vPos,Model.wPos)  
    Model.SurfaceValues = SurfaceValues
  end
end  
# SurfaceFlux
if Model.SurfaceFlux || Model.VerticalDiffusion || Model.SurfaceFluxMom || Model.VerticalDiffusionMom
  if SurfaceScheme == ""
    @show "Warning: No surface scheme"  
  elseif SurfaceScheme == "MOST"
    @show "SurfaceScheme MOST"
    SurfaceFluxValues = Surfaces.MOSurfaceFlux()(Surfaces.Businger(),Phys,Model.RhoPos,Model.uPos,
      Model.vPos,Model.wPos,Model.ThPos,Global.LandUseData.z0M,Global.LandUseData.z0H)
    Model.SurfaceFluxValues = SurfaceFluxValues
  end
end
if Model.SurfaceFlux
  if RhoVPos > 0  
    Model.SurfaceFluxRhs = Surfaces.SurfaceFlux(Phys,Param,Model.ThPos,Model.RhoPos,Model.RhoVPos)
  else  
    Model.SurfaceFluxRhs = Surfaces.SurfaceFlux(Phys,Param,Model.ThPos,Model.RhoPos)
  end  
end

#Vertical Diffusion
if Model.VerticalDiffusion || Model.VerticalDiffusionMom
  if Model.Turbulence
    Model.Eddy = Surfaces.TkeKoefficient()(Param,Phys,Model.TkePos,Model.RhoPos)
  else
    Model.Eddy = Surfaces.SimpleKoefficient()(Param,Phys)
  end
end

#Turbulence
if Model.Turbulence
  Model.TurbulenceSource = Models.TKEModel()(Param,Phys,Model.RhoPos,Model.uPos,
    Model.vPos,Model.ThPos,Model.TkePos)
end
# Damping
if Damping
  Damp = Sources.DampingW()(FTB(H),FTB(StrideDamp),FTB(Relax),
    Model.uPos,Model.vPos,Model.wPos,VelForm,Grid.Form)
  Model.Damp = Damp
end  

# Output
Global.Output.vtk=0
Global.Output.Flat = Flat
Global.Output.H = H
if ModelType == "VectorInvariant" || ModelType == "Advection"
  if State == "Dry" || State == "DryInternalEnergy" || State == "DryTotalEnergy"
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
  elseif State == "ShallowWater" 
    Global.Output.cNames = [
      "Rho",
      "u",
      "v",
      ]
  elseif State == "Moist" || State == "MoistInternalEnergy" || State == "MoistTotalEnergy" || State == "IceInternalEnergy"
    Global.Output.cNames = [
      "Rho",
      "u",
      "v",
      "wB",
      "Th",
#     "ThE",
      "Pres",
      ]
    if TkePos > 0
      push!(Global.Output.cNames,"Tke")
    end
    if RhoVPos > 0
      push!(Global.Output.cNames,"qV")
    end
    if RhoTPos > 0
      push!(Global.Output.cNames,"qT")
    end
    if RhoCPos > 0
      push!(Global.Output.cNames,"qC")
    end
    if VerticalDiffusion
      push!(Global.Output.cNames,"DiffKoeff")
    end
  end    
elseif ModelType == "Conservative"
  Global.Output.cNames = [
    "Rho",
    "Rhou",
    "Rhov",
    "w",
    "Th",
    "Vort",
    ]
end

Global.Output.PrintDays = PrintDays
Global.Output.PrintHours = PrintHours
Global.Output.PrintMinutes = PrintMinutes
Global.Output.PrintSeconds = PrintSeconds
Global.Output.PrintTime = PrintTime
Global.Output.PrintStartTime = PrintStartTime
Global.Output.OrdPrintZ = OrdPrintZ
if OrdPrint == 0
  Global.Output.OrdPrint = CG.OrdPoly
else  
  Global.Output.OrdPrint = OrdPrint
end 
Global.Output.vtkFileName = vtkFileName
Global.vtkCache = Outputs.vtkStruct{FTB}(backend,Global.Output.OrdPrint,Global.Output.OrdPrintZ,Trans,CG,Metric,Global)


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
if ModelType == "VectorInvariant" || ModelType == "Advection"
  DiscType = Val(:VectorInvariant)
elseif ModelType == "Conservative"
  DiscType = Val(:Conservative)
end

Global.ParallelCom.NumberThreadGPU = NumberThreadGPU
Global.ParallelCom.NumberThreadTriGPU = NumberThreadTriGPU
if JuliaGPU == "Metal"
  Global.ParallelCom.NumberThreadGPU = NumberThreadGPU / 2
  Global.ParallelCom.NumberThreadTriGPU = NumberThreadTriGPU / 2
end  
nT = max(7 + NumTr, NumV + NumTr)
Parallels.InitExchangeData3D(backend,FTB,1,nz,nT,Exchange)
Integration.TimeStepper!(U,CGSEM.Fcn!,CGSEM.JacGPU!,
  Trans,CG,Metric,Phys,Exchange,Global,Param,Model.Equation)

