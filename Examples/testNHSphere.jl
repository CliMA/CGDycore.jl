using CGDycore
using MPI
using Base

# Model
parsed_args = CGDycore.parse_commandline()
Problem = parsed_args["Problem"]
ProfRho = parsed_args["ProfRho"]
ProfTheta = parsed_args["ProfTheta"]
PertTh = parsed_args["PertTh"]
ProfVel = parsed_args["ProfVel"]
ProfVelGeo = parsed_args["ProfVelGeo"]
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
RelCloud = parsed_args["RelCloud"]
Rain = parsed_args["Rain"]
Source = parsed_args["Source"]
VerticalDiffusion = parsed_args["VerticalDiffusion"]
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
# CG Element
OrdPoly = parsed_args["OrdPoly"]
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
PrintStartTime = parsed_args["PrintStartTime"]
Flat = parsed_args["Flat"]

Param = CGDycore.Parameters(Problem)

MPI.Init()

# Physical parameters
Phys=CGDycore.PhysParameters()

#ModelParameters
Model = CGDycore.Model()

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
if Model.Equation == "CompressibleMoist"
  Model.RhoVPos = 1
  Model.RhoCPos = 2
end  
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
Model.Source = Source
Model.Microphysics = Microphysics
Model.RelCloud = RelCloud
Model.Rain = Rain
Model.Source = Source
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
Model.HyperDDiv = HyperDDiv

OrdPolyZ = 1

Topography = (TopoS=TopoS,H=H,Rad=Phys.RadEarth)

(CG,Global) = CGDycore.InitSphere(OrdPoly,OrdPolyZ,nz,nPanel,H,GridType,Topography,Decomp,Model,Phys)



# Initial values
U = CGDycore.InitialConditions(CG,Global,Param)

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
      "Vort",
      "Pres",
      ]
  elseif Model.Equation == "CompressibleMoist"  
    Global.Output.cNames = [
      "Rho",
      "u",
      "v",
      "wB",
      "Th",
      "Vort",
      "Pres",
      "Tr1",
      "Tr2",
      ]
  end  
elseif ModelType == "Conservative"
  Global.Output.cNames = [
    "Rho",
    "Rhou",
    "Rhov",
    "w",
    "Th",
    "Vort",
    "Pres",
    ]
end

Global.Output.PrintDays = PrintDays
Global.Output.PrintHours = PrintHours
Global.Output.PrintMinutes = PrintMinutes
Global.Output.PrintSeconds = PrintSeconds
Global.Output.StartAverageDays = StartAverageDays
Global.Output.PrintStartTime = PrintStartTime
Global.Output.OrdPrint = CG.OrdPoly

Global.vtkCache = CGDycore.vtkStruct(Global.Output.OrdPrint,CGDycore.TransSphereX,CG,Global)

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
nT = max(9 + NumTr, NumV + NumTr)
CGDycore.InitExchangeData3D(nz,nT,Global.Exchange)
CGDycore.TimeStepper!(U,CGDycore.Fcn!,CGDycore.TransSphereX,CG,Global,Param,DiscType)
