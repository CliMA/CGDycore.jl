using CGDycore
using MPI
using Base

# Model
parsed_args = CGDycore.parse_commandline()
Problem = parsed_args["Problem"]
ProfRho = parsed_args["ProfRho"]
ProfTheta = parsed_args["ProfTheta"]
ProfVel = parsed_args["ProfVel"]
ProfVelW = parsed_args["ProfVelW"]
ProfTr = parsed_args["ProfTr"]
HorLimit = parsed_args["HorLimit"]
Upwind = parsed_args["Upwind"]
Damping = parsed_args["Damping"]
Relax = parsed_args["Relax"]
StrideDamp = parsed_args["StrideDamp"]
NumV = parsed_args["NumV"]
NumTr = parsed_args["NumTr"]
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
TopoS = parsed_args["TopoS"]
GridType = parsed_args["GridType"]
Coriolis = parsed_args["Coriolis"]
CoriolisType = parsed_args["CoriolisType"]
Microphysics = parsed_args["Microphysics"]
Source = parsed_args["Source"]
VerticalDiffusion = parsed_args["VerticalDiffusion"]
SurfaceFlux = parsed_args["SurfaceFlux"]
# Grid
RadEarth = parsed_args["RadEarth"]
nz = parsed_args["nz"]
nPanel = parsed_args["nPanel"]
H = parsed_args["H"]
Stretch = parsed_args["Stretch"]
OrdPoly = parsed_args["OrdPoly"]
# Viscosity
HyperVisc = parsed_args["HyperVisc"]
HyperDCurl = parsed_args["HyperDCurl"]
HyperDGrad = parsed_args["HyperDGrad"]
HyperDDiv = parsed_args["HyperDDiv"]
# Output
vtkFileName = parsed_args["vtkFileName"]
PrintDays = parsed_args["PrintDays"]
PrintHours = parsed_args["PrintHours"]
PrintMinutes = parsed_args["PrintMinutes"]
PrintSeconds = parsed_args["PrintSeconds"]
PrintTime = parsed_args["PrintTime"]
Flat = parsed_args["Flat"]

Param = CGDycore.Parameters(Problem)

MPI.Init()

OrdPolyZ=1

# Physical parameters
Phys=CGDycore.PhysParameters()

#ModelParameters
Model = CGDycore.Model()
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
if ProfTr == ""
  Model.ProfTr = Problem
else
  Model.ProfTr = ProfTr  
end  
Model.RhoPos=1
Model.uPos=0
Model.vPos=0
Model.wPos=0
Model.HorLimit = HorLimit
Model.Upwind = Upwind
Model.HyperVisc = HyperVisc
Model.HyperDCurl = HyperDCurl # =7.e15
Model.HyperDGrad = HyperDGrad # =7.e15
Model.HyperDDiv = HyperDDiv # =7.e15

if RadEarth == 0.0
  RadEarth = Phys.RadEarth
end  


# Grid
Topography=(TopoS=TopoS,H=H,Rad=RadEarth)

(CG,Global) = CGDycore.InitSphere(OrdPoly,OrdPolyZ,nz,nPanel,H,GridType,Topography,Decomp,Model,Phys)

U = CGDycore.InitialConditionsAdvection(CG,Global,Param)

# Output
  Global.Output.vtkFileName=string(vtkFileName*"_")
  Global.Output.vtk=0
  Global.Output.Flat=Flat
  Global.Output.nPanel=nPanel
  Global.Output.RadPrint=H
  Global.Output.H=H
  Global.Output.cNames = [
    "Rho",
    "Tr1",
]
  Global.Output.PrintDays = PrintDays
  Global.Output.PrintHours = PrintHours
  Global.Output.PrintSeconds = PrintSeconds
  Global.Output.PrintTime = PrintTime
  Global.Output.PrintStartTime = 0
  Global.Output.OrdPrint = CG.OrdPoly
  Global.vtkCache = CGDycore.vtkStruct(Global.Output.OrdPrint,CGDycore.TransSphereX,CG,Global)

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

  CGDycore.TimeStepperAdvectionConv!(U,CGDycore.TransSphereX,CG,Global,Param)
