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
nx = parsed_args["nx"]
ny = parsed_args["ny"]
nz = parsed_args["nz"]
nPanel = parsed_args["nPanel"]
H = parsed_args["H"]
stretch = parsed_args["stretch"]
OrdPoly = parsed_args["OrdPoly"]
Lx = parsed_args["Lx"]
Ly = parsed_args["Ly"]
x0 = parsed_args["x0"]
y0 = parsed_args["y0"]
BoundaryWE = parsed_args["BoundaryWE"]
BoundarySN = parsed_args["BoundarySN"]
BoundaryBT = parsed_args["BoundaryBT"]
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
Flat = parsed_args["Flat"]
#Orography
TopoS = parsed_args["TopoS"]


Param = CGDycore.Parameters(Problem)

MPI.Init()
comm = MPI.COMM_WORLD
Proc = MPI.Comm_rank(comm) + 1
ProcNumber = MPI.Comm_size(comm)
ParallelCom = CGDycore.ParallelCom()
ParallelCom.Proc = Proc
ParallelCom.ProcNumber  = ProcNumber
print("$Proc: \n")
print("$ProcNumber: \n")

OrdPolyZ=1
Parallel = true

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


# Grid
TimeStepper=CGDycore.TimeStepper()


Boundary = CGDycore.Boundary()
Boundary.WE = BoundaryWE
Boundary.SN = BoundarySN
Boundary.BT = BoundaryBT

Topography=(TopoS=TopoS,
            H=H,
            )

Grid=CGDycore.Grid(nz,Topography)
Grid=CGDycore.CartGrid(nx,ny,Lx,Ly,x0,y0,CGDycore.OrientFaceCart,Boundary,Grid)

if Parallel
  CellToProc = CGDycore.Decompose(Grid,ProcNumber)
  SubGrid = CGDycore.ConstructSubGrid(Grid,CellToProc,Proc)

  if stretch
    sigma = 1.0
    lambda = 3.16
    CGDycore.AddStretchICONVerticalGrid!(SubGrid,nz,H,sigma,lambda)
  else
    CGDycore.AddVerticalGrid!(SubGrid,nz,H)
  end
  Exchange = CGDycore.InitExchangeCG(SubGrid,OrdPoly,CellToProc,Proc,ProcNumber,Parallel,HorLimit)
  Output=CGDycore.Output(Topography)
  Global = CGDycore.Global(SubGrid,Model,TimeStepper,ParallelCom,Phys,Output,Exchange,OrdPoly+1,nz,NumV,NumTr,())
  Global.Metric=CGDycore.Metric(OrdPoly+1,OrdPolyZ+1,SubGrid.NumFaces,nz)
  (CG,Global)=CGDycore.DiscretizationCG(OrdPoly,OrdPolyZ,CGDycore.JacobiDG3Neu,Global)
  # Output partition
  nzTemp = Global.Grid.nz
  Global.Grid.nz = 1
  vtkCachePart = CGDycore.vtkInit3D(1,CGDycore.TransCartX,CG,Global)
  CGDycore.unstructured_vtkPartition(vtkCachePart, Global.Grid.NumFaces, Proc, ProcNumber)
  Global.Grid.nz = nzTemp
  AA = CGDycore.HorLimiter(SubGrid)

else
  CellToProc=zeros(0)
  Proc = 0
  ProcNumber = 0
  sigma = 1.0
  lambda = 3.16
  CGDycore.AddStretchICONVerticalGrid!(Grid,nz,H,sigma,lambda)
  Exchange = CGDycore.InitExchange(Grid,OrdPoly,CellToProc,Proc,ProcNumber,Parallel)
  Output=CGDycore.Output(Topography)
  Global = CGDycore.Global(Grid,Model,Phys,Output,Exchange,OrdPoly+1,nz,NumV,NumTr,())
  Global.Metric=CGDycore.Metric(OrdPoly+1,OrdPolyZ+1,Grid.NumFaces,nz)
end  
  (CG,Global)=CGDycore.DiscretizationCG(OrdPoly,OrdPolyZ,CGDycore.JacobiDG3Neu,Global)

Model.HyperVisc = HyperVisc
Model.HyperDCurl = HyperDCurl # =7.e15
Model.HyperDGrad = HyperDGrad # =7.e15
Model.HyperDDiv = HyperDDiv # =7.e15


U = CGDycore.InitialConditionsAdvection(CG,Global,Param)

# Output partition  
  nzTemp = Global.Grid.nz
  Global.Grid.nz = 1
  vtkCachePart = CGDycore.vtkInit2D(1,CGDycore.TransCartX,CG,Global)
  Global.Grid.nz = nzTemp
  CGDycore.unstructured_vtkPartition(vtkCachePart, Global.Grid.NumFaces, Proc, ProcNumber)

# Output
  Output.vtkFileName=string(vtkFileName*"_")
  Output.vtk=0
  Output.Flat=Flat
  Output.nPanel=nPanel
  Output.RadPrint=H
  Output.H=H
  Output.cNames = [
    "Rho",
    "Tr1",
]
  Output.PrintDays = PrintDays
  Output.PrintHours = PrintHours
  Output.PrintSeconds = PrintSeconds
  Output.PrintStartDays = 0
  Output.OrdPrint=CG.OrdPoly
  Global.vtkCache = CGDycore.vtkInit3D(Output.OrdPrint,CGDycore.TransCartX,CG,Global)

  # TimeStepper
  time=[0.0]
  TimeStepper.IntMethod = IntMethod
  TimeStepper.Table = Table
  TimeStepper.dtau = dtau
  TimeStepper.SimDays = SimDays
  TimeStepper.SimHours = SimHours
  TimeStepper.SimMinutes = SimMinutes
  TimeStepper.SimSeconds = SimSeconds
  CGDycore.TimeStepperAdvection!(U,CGDycore.TransCartX,CG,Global,Param)
