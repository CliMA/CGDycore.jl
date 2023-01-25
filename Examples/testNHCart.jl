using CGDycore
using MPI
using Base

# Model
parsed_args = CGDycore.parse_commandline()
Problem = parsed_args["Problem"]
ProfRho = parsed_args["ProfRho"]
ProfTheta = parsed_args["ProfTheta"]
ProfVel = parsed_args["ProfVel"]
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
#Orography
TopoS = parsed_args["TopoS"]
P1 = parsed_args["P1"]
P2 = parsed_args["P2"]
P3 = parsed_args["P3"]

# Parallel
Decomp = parsed_args["Decomp"]
SimDays = parsed_args["SimDays"]
SimHours = parsed_args["SimHours"]
SimMinutes = parsed_args["SimMinutes"]
SimSeconds = parsed_args["SimSeconds"]
dtau = parsed_args["dtau"]
IntMethod = parsed_args["IntMethod"]
Table = parsed_args["Table"]
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
H = parsed_args["H"]
stretch = parsed_args["stretch"]
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
  Model.ProfpBGrd="Isothermal"
  Model.ProfRhoBGrd="Isothermal"
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
  Model.Source = Source
  Model.Microphysics = Microphysics
  Model.Source = Source
  Model.SurfaceFlux = SurfaceFlux
  Model.Thermo = "" #"InternalEnergy" #"" #"TotalEnergy"



# Grid
TimeStepper=CGDycore.TimeStepper()


Boundary = CGDycore.Boundary()
Boundary.WE = BoundaryWE
Boundary.SN = BoundarySN
Boundary.BT = BoundaryBT
Topography=(TopoS=TopoS,
            H=H,
            P1=P1,
            P2=P2,
            P3=P3,
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
  Exchange = CGDycore.InitExchangeCG(SubGrid,OrdPoly,CellToProc,Proc,ProcNumber,Parallel)
  Output=CGDycore.Output(Topography)
  Global = CGDycore.Global(SubGrid,Model,TimeStepper,ParallelCom,Phys,Output,Exchange,OrdPoly+1,nz,NumV,NumTr,())
  Global.Metric=CGDycore.Metric(OrdPoly+1,OrdPolyZ+1,SubGrid.NumFaces,nz)
  (CG,Global)=CGDycore.DiscretizationCG(OrdPoly,OrdPolyZ,CGDycore.JacobiDG3,Global)
  # Output partition
  nzTemp = Global.Grid.nz
  Global.Grid.nz = 1
  vtkCachePart = CGDycore.vtkInit3D(1,CGDycore.TransCartX,CG,Global)
  CGDycore.unstructured_vtkPartition(vtkCachePart, Global.Grid.NumFaces, Proc, ProcNumber)
  Global.Grid.nz = nzTemp

  if TopoS == "EarthOrography"
    SubGrid = CGDycore.StencilFace(SubGrid)
    zS = CGDycore.Orography(CG,Global)
    Output.RadPrint = H
    Output.Flat=false
    nzTemp = Global.Grid.nz
    Global.Grid.nz = 1
    vtkCacheOrography = CGDycore.vtkInit2D(CG.OrdPoly,CGDycore.TransCartX,CG,Global)
    CGDycore.unstructured_vtkOrography(zS,vtkCacheOrography, Global.Grid.NumFaces, CG,  Proc, ProcNumber)
    Global.Grid.nz = nzTemp
  end
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
if TopoS == "EarthOrography"
  (CG,Global)=CGDycore.DiscretizationCG(OrdPoly,OrdPolyZ,CGDycore.JacobiDG3,Global,zS)
else
  (CG,Global)=CGDycore.DiscretizationCG(OrdPoly,OrdPolyZ,CGDycore.JacobiDG3,Global)
end

Model.HyperVisc = HyperVisc
Model.HyperDCurl = HyperDCurl # =7.e15
Model.HyperDGrad = HyperDGrad # =7.e15
Model.HyperDDiv = HyperDDiv # =7.e15


  U = CGDycore.InitialConditions(CG,Global,Param)

# Output partition  
  nzTemp = Global.Grid.nz
  Global.Grid.nz = 1
  vtkCachePart = CGDycore.vtkInit2D(1,CGDycore.TransCartX,CG,Global)
  Global.Grid.nz = nzTemp
  CGDycore.unstructured_vtkPartition(vtkCachePart, Global.Grid.NumFaces, Proc, ProcNumber)

# Output
  Output.vtkFileName=string(Problem*"_")
  Output.vtk=0
  Output.Flat=true
  Output.H=H
  Output.cNames = [
    "Rho",
    "u",
    "v",
    "w",
    "Th",
    "Pres",
]
  Output.PrintDays = PrintDays
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
  CGDycore.TimeStepper!(U,CGDycore.TransCartX,CG,Global,Param)
