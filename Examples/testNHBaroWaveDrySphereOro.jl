using CGDycore
using MPI
using Base

parsed_args = CGDycore.parse_commandline()
Problem = parsed_args["Problem"]
Decomp = parsed_args["Decomp"]
SimDays = parsed_args["SimDays"]
dtau = parsed_args["dtau"]
IntMethod = parsed_args["IntMethod"]
Table = parsed_args["Table"]
TopoS = parsed_args["TopoS"]
stretch = parsed_args["stretch"]
GridType = parsed_args["GridType"]

Base.@kwdef struct ParamStruct
  T0E=310.0
  T0P=240.0
  B=2.0
  K=3.0
  LapseRate=0.005
  U0=-0.5
  PertR=1.0/6.0
  Up=1.0
  PertExpR=0.1
  PertLon=pi/9.0
  PertLat=2.0*pi/9.0
  PertZ=15000.0
  NBr=1.e-2
  DeltaT=1
  ExpDist=5
  T0=300
  TEq=300
  T_init = 315
  lapse_rate = -0.008
  Deep=false
  pert = 0.1
  uMax = 1.0
  vMax = 0.0
  DeltaT_y=0
  DeltaTh_z=-5
  T_equator=315
  T_min=200
  sigma_b=7/10
  z_D=20.0e3
  #      Moist
  q_0 = 0.018                # Maximum specific humidity (default: 0.018)
  q_t = 1.0e-12
end  

MPI.Init()
comm = MPI.COMM_WORLD
Proc = MPI.Comm_rank(comm) + 1
ProcNumber = MPI.Comm_size(comm)
ParallelCom = CGDycore.ParallelCom()
ParallelCom.Proc = Proc
ParallelCom.ProcNumber  = ProcNumber
print("$Proc: \n")
print("$ProcNumber: \n")

OrdPoly = 4
nz = 45
OrdPolyZ=1
nPanel = 16
NF = 6 * nPanel * nPanel
NumV = 5
NumTr = 0
Parallel = true

# Physical parameters
Phys=CGDycore.PhysParameters()

#ModelParameters
Param = ParamStruct()
Model = CGDycore.Model()
# Initial conditions
  Model.Equation="Compressible"
  Model.NumV=NumV
  Model.NumTr=NumTr
  Model.Problem="BaroWaveSphere"
  Model.ProfRho="BaroWaveSphere"
  Model.ProfTheta="BaroWaveSphere"
  Model.ProfVel="BaroWaveSphere"
  Model.ProfpBGrd="Isothermal"
  Model.ProfRhoBGrd="Isothermal"
  Model.RhoPos=1
  Model.uPos=2
  Model.vPos=3
  Model.wPos=4
  Model.ThPos=5
  Model.HorLimit = false
  Model.Source = false
  Model.Upwind = true
  Model.Damping = true
  Model.StrideDamp=10000.0
  Model.Relax = 1.0/100.0
  Model.Coriolis=true
  Model.CoriolisType="Sphere"
  Model.VerticalDiffusion = false
  Model.Thermo = "" #"InternalEnergy" #"" #"TotalEnergy"

# Grid
H = 30000.0
Topography=(TopoS=TopoS,H=H,Rad=Phys.RadEarth)
TimeStepper=CGDycore.TimeStepper()

Grid=CGDycore.Grid(nz,Topography)
if GridType == "Helix"
  Grid=CGDycore.InputGridH("Grid/mesh_H12_no_pp.nc",
  CGDycore.OrientFaceSphere,Phys.RadEarth,Grid)
elseif GridType == "SQuadGen"
  Grid=CGDycore.InputGrid("Grid/baroclinic_wave_2deg_x4.g",
  CGDycore.OrientFaceSphere,Phys.RadEarth,Grid)
elseif GridType == "CubedSphere"
  Grid=CGDycore.CubedGrid(nPanel,CGDycore.OrientFaceSphere,Phys.RadEarth,Grid)
end
if Parallel
  if Decomp == "Hilbert"
    CGDycore.HilbertFaceSphere!(Grid)
    CellToProc = CGDycore.Decompose(Grid,ProcNumber)
  elseif Decomp == "EqualArea"
    CellToProc = CGDycore.DecomposeEqualArea(Grid,ProcNumber)
  else
    CellToProc = ones(Int,Grid.NumFaces)
    println(" False Decomp method ")
  end
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
  (CG,Global)=CGDycore.DiscretizationCG(OrdPoly,OrdPolyZ,CGDycore.JacobiSphere3,Global)
  # Output partition
  nzTemp = Global.Grid.nz
  Global.Grid.nz = 1
  vtkCachePart = CGDycore.vtkInit3D(1,CGDycore.TransSphereX,CG,Global)
  CGDycore.unstructured_vtkPartition(vtkCachePart, Global.Grid.NumFaces, Proc, ProcNumber)
  Global.Grid.nz = nzTemp

  if TopoS == "EarthOrography"
    SubGrid = CGDycore.StencilFace(SubGrid)
    zS = CGDycore.Orography(CG,Global)
    Output.RadPrint = H
    Output.Flat=false
    nzTemp = Global.Grid.nz
    Global.Grid.nz = 1
    vtkCacheOrography = CGDycore.vtkInit2D(CG.OrdPoly,CGDycore.TransSphereX,CG,Global)
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
  (CG,Global)=CGDycore.DiscretizationCG(OrdPoly,OrdPolyZ,CGDycore.JacobiSphere3,Global,zS)
else
  (CG,Global)=CGDycore.DiscretizationCG(OrdPoly,OrdPolyZ,CGDycore.JacobiSphere3,Global)
end

Model.HyperVisc=true
Model.HyperDCurl=7.e15
Model.HyperDGrad=7.e15
Model.HyperDDiv=7.e15


  U = CGDycore.InitialConditions(CG,Global,Param)

# Output partition  
  nzTemp = Global.Grid.nz
  Global.Grid.nz = 1
  vtkCachePart = CGDycore.vtkInit2D(1,CGDycore.TransSphereX,CG,Global)
  Global.Grid.nz = nzTemp
  CGDycore.unstructured_vtkPartition(vtkCachePart, Global.Grid.NumFaces, Proc, ProcNumber)

# Output
  Output.vtkFileName=string(Problem*"_")
  Output.vtk=0
  Output.Flat=true
  Output.nPanel=nPanel
  Output.RadPrint=H
  Output.H=H
  Output.cNames = [
    "Rho",
    "u",
    "v",
    "w",
    "Th",
    "Pres",
]
  Output.PrintDay = 0.5
  Output.PrintStartDay = 0
  Output.OrdPrint=CG.OrdPoly
  Global.vtkCache = CGDycore.vtkInit3D(Output.OrdPrint,CGDycore.TransSphereX,CG,Global)

  # TimeStepper
  time=[0.0]
  TimeStepper.IntMethod = IntMethod
  TimeStepper.Table = Table
  TimeStepper.dtau = dtau
  TimeStepper.SimDays = SimDays
  CGDycore.TimeStepper!(U,CGDycore.TransSphereX,CG,Global,Param)
