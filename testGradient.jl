using CGDycore
using MPI
using Base

# Model
parsed_args = CGDycore.parse_commandline()
Problem = parsed_args["Problem"]
ProfRho = parsed_args["ProfRho"]
ProfTheta = parsed_args["ProfTheta"]
ProfVel = parsed_args["ProfVel"]
ProfpBGrd = parsed_args["ProfpBGrd"]
ProfRhoBGrd = parsed_args["ProfRhoBGrd"]
ProfTest = parsed_args["ProfTest"]
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
PrintTime = parsed_args["PrintTime"]

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
  Model.ProfpBGrd = ProfpBGrd
  Model.ProfRhoBGrd = ProfRhoBGrd
  Model.RefProfile = RefProfile
  Model.ProfTest = ProfTest
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
  Model.ModelType = ModelType



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
            P4=P4,
            )

  Grid=CGDycore.Grid(nz,Topography)
  Grid=CGDycore.CartGrid(nx,ny,Lx,Ly,x0,y0,CGDycore.OrientFaceCart,Boundary,Grid)

  CellToProc=zeros(0)
  Proc = 0
  ProcNumber = 0
  if stretch
    sigma = 1.0
    lambda = 3.16
    CGDycore.AddStretchICONVerticalGrid!(Grid,nz,H,sigma,lambda)
  else
    CGDycore.AddVerticalGrid!(Grid,nz,H)
  end
  Exchange = CGDycore.InitExchangeCG()
  Output=CGDycore.Output(Topography)
  Global = CGDycore.Global(Grid,Model,TimeStepper,ParallelCom,Phys,Output,Exchange,OrdPoly+1,nz,NumV,NumTr,())
  Global.Metric=CGDycore.Metric(OrdPoly+1,OrdPolyZ+1,Grid.NumFaces,nz)
  if TopoS == "EarthOrography"
    (CG,Global)=CGDycore.DiscretizationCG(OrdPoly,OrdPolyZ,CGDycore.JacobiDG3Neu,Global,zS)
  else
   (CG,Global)=CGDycore.DiscretizationCG(OrdPoly,OrdPolyZ,CGDycore.JacobiDG3Neu,Global)
  end

  Model.HyperVisc = HyperVisc
  Model.HyperDCurl = HyperDCurl # =7.e15
  Model.HyperDGrad = HyperDGrad # =7.e15
  Model.HyperDDiv = HyperDDiv # =7.e15


  RhoPos = Model.RhoPos
  uPos = Model.uPos
  vPos = Model.vPos
  wPos = Model.wPos
  ThPos = Model.ThPos
  cC = zeros(Float64,nz,CG.NumG,2)
  Grad = zeros(Float64,nz,CG.NumG,NumV)
  GradEx = zeros(Float64,nz,CG.NumG,NumV)
  cC[:,:,1] = CGDycore.Project(CGDycore.fScalar,0.0,CG,Global,Param)
  GradEx[:,:,uPos],GradEx[:,:,vPos] = CGDycore.ProjectVec(CGDycore.fGrad12,0.0,CG,Global,Param)
  GradEx[:,:,wPos] = CGDycore.ProjectW(CGDycore.fGrad3,0.0,CG,Global,Param)
  @views @. GradEx[:,:,RhoPos] = 1.0
  @views @. GradEx[:,:,ThPos] = cC[:,:,1]



  Global.Cache=CGDycore.CacheCreate(CG.OrdPoly+1,Global.Grid.NumFaces,Global.Grid.NumGhostFaces,CG.NumG,nz,
    NumV,NumTr)

  @views CGDycore.Fcn!(Grad,cC[:,:,1],CG,Global,Param,Val(:TestGrad))
# for iG = 1 : CG.NumG
#   @show GradEx[1:4,iG,uPos]
#   @show Grad[1:4,iG,uPos]
#   @show GradEx[1:4,iG,wPos]
#   @show Grad[1:4,iG,wPos]
# end
  @show sum(abs.(GradEx[:,:,uPos:wPos]))
  @show sum(abs.(GradEx[:,:,uPos:wPos]-Grad[:,:,uPos:wPos]))



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
  ]
  Output.OrdPrint=CG.OrdPoly

  @views @. Grad[:,:,RhoPos] = 1.0
  @views @. Grad[:,:,ThPos] = cC[:,:,1]
  Global.vtkCache = CGDycore.vtkInit3D(Output.OrdPrint,CGDycore.TransCartX,CG,Global)
  CGDycore.unstructured_vtkSphere(Grad,CGDycore.TransCartX,CG,Global,Proc,ProcNumber)
  CGDycore.unstructured_vtkSphere(GradEx,CGDycore.TransCartX,CG,Global,Proc,ProcNumber)


  cC[:,:,2] = CGDycore.Project(CGDycore.fScalar,0.0,CG,Global,Param)
  Global.Model.ProfTest = "TestGrad2"
  GradEx[:,:,uPos],GradEx[:,:,vPos] = CGDycore.ProjectVec(CGDycore.fGrad12,0.0,CG,Global,Param)
  GradEx[:,:,wPos] = CGDycore.ProjectW(CGDycore.fGrad3,0.0,CG,Global,Param)
  CGDycore.Fcn!(Grad,cC,CG,Global,Param,Val(:TestFunCGrad))
# for iG = 1 : CG.NumG
#   @show GradEx[1:4,iG,uPos]
#   @show Grad[1:4,iG,uPos]
#   @show GradEx[1:4,iG,wPos]
#   @show Grad[1:4,iG,wPos]
# end
  @show sum(abs.(GradEx[:,:,uPos:wPos]))
  @show sum(abs.(GradEx[:,:,uPos:wPos]-Grad[:,:,uPos:wPos]))
  CGDycore.unstructured_vtkSphere(Grad,CGDycore.TransCartX,CG,Global,Proc,ProcNumber)
  CGDycore.unstructured_vtkSphere(GradEx,CGDycore.TransCartX,CG,Global,Proc,ProcNumber)
