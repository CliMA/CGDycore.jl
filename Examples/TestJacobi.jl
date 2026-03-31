import CGDycore:
  Parameters, Thermodynamics, Examples, Sources, Parallels, Models, Grids, Surfaces,  Outputs, Integration, FiniteElements, DGSEM, CGSEM, DyCore, IMEXRosenbrock
using MPI
using Base
using CUDA
#using AMDGPU
#using Metal
using KernelAbstractions
#using StaticArrays
using ArgParse


#=
# Model
parsed_args = Parameters.parse_commandline()
Problem = parsed_args["Problem"]
Discretization = parsed_args["Discretization"]
VelocityForm = parsed_args["VelocityForm"]
FluxDG = parsed_args["FluxDG"]
InterfaceFluxDG = parsed_args["InterfaceFluxDG"]
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
pAuxPos = parsed_args["pAuxPos"]
GPAuxPos = parsed_args["GPAuxPos"]
NumV = parsed_args["NumV"]
NumAux = parsed_args["NumAux"]
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
dtauSmall = parsed_args["dtauSmall"]
IntMethod = parsed_args["IntMethod"]
IntMethodFast = parsed_args["IntMethodFast"]
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
ns = parsed_args["ns"]
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
OrdPolyZ = parsed_args["OrdPolyZ"]
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
# Examples
aC = parsed_args["LengthOfAgnesiHill"]
=#

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
#    AMDGPUS.allowscalar(false)
#  elseif JuliaGPU == "Metal"
#    backend = MetalBackend()
#    Metal.allowscalar(true)
  end
else
  backend = CPU()
end
FloatTypeBackend = "Float64"
if FloatTypeBackend == "Float64"
  FTB = Float64
elseif FloatTypeBackend == "Float32"
  FTB = Float32
else
  @show "False FloatTypeBackend"
  stop
end

Problem = "BaroWaveDrySphere"

Param = Examples.Parameters(FTB,Problem)

KernelAbstractions.synchronize(backend)

Parallel = true

# Physical parameters
Phys = DyCore.PhysParameters{FTB}()

RefineLevel = 6
nz = 4
nQuad = 3
nQuadM = 3 #2
nQuadS = 3 #3
Decomp = "EqualArea"
nLat = 0
nLon = 0
LatB = 0.0
H=10000.0
Lx=20000.0
Ly=2000.0
x0=0.0
y0=0.0
nz=16
nx=32
ny=2
BoundaryWE="Period" 
BoundarySN="Period" 
BoundaryBT="" 
TopoS = ""
P1 = 0.0
P2 = 0.0
P3 = 0.0
P4 = 0.0
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

#Quad
GridType = "CubedSphere"
#GridType = "TriangularSphere"
RefineLevel = 3
nPanel = 20
#GridType = "HealPix"
ns = 57
OrdPoly = 4
OrdPolyZ = 4

#Grid construction
RadEarth = Phys.RadEarth
Discretization = "DiscretizationDG"
#ModelParameters
Model = DyCore.ModelStruct{FTB}()

#Grid, CellToProc = Grids.InitGridSphere(backend,FTB,OrdPoly,nz,nPanel,RefineLevel,ns,nLon,nLat,LatB,
#    GridType,Decomp,RadEarth,Model,ParallelCom;Discretization=Discretization,ChangeOrient=2)

Grid, CellToProc = Grids.InitGridCart(backend,FTB,OrdPoly,nx,ny,Lx,Ly,x0,y0,Boundary,nz,Model,ParallelCom;
  order=true,GridType="Quad",Discretization="CG",ChangeOrient=3)


z = zeros(nz+1)
z[1] = 0.0
for i = 1 : nz
  z[i+1] = z[i] + 300.0
end  
@. Grid.z = z
Grid.nz = nz

Grid.AdaptGrid = Grids.AdaptGrid(FTB,"Sleve",FTB(z[end]))

OrdPrint = 1
OrdPrintZ = 1
TimeStepper = DyCore.TimeStepperStruct{FTB}(backend)

Output = DyCore.OutputStruct()
DoF = (OrdPoly + 1) * (OrdPoly + 1)
NumV = 5
NumTr = 0
Global = DyCore.GlobalStruct{FTB}(backend,Grid,Model,TimeStepper,ParallelCom,Output,DoF,nz,
    NumV,NumTr)
DGMethod = "Kubatko2LGL"
if Grid.Type == Grids.Quad()
  DG = FiniteElements.DGQuad{FTB}(backend,OrdPoly,OrdPolyZ,OrdPrint,OrdPrintZ,Grid,ParallelCom.Proc)
else
  DGMethod = "Kubatko2LGL"
  DG = FiniteElements.DGTri{FTB}(backend,DGMethod,OrdPolyZ,OrdPrint,OrdPrintZ,Grid,ParallelCom.Proc)  
end  


NF = Grid.NumFaces
@show Grid.Type
if Grid.Type == Grids.Tri()
  nP = 3
elseif Grid.Type == Grids.Quad()
  nP = 4
end  
F = zeros(nP,3,NF)
for iF = 1 : NF
  for i = 1 : nP
    iN = Grid.Faces[iF].N[i]  
    F[i,1,iF] = Grid.Nodes[iN].P.x
    F[i,2,iF] = Grid.Nodes[iN].P.y
    F[i,3,iF] = Grid.Nodes[iN].P.z
  end
end  
@show F[1,:,1]
@show F[2,:,1]
@show F[3,:,1]
@show size(F)


Metric = Grids.MetricStruct{FTB}(backend,DG.DoF,DG.OrdPolyZ+1,Grid.NumFaces,nz,DG.NumG)

zS = zeros(DG.DoF,NF)
Grids.JacobiDG3GPU!(Grid.AdaptGrid,Metric.X,Metric.dXdxI,Metric.J,Metric.Rotate,DG,F,z,zS,Grid.Rad,Grid.Type,Grid.Form)
@show F[1,:,1]
@show F[2,:,1]
@show F[3,:,1]
@show size(F)
NumberThreadGPU = 256
MetricNeu = FiniteElements.Metric!(backend,FTB,DG,Grid,NumberThreadGPU,zS,Grid.Type)
@show Metric.dXdxI[1,:,1,1,1,1]
@show MetricNeu.dXdxI[1,:,1,1,1,1]
@show Metric.dXdxI[2,:,1,1,1,1]
@show MetricNeu.dXdxI[2,:,1,1,1,1]
@show Metric.dXdxI[3,:,1,1,1,1]
@show MetricNeu.dXdxI[3,:,1,1,1,1]
@show "-----------------------"

#@show Metric.dXdxI[1, :, 5, 1, 4, 1811]
#@show MetricNeu.dXdxI[1, :, 5, 1, 4, 1811]
#@show Metric.dXdxI[2, :, 5, 1, 4, 1811]
#@show MetricNeu.dXdxI[2, :, 5, 1, 4, 1811]
#@show Metric.dXdxI[3, :, 5, 1, 4, 1811]
#@show MetricNeu.dXdxI[3, :, 5, 1, 4, 1811]
#@show Metric.X[1,1,:,1,1]
#@show MetricNeu.X[1,1,:,1,1]
@show sum(abs.(MetricNeu.X-Metric.X))
@show sum(abs.(MetricNeu.dXdxI-Metric.dXdxI))
@show findmax(abs.(MetricNeu.dXdxI-Metric.dXdxI))
@show findmin(abs.(MetricNeu.dXdxI-Metric.dXdxI))
@show sum(abs.(MetricNeu.J-Metric.J))
@show findmax(abs.(MetricNeu.J-Metric.J))
@show Metric.J[17, 4, 15, 32],MetricNeu.J[17, 4, 15, 32]
@show sum(abs.(MetricNeu.Rotate-Metric.Rotate))
@show sum(abs.(MetricNeu.Rotate))
@show sum(abs.(Metric.Rotate))
stop
