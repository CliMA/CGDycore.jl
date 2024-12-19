import CGDycore:
  Examples, Parallels, Models, Grids, Outputs, Integration,  GPU, DyCore, FEMSei, FiniteVolumes
using MPI
using Base
using CUDA
using AMDGPU
using Metal
using KernelAbstractions
using StaticArrays
using ArgParse
using LinearAlgebra

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
Relax = parsed_args["Relax"]
StrideDamp = parsed_args["StrideDamp"]
Geos = parsed_args["Geos"]
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
HyperDDivW = parsed_args["HyperDDivW"]
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
Flat = false #testen
Device = "CPU"
FloatTypeBackend = "Float64"

if Device == "CPU" 
  backend = CPU()
elseif Device == "GPU" 
  if GPUType == "CUDA"
    backend = CUDABackend()
    CUDA.allowscalar(true)
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

MPI.Init()
comm = MPI.COMM_WORLD
Proc = MPI.Comm_rank(comm) + 1
ProcNumber = MPI.Comm_size(comm)
ParallelCom = DyCore.ParallelComStruct()
ParallelCom.Proc = Proc
ParallelCom.ProcNumber  = ProcNumber

# Physical parameters
Phys = DyCore.PhysParameters{FTB}()

#ModelParameters
Model = DyCore.ModelStruct{FTB}()

RefineLevel = 6
nz = 1
nQuad = 3
nQuadM = 3 #2
nQuadS = 3 #3
Decomp = "EqualArea"
nLat = 0
nLon = 0
LatB = 0.0

#Quad
GridType = "CubedSphere"
nPanel = 30
#GridType = "HealPix"
ns = 57
OrdPoly = 3
OrdPolyZ = 2

#Grid construction
RadEarth = Phys.RadEarth
Grid, Exchange = Grids.InitGridSphere(backend,FTB,OrdPoly,nz,nPanel,RefineLevel,ns,
  nLat,nLon,LatB,GridType,Decomp,RadEarth,Model,ParallelCom)
H = 20000.0
Grid.AdaptGrid = Grids.AdaptGrid(FTB,"Sleve",FTB(H))

print("Which Problem do you want so solve? \n")
print("1 - GalewskiSphere\n\
       2 - HaurwitzSphere\n")

#text = readline() 
#a = parse(Int,text)
a = 1
if  a == 1
    Problem = "GalewskiSphere"
    Param = Examples.Parameters(FTB,Problem)
    GridLengthMin,GridLengthMax = Grids.GridLength(Grid)
    cS = sqrt(Phys.Grav * Param.H0G)
    dtau = GridLengthMin / cS / sqrt(2) * .5 
    EndTime = 24 * 3600 * 4# One day
    PrintTime = 3600
    nAdveVel = round(EndTime / dtau)
    dtau = EndTime / nAdveVel
    nprint = ceil(PrintTime/dtau)
    GridTypeOut = GridType*"NonLinShallowGal"
    @show GridLengthMin,GridLengthMax
    @show nAdveVel
    @show dtau
    @show nprint
elseif  a == 2
    Problem = "HaurwitzSphere"
    Param = Examples.Parameters(FTB,Problem)
    GridLengthMin,GridLengthMax = Grids.GridLength(Grid)
    cS = sqrt(Phys.Grav * Param.h0)
    dtau = GridLengthMin / cS / sqrt(2) * .3 
    EndTime = 24 * 3600 # One day
    PrintTime = 3600
    nAdveVel = round(EndTime / dtau)
    dtau = EndTime / nAdveVel
    nprint = ceil(PrintTime/dtau)
    GridTypeOut = GridType*"NonLinShallowHaurwitz"
    @show GridLengthMin,GridLengthMax
    @show nAdveVel
    @show dtau
    @show nprint
else 
    print("Error")
end
println("The chosen Problem is ") 
Examples.InitialProfile!(Model,Problem,Param,Phys)

nz = Grid.nz
Proc = ParallelCom.Proc
ProcNumber = ParallelCom.ProcNumber

TimeStepper = DyCore.TimeStepperStruct{FTB}(backend)

Output = DyCore.OutputStruct()
DoF = (OrdPoly + 1) * (OrdPoly + 1)
Global = DyCore.GlobalStruct{FTB}(backend,Grid,Model,TimeStepper,ParallelCom,Output,DoF,nz,
  Model.NumV,Model.NumTr)

DG = DyCore.DGQuad{FTB}(backend,OrdPoly,OrdPolyZ,Global.Grid)


zS = zeros(1)
Grids.AddVerticalGrid!(Grid,nz,H)
DG, Metric = DyCore.DiscretizationDG(backend,FTB,Grids.JacobiSphere3,DG,Exchange,Global,zS)

Global.Output.Flat = true
Global.Output.OrdPrint = OrdPoly
Global.Output.OrdPrintZ = OrdPolyZ
Global.Output.RadPrint = H
Global.Output.H = H
Global.vtkCache = Outputs.vtkStruct{FTB}(backend,Global.Output.OrdPrint,Grids.TransSphereX!,DG,Metric,Global)


OP = OrdPoly + 1
OPZ = OrdPolyZ + 1
NE = Grid.NumEdges
Model.NumV = 5
Model.RhoPos = 1
Model.uPos = 2
Model.vPos = 3
Model.wPos = 4
Model.ThPos = 5
Examples.InitialProfile!(Model,Problem,Param,Phys)
@show size(Metric.X)
U = GPU.InitialConditions(backend,FTB,DG,Metric,Phys,Global,Model.InitialProfile,Param)
#Output
Global.Output.vtkFileName = "TestDG"
Global.Output.vtk = 0
Global.Output.nPanel = nPanel
Global.Output.cNames =
    [
      "Rho",
      "u",
      "v",
#     "wB",
#     "Th",
#     "Pres",
      ]
Outputs.unstructured_vtkSphere(U,Grids.TransSphereX,DG,Metric,Phys,Global,Proc,ProcNumber)      

F = similar(U)
NzG = min(div(NumberThreadGPU,OP*OPZ),nz)
group = (OPZ,OP,NzG)
ndrange = (OPZ,OP,nz,NE)
KRiemannNonLinH! = GPU.RiemannNonLinH!(backend,group)

KRiemannNonLinH!(F,U,Metric.NH,Metric.T1H,Metric.T1H,Metric.VolSurfH,
  Grid.EF,Grid.FE,DG.Glob,ndrange=ndrange)
A = 3
