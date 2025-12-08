import CGDycore:
  Examples, Parallels, Models, Grids, Outputs, Integration,  GPUS, DyCore, FEM, FiniteVolumes
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
RefineLevel = parsed_args["RefineLevel"]
nLon = parsed_args["nLon"]
nLat = parsed_args["nLat"]
LatB = parsed_args["LatB"]
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
ParallelCom.ProcNumber = ProcNumber

# Physical parameters
Phys = DyCore.PhysParameters{FTB}()

#ModelParameters
Model = DyCore.ModelStruct{FTB}()

#Problem = "GalewskySphere"
#Problem = "HaurwitzSphere"
RadEarth = Phys.RadEarth
dtau = 6
nAdveVel = 16000
#=
nAdveVel = 5000
Problem = "LinearBlob"
Fac = 1.0
RadEarth = 1.0 * Fac
dtau = 0.0001 * Fac
nAdveVel = 16000
PrintStp = 800
=#
Flat = false


Param = Examples.Parameters(FTB,Problem)
Examples.InitialProfile!(backend,FTB,Model,Problem,Param,Phys)

RefineLevel = 6
nz = 1
nPanel = 100
nQuad = 10
ns = 120
Decomp = ""
Decomp = "EqualArea"
Model.HorLimit = true
OrdPoly = 1

#TRI
#GridType = "TriangularSphere"
#GridType = "DelaunaySphere"
GridType = "CubedSphere"
#GridType = "HealPix"
#GridType = "MPAS"
Grid, Exchange = Grids.InitGridSphere(backend,FTB,OrdPoly,nz,nPanel,RefineLevel,ns,nLon,nLat,LatB,GridType,Decomp,RadEarth,
# Model,ParallelCom;order=false,ChangeOrient=2)
  Model,ParallelCom;order=false)

vtkSkeletonMesh = Outputs.vtkStruct{Float64}(backend,Grid,Grid.NumFaces,Flat)
vtkSkeletonMeshGhost = Outputs.vtkStruct{Float64}(backend,Grid,Grid.NumFaces+Grid.NumFacesG,Flat)

MetricFV = FiniteVolumes.MetricFiniteVolume(backend,FTB,Grid)

Cache = FiniteVolumes.CacheFV(backend,FTB,MetricFV,Grid)

U = zeros(FTB,Grid.NumFaces+Grid.NumFacesG+Grid.NumEdges)
r = similar(U)
UNew = similar(U)
pPosS = 1
pPosEI = Grid.NumFaces
pPosE = Grid.NumFaces + Grid.NumFacesG
uPosS = pPosE + 1
uPosE = pPosE + Grid.NumEdges
@views Up = U[pPosS:pPosE]
@views UpI = U[pPosS:pPosEI]
@views Uu = U[uPosS:uPosE]
@views rp = r[pPosS:pPosE]
@views rpI = r[pPosS:pPosEI]
@views ru = r[uPosS:uPosE]

FiniteVolumes.ProjectFace!(backend,FTB,UpI,Grid,Model.InitialProfile)
FiniteVolumes.ProjectEdge!(backend,FTB,Uu,Grid,Model.InitialProfile)
FileNumber = 0
VelCa = zeros(Grid.NumFaces,Grid.Dim)
VelSp = zeros(Grid.NumFaces,2)
FiniteVolumes.ConvertVelocityCart!(backend,FTB,VelCa,Uu,Grid)
FiniteVolumes.ConvertVelocitySp!(backend,FTB,VelSp,Uu,Grid)
SpeciesList = ["h","uC","vC","wC","uS","vS"]
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType*"FV", Proc, ProcNumber, [UpI VelCa VelSp], FileNumber, SpeciesList)
FileNumber += 1

time = 0.0

nAdveVel = 17280/2
dtau = 30.0
PrintStp = 1440/2
PrintStp = 200
nAdveVel=8000
#FiniteVolumes.FcnFV!(r,U,MetricFV,Grid,Cache,Phys)
#    FiniteVolumes.ConvertVelocityCart!(backend,FTB,VelCa,ru,Grid)
#    FiniteVolumes.ConvertVelocitySp!(backend,FTB,VelSp,ru,Grid)
#    Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType*"FVKC", Proc, ProcNumber, [rpI VelCa VelSp], FileNumber,SpeciesList)
#    stop
for i = 1 : nAdveVel
  FiniteVolumes.FcnFV!(r,U,MetricFV,Grid,Cache,Phys)
  @. UNew = U + 1/3 * dtau * r  
  FiniteVolumes.FcnFV!(r,UNew,MetricFV,Grid,Cache,Phys)
  @. UNew = U + 0.5 * dtau * r  
  FiniteVolumes.FcnFV!(r,UNew,MetricFV,Grid,Cache,Phys)
  @. U = U + dtau * r  
  if mod(i,PrintStp) == 0
    @show dtau*i  
    FiniteVolumes.ConvertVelocityCart!(backend,FTB,VelCa,Uu,Grid)
    FiniteVolumes.ConvertVelocitySp!(backend,FTB,VelSp,Uu,Grid)
    Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType*"FV", Proc, ProcNumber, [UpI VelCa VelSp], FileNumber,SpeciesList)
    global FileNumber += 1
  end  

  #time = time + dtau
end


FileNumber += 1
FiniteVolumes.ConvertVelocityCart!(backend,FTB,VelCa,Uu,Grid)
FiniteVolumes.ConvertVelocitySp!(backend,FTB,VelSp,Uu,Grid)
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType*"FV", Proc, ProcNumber, [UpI VelCa VelSp], FileNumber,SpeciesList)


#= ToDO
Gleichungen testen mit @show
mehr Quadraturregeln für Dreiecke höhere ordnungen
RT1, DG1
Wiki von github
=#
