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

RefineLevel = 5
RadEarth = 1.0
nz = 1
nPanel = 40
nQuad = 3
Decomp = ""
Decomp = "EqualArea"

#TRI
GridType = "DelaunaySphere"
#GridType = "CubedSphere"
Grid, Exchange = Grids.InitGridSphere(backend,FTB,OrdPoly,nz,nPanel,RefineLevel,GridType,Decomp,RadEarth,
  Model,ParallelCom;order=false)
vtkSkeletonMesh = Outputs.vtkStruct{Float64}(backend,Grid)
p = ones(Grid.NumFaces,1)
for i = 1 : Grid.NumFaces
  p[i] = i
end  
FileNumber = 0
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType, Proc, ProcNumber, p, FileNumber)

KiteGrid = Grids.Grid2KiteGrid(backend,FTB,Grid,Grids.OrientFaceSphere)
vtkSkeletonKite = Outputs.vtkStruct{Float64}(backend,KiteGrid)
pKite = ones(KiteGrid.NumFaces,1)
FileNumber += 1
Outputs.vtkSkeleton!(vtkSkeletonKite, "KiteGrid", Proc, ProcNumber, pKite, FileNumber)


CG1KiteP = FEMSei.CG1KitePrimalStruct{FTB}(Grids.Quad(),backend,KiteGrid)
CG1KiteP.M = FEMSei.MassMatrix(backend,FTB,CG1KiteP,KiteGrid,1,FEMSei.Jacobi)
CG1KiteD = FEMSei.CG1KiteDualStruct{FTB}(Grids.Quad(),backend,KiteGrid)
CG1KiteD.M = FEMSei.MassMatrix(backend,FTB,CG1KiteD,KiteGrid,1,FEMSei.Jacobi)
Div = FEMSei.DivMatrix(backend,FTB,CG1KiteD,CG1KiteP,KiteGrid,3,FEMSei.Jacobi)
#Grad = FEMSei.GradMatrix(backend,FTB,CG1KiteP,CG1KiteD,KiteGrid,3,FEMSei.Jacobi)


U = zeros(FTB,CG1KiteP.NumG+CG1KiteD.NumG)
UNew = zeros(FTB,CG1KiteP.NumG+CG1KiteD.NumG)
RhoPos = 1
uPos = CG1KiteP.NumG + 1


FileNumber += 1
@views FEMSei.Project!(backend,FTB,U[RhoPos:CG1KiteP.NumG],CG1KiteP,KiteGrid,nQuad,FEMSei.Jacobi,FEMSei.fp)
VelCa = zeros(KiteGrid.NumFaces,Grid.Dim)
VelSp = zeros(KiteGrid.NumFaces,2)
@views pM = FEMSei.ComputeScalar(backend,FTB,CG1KiteP,KiteGrid,U[RhoPos:CG1KiteP.NumG])
Outputs.vtkSkeleton!(vtkSkeletonKite, "KiteGrid", Proc, ProcNumber, [pM VelCa VelSp], FileNumber)

nAdveVel = 1000
dtau = 0.001
time = 0.0
r = similar(U)
UNew = similar(U)

for i = 1 : nAdveVel
  @show i  
# @views mul!(r[RhoPos:uPos-1],Div,U[uPos:end])
  @views r[RhoPos:uPos-1] = Div * U[uPos:end]
  @views r[RhoPos:uPos-1] = CG1KiteP.M\r[RhoPos:uPos-1]
# @views mul!(r[uPos:end],Div',U[RhoPos:uPos-1])
  @views r[uPos:end] = Div' * U[RhoPos:uPos-1]
  @views r[uPos:end] = -CG1KiteD.M\r[uPos:end]

  @. UNew = U + 0.5 * dtau * r

# @views mul!(r[RhoPos:uPos-1],Div,UNew[uPos:end])
  @views r[RhoPos:uPos-1] = Div * UNew[uPos:end]
  @views r[RhoPos:uPos-1] = CG1KiteP.M\r[RhoPos:uPos-1]
# @views mul!(r[uPos:end],Div',UNew[RhoPos:uPos-1])
  @views r[uPos:end] = Div' * UNew[RhoPos:uPos-1]
  @views r[uPos:end] = -CG1KiteD.M\r[uPos:end]

  @. U = U + dtau * r

end
FileNumber += 1
@views pM = FEMSei.ComputeScalar(backend,FTB,CG1KiteP,KiteGrid,U[RhoPos:uPos-1])
@views FEMSei.ConvertVelocityCart!(backend,FTB,VelCa,U[uPos:end],CG1KiteD,KiteGrid,FEMSei.Jacobi)
@views FEMSei.ConvertVelocitySp!(backend,FTB,VelSp,U[uPos:end],CG1KiteD,KiteGrid,FEMSei.Jacobi)
Outputs.vtkSkeleton!(vtkSkeletonKite, "KiteGrid", Proc, ProcNumber, [pM VelCa VelSp], FileNumber)


nothing
