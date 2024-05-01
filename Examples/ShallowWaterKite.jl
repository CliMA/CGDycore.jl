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

RefineLevel = 4
nz = 1
nPanel = 25
nQuad = 3
Decomp = ""
Decomp = "EqualArea"
Problem = "GalewskiSphere"
RadEarth = Phys.RadEarth
dtau = 50
nAdveVel = 60
Problem = "LinearBlob"
RadEarth = 1.0
dtau = 0.00025
nAdveVel = 8000
Param = Examples.Parameters(FTB,Problem)
Examples.InitialProfile!(Model,Problem,Param,Phys)

#TRI
#GridType = "DelaunaySphere"
GridType = "CubedSphere"
#GridType = "TriangularSphere" # Achtung Orientierung
Grid, Exchange = Grids.InitGridSphere(backend,FTB,OrdPoly,nz,nPanel,RefineLevel,GridType,Decomp,RadEarth,
  Model,ParallelCom;order=false)
vtkSkeletonMesh = Outputs.vtkStruct{Float64}(backend,Grid)
p = ones(Grid.NumFaces,1)
for i = 1 : Grid.NumFaces
  p[i] = i
end  
global FileNumber = 0
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType, Proc, ProcNumber, p, FileNumber)

KiteGrid = Grids.Grid2KiteGrid(backend,FTB,Grid,Grids.OrientFaceSphere)
vtkSkeletonKite = Outputs.vtkStruct{Float64}(backend,KiteGrid)


CG1KiteP = FEMSei.CG1KitePrimalStruct{FTB}(Grids.Quad(),backend,KiteGrid)
CG1KiteP.M = FEMSei.MassMatrix(backend,FTB,CG1KiteP,KiteGrid,nQuad,FEMSei.Jacobi!)
FpM = lu(CG1KiteP.M)

CG1KiteDHDiv = FEMSei.CG1KiteDualHDiv{FTB}(Grids.Quad(),backend,KiteGrid)
CG1KiteDHDiv.M = FEMSei.MassMatrix(backend,FTB,CG1KiteDHDiv,KiteGrid,nQuad,FEMSei.Jacobi!)
FuM = lu(CG1KiteDHDiv.M)
Div = FEMSei.DivMatrix(backend,FTB,CG1KiteDHDiv,CG1KiteP,KiteGrid,nQuad,FEMSei.Jacobi!)
Grad = FEMSei.GradMatrix(backend,FTB,CG1KiteP,CG1KiteDHDiv,KiteGrid,nQuad,FEMSei.Jacobi!)
@show sum(abs.(Div'-Grad)),sum(abs.(Div'+Grad))

CG1KiteDHCurl = FEMSei.CG1KiteDualHCurl{FTB}(Grids.Quad(),backend,KiteGrid)
CG1KiteDHCurl.M = FEMSei.MassMatrix(backend,FTB,CG1KiteDHCurl,KiteGrid,nQuad,FEMSei.Jacobi!)
Curl = FEMSei.CurlMatrix(backend,FTB,CG1KiteDHCurl,CG1KiteP,KiteGrid,nQuad,FEMSei.Jacobi!)

pPosS = 1
pPosE = CG1KiteP.NumG
uPosS = CG1KiteP.NumG + 1
uPosE = CG1KiteP.NumG + CG1KiteDHDiv.NumG
U = zeros(FTB,CG1KiteP.NumG+CG1KiteDHDiv.NumG)
@views Up = U[pPosS:pPosE]
@views Uu = U[uPosS:uPosE]
UNew = zeros(FTB,CG1KiteP.NumG+CG1KiteDHDiv.NumG)
@views UNewp = U[pPosS:pPosE]
@views UNewu = U[uPosS:uPosE]
F = zeros(FTB,CG1KiteP.NumG+CG1KiteDHDiv.NumG)
@views Fp = F[pPosS:pPosE]
@views Fu = F[uPosS:uPosE]

pM = zeros(KiteGrid.NumFaces)
VelCa = zeros(KiteGrid.NumFaces,Grid.Dim)
VelSp = zeros(KiteGrid.NumFaces,2)

FEMSei.Project!(backend,FTB,Uu,CG1KiteDHDiv,KiteGrid,nQuad, FEMSei.Jacobi!,Model.InitialProfile)
FEMSei.Project!(backend,FTB,Up,CG1KiteP,KiteGrid,nQuad, FEMSei.Jacobi!,Model.InitialProfile)
pM = FEMSei.ComputeScalar(backend,FTB,CG1KiteP,KiteGrid,Up)
FEMSei.ConvertVelocityCart!(backend,FTB,VelCa,Uu,CG1KiteDHDiv,KiteGrid,FEMSei.Jacobi!)
FEMSei.ConvertVelocitySp!(backend,FTB,VelSp,Uu,CG1KiteDHDiv,KiteGrid,FEMSei.Jacobi!)
Outputs.vtkSkeleton!(vtkSkeletonKite, "KiteGrid", Proc, ProcNumber, [pM VelCa VelSp], FileNumber)
FileNumber += 1


UU = similar(U)
@views UUp = UU[pPosS:pPosE]
@. UUp = 1
@views UUu = UU[uPosS:uPosE]
@. UUu = 0

UCurl = similar(U)
@views qVort = UU[pPosS:pPosE]
@views UCurl = UU[uPosS:uPosE]

time = 0.0


for i = 1 : nAdveVel
  @. F = 0  
# FEMSei.DivRhs!(backend,FTB,Fp,Uu,CG1KiteDHDiv,UUp,CG1KiteP,CG1KiteP,Grid,Grids.Quad(),nQuad,FEMSei.Jacobi!)
  mul!(Fp,Div,Uu)
  ldiv!(FpM,Fp)

# @views FEMSei.GradRhs!(backend,FTB,Fu,Up,CG1KiteP,
# UUu,CG1KiteDHDiv,CG1KiteDHDiv,KiteGrid,Grids.Quad(),nQuad,FEMSei.Jacobi!)
# @views FEMSei.ProjectHDivHCurl!(backend,FTB,UCurl,CG1KiteDHCurl,Uu,CG1KiteDHDiv,
#   Grid,Grids.Quad(),nQuad,FEMSei.Jacobi!)
# mul!(qVort,Curl,UCurl)
# @. qVort = 0
# @time @views FEMSei.CrossRhs!(backend,FTB,Fu,qVort,CG1KiteP,
# Uu,CG1KiteDHDiv,CG1KiteDHDiv,KiteGrid,Grids.Quad(),nQuad,FEMSei.Jacobi!)
  mul!(Fu,Div',Up)
  ldiv!(FuM,Fu)

  @. Fu = -Fu
  @. UNew = U + 0.5 * dtau * F

  @. F = 0  
# FEMSei.DivRhs!(backend,FTB,Fp,UNewu,CG1KiteDHDiv,UUp,CG1KiteP,CG1KiteP,Grid,Grids.Quad(),nQuad,FEMSei.Jacobi!)
  mul!(Fp,Div,UNewu)
  ldiv!(FpM,Fp)

# @views FEMSei.GradRhs!(backend,FTB,Fu,UNewp,CG1KiteP,
# UUu,CG1KiteDHDiv,CG1KiteDHDiv,KiteGrid,Grids.Quad(),nQuad,FEMSei.Jacobi!)
# @views FEMSei.ProjectHDivHCurl!(backend,FTB,UCurl,CG1KiteDHCurl,UNewu,CG1KiteDHDiv,
#   Grid,Grids.Quad(),nQuad,FEMSei.Jacobi!)
# mul!(qVort,Curl,UCurl)
# @. qVort = 0
# @time @views FEMSei.CrossRhs!(backend,FTB,Fu,qVort,CG1KiteP,
# Uu,CG1KiteDHDiv,CG1KiteDHDiv,KiteGrid,Grids.Quad(),nQuad,FEMSei.Jacobi!)
  mul!(Fu,Div',UNewp)
  ldiv!(FuM,Fu)

  @. Fu = -Fu
  @. U = U + dtau * F

  if mod(i,400) == 0
    global FileNumber += 1
    pM = FEMSei.ComputeScalar(backend,FTB,CG1KiteP,KiteGrid,Up)
    FEMSei.ConvertVelocityCart!(backend,FTB,VelCa,Uu,CG1KiteDHDiv,KiteGrid,FEMSei.Jacobi!)
    FEMSei.ConvertVelocitySp!(backend,FTB,VelSp,Uu,CG1KiteDHDiv,KiteGrid,FEMSei.Jacobi!)
    Outputs.vtkSkeleton!(vtkSkeletonKite, "KiteGrid", Proc, ProcNumber, [pM VelCa VelSp], FileNumber)
  end
end
pM = FEMSei.ComputeScalar(backend,FTB,CG1KiteP,KiteGrid,Up)
FEMSei.ConvertVelocityCart!(backend,FTB,VelCa,Uu,CG1KiteDHDiv,KiteGrid,FEMSei.Jacobi!)
FEMSei.ConvertVelocitySp!(backend,FTB,VelSp,Uu,CG1KiteDHDiv,KiteGrid,FEMSei.Jacobi!)
Outputs.vtkSkeleton!(vtkSkeletonKite, "KiteGrid", Proc, ProcNumber, [pM VelCa VelSp], FileNumber)
FileNumber += 1

nothing






