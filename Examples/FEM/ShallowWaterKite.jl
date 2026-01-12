import CGDycore:
  Examples, Parallels, Models, Grids, Outputs, Integration,  GPU, DyCore, FEM, FiniteVolumes
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
ns = parsed_args["ns"]
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
ParallelCom.ProcNumber  = ProcNumber

# Physical parameters
Phys = DyCore.PhysParameters{FTB}()

#ModelParameters
Model = DyCore.ModelStruct{FTB}()

RefineLevel = 5
nz = 1
nPanel = 20
nQuad = 3
Decomp = ""
Decomp = "EqualArea"
#Problem = "GalewskySphere"
RadEarth = Phys.RadEarth
#Problem = "LinearBlob"
@show Problem
Param = Examples.Parameters(FTB,Problem)
Examples.InitialProfile!(backend,FTB,Model,Problem,Param,Phys)

#TRI
#GridType = "DelaunaySphere"
GridType = "CubedSphere"
#GridType = "TriangularSphere" # Achtung Orientierung
Grid, Exchange = Grids.InitGridSphere(backend,FTB,OrdPoly,nz,nPanel,RefineLevel,ns,
  nLat,nLon,LatB,GridType,Decomp,RadEarth,Model,ParallelCom;ChangeOrient=2)
vtkSkeletonMesh = Outputs.vtkStruct{Float64}(backend,Grid,Grid.NumFaces,Flat)
p = ones(Grid.NumFaces,1)
for i = 1 : Grid.NumFaces
  p[i] = i
end  
global FileNumber = 0
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType, Proc, ProcNumber, p, FileNumber, ["p"])

KiteGrid = Grids.Grid2KiteGrid(backend,FTB,Grid,Grids.OrientFaceSphere)
vtkSkeletonKite = Outputs.vtkStruct{Float64}(backend,KiteGrid,KiteGrid.NumFaces,Flat)


CG1KiteP = FEM.CG1KitePrimalStruct{FTB}(Grids.Quad(),backend,KiteGrid)
CG1KiteP.M = FEM.MassMatrix(backend,FTB,CG1KiteP,KiteGrid,2,FEM.Jacobi!)
CG1KiteP.LUM = lu(CG1KiteP.M)

CG1KiteDHDiv = FEM.CG1KiteDualHDiv{FTB}(Grids.Quad(),backend,KiteGrid)
CG1KiteDHDiv.M = FEM.MassMatrix(backend,FTB,CG1KiteDHDiv,KiteGrid,nQuad,FEM.Jacobi!)
CG1KiteDHDiv.LUM = lu(CG1KiteDHDiv.M)

Div = FEM.DivMatrix(backend,FTB,CG1KiteDHDiv,CG1KiteP,KiteGrid,nQuad,FEM.Jacobi!)
Grad = FEM.GradMatrix(backend,FTB,CG1KiteP,CG1KiteDHDiv,KiteGrid,nQuad,FEM.Jacobi!)
@show sum(abs.(Div'-Grad)),sum(abs.(Div'+Grad))

CG1KiteDHCurl = FEM.CG1KiteDualHCurl{FTB}(Grids.Quad(),backend,KiteGrid)
CG1KiteDHCurl.M = FEM.MassMatrix(backend,FTB,CG1KiteDHCurl,KiteGrid,nQuad,FEM.Jacobi!)
Curl = FEM.CurlMatrix(backend,FTB,CG1KiteDHCurl,CG1KiteP,KiteGrid,nQuad,FEM.Jacobi!)

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
F1 = zeros(FTB,CG1KiteP.NumG+CG1KiteDHDiv.NumG)
@views F1p = F1[pPosS:pPosE]
@views F1u = F1[uPosS:uPosE]

pM = zeros(KiteGrid.NumFaces)
VelCa = zeros(KiteGrid.NumFaces,Grid.Dim)
VelSp = zeros(KiteGrid.NumFaces,2)

FEM.Project!(backend,FTB,Uu,CG1KiteDHDiv,KiteGrid,nQuad, FEM.Jacobi!,Model.InitialProfile)
FEM.Project!(backend,FTB,Up,CG1KiteP,KiteGrid,nQuad, FEM.Jacobi!,Model.InitialProfile)
pM = FEM.ComputeScalar(backend,FTB,CG1KiteP,KiteGrid,Up)
FEM.ConvertVelocityCart!(backend,FTB,VelCa,Uu,CG1KiteDHDiv,KiteGrid,FEM.Jacobi!)
FEM.ConvertVelocitySp!(backend,FTB,VelSp,Uu,CG1KiteDHDiv,KiteGrid,FEM.Jacobi!)
Outputs.vtkSkeleton!(vtkSkeletonKite, "KiteGrid", Proc, ProcNumber, [pM VelCa VelSp], FileNumber,["p","uC","vC","wC","uS","vS"])
FileNumber += 1


UU = similar(U)
@views qVort = UU[pPosS:pPosE]
@views UCurl = UU[uPosS:uPosE]

time = 0.0
dtau = 100
nAdveVel = 1

qVort = similar(Fp)
for i = 1 : nAdveVel
  @show i  
  @. F = 0  
  FEM.DivRhs!(backend,FTB,Fp,Uu,CG1KiteDHDiv,Up,CG1KiteP,CG1KiteP,KiteGrid,Grids.Quad(),nQuad,FEM.Jacobi!)

  @show maximum(Uu),minimum(Uu)
  FEM.CurlVel!(qVort,CG1KiteP,Uu,CG1KiteDHDiv,nQuad,Grids.Quad(),KiteGrid,FEM.Jacobi!)
  @show maximum(qVort),minimum(qVort)
pM1 = FEM.ComputeScalar(backend,FTB,CG1KiteP,KiteGrid,qVort)
  FEM.CrossRhs!(backend,FTB,Fu,qVort,CG1KiteP,Uu,CG1KiteDHDiv,CG1KiteDHDiv,KiteGrid,Grids.Quad(),nQuad,FEM.Jacobi!)
  FEM.GradRhs!(backend,FTB,Fu,Up,CG1KiteP,Uu,CG1KiteDHDiv,CG1KiteDHDiv,KiteGrid,Grids.Quad(),nQuad,FEM.Jacobi!)
  ldiv!(CG1KiteDHDiv.LUM,Fu)

  @. UNew = U + 1/3 * dtau * F
FEM.ConvertVelocityCart!(backend,FTB,VelCa,UNewu,CG1KiteDHDiv,KiteGrid,FEM.Jacobi!)
FEM.ConvertVelocitySp!(backend,FTB,VelSp,UNewu,CG1KiteDHDiv,KiteGrid,FEM.Jacobi!)
Outputs.vtkSkeleton!(vtkSkeletonKite, "KiteGrid", Proc, ProcNumber, [pM1 VelCa VelSp], FileNumber,["p","uC","vC","wC","uS","vS"])
  @show FileNumber
  stop


    @views FEM.GradRhs!(backend,FTB,Fu,Up,CG1KiteP,
  UUu,CG1KiteDHDiv,CG1KiteDHDiv,KiteGrid,Grids.Quad(),nQuad,FEM.Jacobi!)
# @views FEM.ProjectHDivHCurl!(backend,FTB,UCurl,CG1KiteDHCurl,Uu,CG1KiteDHDiv,
#   Grid,Grids.Quad(),nQuad,FEM.Jacobi!)
# mul!(qVort,Curl,UCurl)
# @. qVort = 0
# @time @views FEM.CrossRhs!(backend,FTB,Fu,qVort,CG1KiteP,
# Uu,CG1KiteDHDiv,CG1KiteDHDiv,KiteGrid,Grids.Quad(),nQuad,FEM.Jacobi!)

  @. F = 0  
  FEM.DivRhs!(backend,FTB,Fp,CG1KiteP,UNewu,CG1KiteDHDiv,KiteGrid,Grids.Quad(),nQuad,FEM.Jacobi!)
  ldiv!(CG1KiteP.LUM,Fp)

  FEM.CrossRhs!(backend,FTB,Fu,CG1KiteDHDiv,UNewu,CG1KiteDHDiv,KiteGrid,Grids.Quad(),nQuad,FEM.Jacobi!)
  FEM.GradRhs!(backend,FTB,Fu,CG1KiteDHDiv,UNewp,CG1KiteP,KiteGrid,Grids.Quad(),nQuad,FEM.Jacobi!)
  ldiv!(CG1KiteDHDiv.LUM,Fu)

  @. UNew = U + 0.5 * dtau * F

  @. F = 0  
  FEM.DivRhs!(backend,FTB,Fp,CG1KiteP,UNewu,CG1KiteDHDiv,KiteGrid,Grids.Quad(),nQuad,FEM.Jacobi!)
  ldiv!(CG1KiteP.LUM,Fp)

  FEM.CrossRhs!(backend,FTB,Fu,CG1KiteDHDiv,UNewu,CG1KiteDHDiv,KiteGrid,Grids.Quad(),nQuad,FEM.Jacobi!)
  @views FEM.GradRhs!(backend,FTB,Fu,CG1KiteDHDiv,UNewp,CG1KiteP,KiteGrid,Grids.Quad(),nQuad,FEM.Jacobi!)
  ldiv!(CG1KiteDHDiv.LUM,Fu)

  @. U = U + dtau * F

  if mod(i,1) == 0
    global FileNumber += 1
    pM1 = FEM.ComputeScalar(backend,FTB,CG1KiteP,KiteGrid,Up)
    FEM.ConvertVelocityCart!(backend,FTB,VelCa,Uu,CG1KiteDHDiv,KiteGrid,FEM.Jacobi!)
    FEM.ConvertVelocitySp!(backend,FTB,VelSp,Uu,CG1KiteDHDiv,KiteGrid,FEM.Jacobi!)
    Outputs.vtkSkeleton!(vtkSkeletonKite, "KiteGrid", Proc, ProcNumber, [pM1 VelCa VelSp], FileNumber,["p","uC","vC","wC","uS","vS"])
  end
end
pM = FEM.ComputeScalar(backend,FTB,CG1KiteP,KiteGrid,Up)
FEM.ConvertVelocityCart!(backend,FTB,VelCa,Uu,CG1KiteDHDiv,KiteGrid,FEM.Jacobi!)
FEM.ConvertVelocitySp!(backend,FTB,VelSp,Uu,CG1KiteDHDiv,KiteGrid,FEM.Jacobi!)
Outputs.vtkSkeleton!(vtkSkeletonKite, "KiteGrid", Proc, ProcNumber, [pM VelCa VelSp], FileNumber,["p","uC","vC","wC","uS","vS"])
FileNumber += 1

nothing






