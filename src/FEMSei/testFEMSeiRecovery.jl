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
nPanel =  80
#GridType = "HealPix"
ns = 10

print("Which Problem do you want so solve? \n")
print("1 - GalewskiSphere\n\
       2 - HaurwitzSphere\n\
       3 - LinearBlob\n\
       4 - AdvectionSpherical\n")
text = readline() 
a = parse(Int,text)
if  a == 1
    Problem = "GalewskiSphere"
    Param = Examples.Parameters(FTB,Problem)
    RadEarth = Phys.RadEarth
    dtau = 30
    nAdveVel = 100 #ceil(Int,6*24*3600/dtau)
    GridTypeOut = GridType*"NonLinShallowGal"
    @show nAdveVel
elseif  a == 2
    Problem = "HaurwitzSphere"
    Param = Examples.Parameters(FTB,Problem)
    RadEarth = Phys.RadEarth
    dtau = 30 #g=9.81, H=8000
    nAdveVel = ceil(Int,6*24*3600/dtau)
    GridTypeOut = GridType*"NonLinShallowHaurwitz"
    @show nAdveVel
elseif  a == 3
    Problem = "LinearBlob"
    Param = Examples.Parameters(FTB,Problem)
    RadEarth = 1.0
    dtau = 0.00025
    nAdveVel = 100
    GridTypeOut = GridType*"NonLinShallowBlob"
    @show nAdveVel
elseif  a == 4
    Problem = "AdvectionSphereSpherical"
    Param = Examples.Parameters(FTB,Problem)
    RadEarth = 1.0
    dtau = 2*pi*RadEarth/4/nPanel/Param.uMax*0.7
    @show dtau  # 0.0004581489286485114 #in s = 2*pi*Rad / 4*nPanel / param.uMax * cFL (ca. 0.7) bei RK2 (RK3 1.7)
    nAdveVel = 100
    nprint = 10
    GridTypeOut = GridType*"Advec"
    @show nAdveVel
else 
    print("Error")
end
println("The chosen Problem is ") 
Examples.InitialProfile!(Model,Problem,Param,Phys)

#Grid construction
Grid, Exchange = Grids.InitGridSphere(backend,FTB,OrdPoly,nz,nPanel,RefineLevel,ns,nLat,nLon,LatB,GridType,Decomp,RadEarth,
  Model,ParallelCom)
for iE = 1 : Grid.NumEdges
  Grids.PosEdgeInFace!(Grid.Edges[iE],Grid.Edges,Grid.Faces)
end
vtkSkeletonMesh = Outputs.vtkStruct{Float64}(backend,Grid,Grid.NumFaces,Flat)

#finite elements
VecDG = FEMSei.VecDG0Struct{FTB}(Grids.Quad(),backend,Grid)
DG = FEMSei.DG0Struct{FTB}(Grids.Quad(),backend,Grid)
CG = FEMSei.CG1Struct{FTB}(Grids.Quad(),backend,Grid)
RT = FEMSei.RT0Struct{FTB}(Grids.Quad(),backend,Grid)

#massmatrix und LU-decomposition
VecDG.M = FEMSei.MassMatrix(backend,FTB,VecDG,Grid,nQuadM,FEMSei.Jacobi!) 
VecDG.LUM = lu(VecDG.M)
DG.M = FEMSei.MassMatrix(backend,FTB,DG,Grid,nQuadM,FEMSei.Jacobi!)
DG.LUM = lu(DG.M)
CG.M = FEMSei.MassMatrix(backend,FTB,CG,Grid,nQuadM,FEMSei.Jacobi!)
CG.LUM = lu(CG.M)
RT.M = FEMSei.MassMatrix(backend,FTB,RT,Grid,nQuadM,FEMSei.Jacobi!)
RT.LUM = lu(RT.M)

#stiffmatrix
Rhs = zeros(FTB,DG.NumG)
uHDiv = zeros(FTB,RT.NumG)
uVecDG = zeros(FTB,VecDG.NumG)
cDG = zeros(FTB,DG.NumG)

#variables for edges
VelCa = zeros(Grid.NumFaces,Grid.Dim)
VelSp = zeros(Grid.NumFaces,2)
#scalar heightmap
h = zeros(FTB,DG.NumG) #size 38400
#velocity field in HDiv-Form
u = zeros(FTB,RT.NumG)

QuadOrd=2
hCG = zeros(FTB,CG.NumG)
#calculation of hCG
FEMSei.Project!(backend,FTB,hCG,CG,Grid,nQuad,FEMSei.Jacobi!,Model.InitialProfile)
#projection from CG1 to DG1
FEMSei.ProjectScalarScalar!(backend,FTB,cDG,DG,hCG,CG,Grid,Grids.Quad(),nQuad,FEMSei.Jacobi!)
stop
#calculation of u
FEMSei.Project!(backend,FTB,u,RT,Grid,nQuad,FEMSei.Jacobi!,Model.InitialProfile)
#calculation of Tracer
FEMSei.ProjectTr!(backend,FTB,cDG,DG,Grid,nQuad,FEMSei.Jacobi!,Model.InitialProfile)
#calculation of h
FEMSei.Project!(backend,FTB,h,DG,Grid,nQuad,FEMSei.Jacobi!,Model.InitialProfile)
FEMSei.ConvertVelocitySp!(backend,FTB,VelSp,u,RT,Grid,FEMSei.Jacobi!)
#print cDG and u
FileNumber=0
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType, Proc, ProcNumber, [cDG VelSp], FileNumber)
#runge-kutta steps
time = 0.0
cDGNew = similar(cDG)

for i = 1 : nAdveVel
  @show i  
  @. Rhs = 0
  FEMSei.DivMomentumScalar!(backend,FTB,Rhs,u,RT,cDG,DG,DG,Grid,Grids.Quad(),nQuad,FEMSei.Jacobi!)
  ldiv!(DG.LUM,Rhs)
  @. cDGNew = cDG + 0.5 * dtau * Rhs
  @. Rhs = 0
  FEMSei.DivMomentumScalar!(backend,FTB,Rhs,u,RT,cDGNew,DG,DG,Grid,Grids.Quad(),nQuad,FEMSei.Jacobi!)
  ldiv!(DG.LUM,Rhs)
  @. cDG = cDG + dtau * Rhs
  if mod(i,nprint) == 0 
    global FileNumber += 1
    Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType, Proc, ProcNumber, [cDG VelSp], FileNumber)
  end
end
@show "finished"


#=
#Berechnung von hu
FEMSei.ProjecthScalaruHDivHDiv!(backend,FTB,hu,RT,h,DG,u,RT,Grid,Grids.Quad(),nQuadS,FEMSei.Jacobi!)
FEMSei.ConvertVelocitySp!(backend,FTB,VelSp,hu,RT,Grid,FEMSei.Jacobi!)
#Ausgabe von hu
FileNumber=1001
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType, Proc, ProcNumber, [h VelSp], FileNumber)

FEMSei.ProjectVectorScalarVectorHDiv(backend,FTB,uDG,VecDG,h,DG,hu,RT,Grid,Grids.Quad(),nQuadS,FEMSei.Jacobi!)

FEMSei.ConvertVelocityCart!(backend,FTB,VelCa,uDG,VecDG,Grid,FEMSei.Jacobi!)

FEMSei.ConvertVelocitySp!(backend,FTB,VelSp,uDG,VecDG,Grid,FEMSei.Jacobi!) 

#Ausgabe von uDG
#Ausgabeliste 
FileNumber=1002
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType, Proc, ProcNumber, [h VelCa VelSp], FileNumber)

FEMSei.DivMomentumScalar!(backend,FTB,Rhs,uHDiv,RT,cDG,DG,DG,Grid,Grids.Quad(),nQuad,FEMSei.Jacobi!)

stop
ModelFEM = FEMSei.ModelFEM(backend,FTB,ND,RT,DG,Grid,nQuadM,nQuadS,FEMSei.Jacobi!)

pPosS = ModelFEM.pPosS
pPosE = ModelFEM.pPosE
uPosS = ModelFEM.uPosS
uPosE = ModelFEM.uPosE
U = zeros(FTB,ModelFEM.DG.NumG+ModelFEM.RT.NumG)
@views Up = U[pPosS:pPosE]
@views Uu = U[uPosS:uPosE]

FEMSei.Project!(backend,FTB,Uu,ModelFEM.RT,Grid,nQuad, FEMSei.Jacobi!,Model.InitialProfile)
FEMSei.Project!(backend,FTB,Up,ModelFEM.DG,Grid,nQuad, FEMSei.Jacobi!,Model.InitialProfile)

FEMSei.TimeStepper(backend,FTB,U,dtau,FEMSei.FcnNonLinShallow!,ModelFEM,Grid,nQuadM,nQuadS,FEMSei.Jacobi!,nAdveVel,GridTypeOut,Proc,ProcNumber)
=#
