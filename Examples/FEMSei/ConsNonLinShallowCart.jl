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
SimTime = parsed_args["SimTime"]
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
GridForm = parsed_args["GridForm"]
GridType = parsed_args["GridType"]
RadEarth = parsed_args["RadEarth"]
# flat
nx = parsed_args["nx"]
ny = parsed_args["ny"]
nz = parsed_args["nz"]
H = parsed_args["H"]
Lx = parsed_args["Lx"]
Ly = parsed_args["Ly"]
x0 = parsed_args["x0"]
y0 = parsed_args["y0"]
BoundaryWE = parsed_args["BoundaryWE"]
BoundarySN = parsed_args["BoundarySN"]
BoundaryBT = parsed_args["BoundaryBT"]
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
PrintTime = parsed_args["PrintTime"]
PrintStartTime = parsed_args["PrintStartTime"]
Flat = parsed_args["Flat"]

# Device
Device = parsed_args["Device"]
GPUType = parsed_args["GPUType"]
FloatTypeBackend = parsed_args["FloatTypeBackend"]
NumberThreadGPU = parsed_args["NumberThreadGPU"]

# Finite elements
k = parsed_args["OrderFEM"]

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

#Grid construction
RadEarth = Phys.RadEarth

Boundary = Grids.Boundary()
Boundary.WE = BoundaryWE
Boundary.SN = BoundarySN
Boundary.BT = BoundaryBT

Grid, Exchange = Grids.InitGridCart(backend,FTB,OrdPoly,nx,ny,Lx,Ly,x0,y0,Boundary,nz,Model,ParallelCom)
@show Grid.Form
#Grid, Exchange = Grids.CartGrid(backend,FTB,nx,ny,Lx,Ly,x0,y0,Grids.OrientFaceCart,Boundary,nz)

Param = Examples.Parameters(FTB,Problem)

  Problem == "BickleyJet"
  GridLengthMin,GridLengthMax = Grids.GridLength(Grid)
  #cS = sqrt(Phys.Grav * Param.H0G)
  dtau = 0.2 * 2π / 64
  EndTime = SimTime + 3600*24*SimDays + 3600 * SimHours + 60 * SimMinutes + SimSeconds
  nAdveVel = round(EndTime / dtau)
  dtau = EndTime / nAdveVel
  PrintT = PrintTime + 3600*24*PrintDays + 3600 * PrintHours + 60 * PrintMinutes + PrintSeconds
  nprint = ceil(PrintT/dtau)
  FileNameOutput = GridType*"BickleyJet"
  FileNameOutput = "Flat/"*GridType*"BickleyJet"
  @show GridLengthMin,GridLengthMax
  @show nAdveVel
  @show dtau
  @show nprint

Examples.InitialProfile!(backend,FTB,Model,Problem,Param,Phys)

#Output
vtkSkeletonMesh = Outputs.vtkStruct{Float64}(backend,Grid,Grid.NumFaces,Flat)

#Quadrature rules
if Grid.Type == Grids.Quad()
  nQuad = 4
  nQuadM = 4
  nQuadS = 4
elseif Grid.Type == Grids.Tri()
  nQuad = 4
  nQuadM = 4
  nQuadS = 4
end  

#Finite elements
DG = FEMSei.DGStruct{FTB}(backend,k,Grid.Type,Grid)
VecDG = FEMSei.VecDGStruct{FTB}(backend,k,Grid.Type,Grid)
RT = FEMSei.RTStruct{FTB}(backend,k,Grid.Type,Grid)
ND = FEMSei.NDStruct{FTB}(backend,k,Grid.Type,Grid)

#massmatrix und LU-decomposition
DG.M = FEMSei.MassMatrix(backend,FTB,DG,Grid,nQuadM,FEMSei.JacobiCart!)
DG.LUM = lu(DG.M)
VecDG.M = FEMSei.MassMatrix(backend,FTB,VecDG,Grid,nQuadM,FEMSei.JacobiCart!)
VecDG.LUM = lu(VecDG.M)
RT.M = FEMSei.MassMatrix(backend,FTB,RT,Grid,nQuadM,FEMSei.JacobiCart!)
RT.LUM = lu(RT.M)
ND.M = FEMSei.MassMatrix(backend,FTB,ND,Grid,nQuadM,FEMSei.JacobiCart!)
ND.LUM = lu(ND.M)
Curl = FEMSei.CurlMatrix(backend,FTB,ND,DG,Grid,nQuad,FEMSei.JacobiCart!)


#Runge-Kutta steps
time = 0.0

hPosS = 1
hPosE = DG.NumG
huPosS = DG.NumG + 1
huPosE = DG.NumG + RT.NumG

U = zeros(FTB,DG.NumG+RT.NumG)
@views Uh = U[hPosS:hPosE]  
@views Uhu = U[huPosS:huPosE]  
UNew = zeros(FTB,DG.NumG+RT.NumG)
@views UNewh = UNew[hPosS:hPosE]  
@views UNewhu = UNew[huPosS:huPosE]  
F = zeros(FTB,DG.NumG+RT.NumG)
@views Fh = F[hPosS:hPosE]  
@views Fhu = F[huPosS:huPosE]  
uRec = zeros(FTB,VecDG.NumG)
cName = ["h";"Vort";"uS";"vS";"wS"]

FEMSei.InterpolatehRT!(Uhu,RT,FEMSei.JacobiCart!,Grid,Grid.Type,nQuad,Model.InitialProfile)
@show maximum(Uhu),minimum(Uhu)
FEMSei.InterpolateDG!(Uh,DG,FEMSei.JacobiCart!,Grid,Grid.Type,Model.InitialProfile)
# Output of the initial values
FileNumber=0
VelCart = zeros(Grid.NumFaces,3)
hout = zeros(Grid.NumFaces)
Vort = zeros(Grid.NumFaces)
FEMSei.ConvertScalarVelocityCart!(backend,FTB,VelCart,Uhu,RT,Uh,DG,Grid,FEMSei.JacobiCart!)
FEMSei.ConvertScalar!(backend,FTB,hout,Uh,DG,Grid,FEMSei.JacobiCart!)
FEMSei.Vorticity!(backend,FTB,Vort,DG,Uhu,RT,Uh,DG,ND,Curl,Grid,Grid.Type,nQuad,FEMSei.JacobiCart!)
Outputs.vtkSkeleton!(vtkSkeletonMesh, FileNameOutput, Proc, ProcNumber, [hout Vort VelCart;] ,FileNumber,cName)

nAdveVel = 1
nPrint = 1
for i = 1 : nAdveVel
  @show i,(i-1)*dtau/3600 
  @. F = 0  
  # Tendency h
  FEMSei.DivRhs!(backend,FTB,Fh,DG,Uhu,RT,Grid,Grid.Type,nQuad,FEMSei.JacobiCart!)
  ldiv!(DG.LUM,Fh)
  # Tendency hu
  FEMSei.InterpolateScalarHDivVecDG!(backend,FTB,uRec,VecDG,Uh,DG,Uhu,RT,Grid,
    Grid.Type,nQuad,FEMSei.JacobiCart!)
  FEMSei.DivMomentumVector!(backend,FTB,Fhu,RT,Uhu,RT,uRec,VecDG,Grid,Grid.Type,nQuad,FEMSei.JacobiCart!)
  FEMSei.CrossRhs!(backend,FTB,Fhu,RT,Uhu,RT,Grid,Grid.Type,nQuad,FEMSei.JacobiCart!)
  FEMSei.GradHeightSquared!(backend,FTB,Fhu,RT,Uh,DG,Grid,Grid.Type,nQuad,FEMSei.JacobiCart!)
  ldiv!(RT.LUM,Fhu)
  @. UNew = U + 1 / 3 * dtau * F

  @. F = 0  
  # Tendency h
  FEMSei.DivRhs!(backend,FTB,Fh,DG,UNewhu,RT,Grid,Grid.Type,nQuad,FEMSei.JacobiCart!)
  ldiv!(DG.LUM,Fh)
  # Tendency hu
  FEMSei.InterpolateScalarHDivVecDG!(backend,FTB,uRec,VecDG,UNewh,DG,UNewhu,RT,Grid,
    Grid.Type,nQuad,FEMSei.JacobiCart!)
  FEMSei.DivMomentumVector!(backend,FTB,Fhu,RT,UNewhu,RT,uRec,VecDG,Grid,Grid.Type,nQuad,FEMSei.JacobiCart!)
  FEMSei.CrossRhs!(backend,FTB,Fhu,RT,UNewhu,RT,Grid,Grid.Type,nQuad,FEMSei.JacobiCart!)
  FEMSei.GradHeightSquared!(backend,FTB,Fhu,RT,UNewh,DG,Grid,Grid.Type,nQuad,FEMSei.JacobiCart!)
  ldiv!(RT.LUM,Fhu)
  @. UNew = U + 0.5 * dtau * F

  @. F = 0  
  # Tendency h
  FEMSei.DivRhs!(backend,FTB,Fh,DG,UNewhu,RT,Grid,Grid.Type,nQuad,FEMSei.JacobiCart!)
  ldiv!(DG.LUM,Fh)
  # Tendency hu
  FEMSei.InterpolateScalarHDivVecDG!(backend,FTB,uRec,VecDG,UNewh,DG,UNewhu,RT,Grid,
    Grid.Type,nQuad,FEMSei.JacobiCart!)
  FEMSei.DivMomentumVector!(backend,FTB,Fhu,RT,UNewhu,RT,uRec,VecDG,Grid,Grid.Type,nQuad,FEMSei.JacobiCart!)
  FEMSei.CrossRhs!(backend,FTB,Fhu,RT,UNewhu,RT,Grid,Grid.Type,nQuad,FEMSei.JacobiCart!)
  FEMSei.GradHeightSquared!(backend,FTB,Fhu,RT,UNewh,DG,Grid,Grid.Type,nQuad,FEMSei.JacobiCart!)
  ldiv!(RT.LUM,Fhu)
  @. U = U + dtau * F
    
  # Output
  if mod(i,nprint) == 0 
    global FileNumber += 1
    FEMSei.ConvertScalarVelocitySp!(backend,FTB,VelSp,Uhu,RT,Uh,DG,Grid,FEMSei.JacobiCart!)
    FEMSei.ConvertScalar!(backend,FTB,hout,Uh,DG,Grid,FEMSei.JacobiCart!)
    FEMSei.Vorticity!(backend,FTB,Vort,DG,Uhu,RT,Uh,DG,ND,Curl,Grid,Grid.Type,nQuad,FEMSei.JacobiCart!)
    Outputs.vtkSkeleton!(vtkSkeletonMesh, FileNameOutput, Proc, ProcNumber, [hout Vort VelSp] ,FileNumber,cName)
  end
end
@show "finished"

