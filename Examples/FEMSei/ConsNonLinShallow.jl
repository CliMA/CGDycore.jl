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
PrintTime = parsed_args["PrintTime"]
PrintStartTime = parsed_args["PrintStartTime"]
Flat = parsed_args["Flat"]
vtkFileName = parsed_args["vtkFileName"]

# Device
Device = parsed_args["Device"]
GPUType = parsed_args["GPUType"]
FloatTypeBackend = parsed_args["FloatTypeBackend"]
NumberThreadGPU = parsed_args["NumberThreadGPU"]

# Finite elements
k = parsed_args["OrderFEM"]

# Grid Output Refine
ref = parsed_args["RefineOutput"]

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
Grid, Exchange = Grids.InitGridSphere(backend,FTB,OrdPoly,nz,nPanel,RefineLevel,ns,
  nLat,nLon,LatB,GridType,Decomp,RadEarth,Model,ParallelCom;ChangeOrient=3)

Param = Examples.Parameters(FTB,Problem)

if Problem == "GalewskySphere"
  GridLengthMin,GridLengthMax = Grids.GridLength(Grid)
  cS = sqrt(Phys.Grav * Param.H0G)
  dtau = GridLengthMin / cS / sqrt(2) * .05 / (k + 1)
  EndTime = SimTime + 3600*24*SimDays + 3600 * SimHours + 60 * SimMinutes + SimSeconds
  nAdveVel::Int = round(EndTime / dtau)
  dtau = EndTime / nAdveVel
  PrintT = PrintTime + 3600*24*PrintDays + 3600 * PrintHours + 60 * PrintMinutes + PrintSeconds
  nprint::Int = ceil(PrintT/dtau)
  FileNameOutput = vtkFileName
  @show GridLengthMin,GridLengthMax
  @show nAdveVel
  @show dtau
  @show nprint
elseif Problem == "HaurwitzSphere"
  GridLengthMin,GridLengthMax = Grids.GridLength(Grid)
  cS = sqrt(Phys.Grav * Param.h0)
  dtau = GridLengthMin / cS / sqrt(2) * .1 / (k + 1)
  EndTime = SimTime + 3600*24*SimDays + 3600 * SimHours + 60 * SimMinutes + SimSeconds
  nAdveVel::Int = round(EndTime / dtau)
  dtau = EndTime / nAdveVel
  PrintT = PrintTime + 3600*24*PrintDays + 3600 * PrintHours + 60 * PrintMinutes + PrintSeconds
  nprint::Int = ceil(PrintT/dtau)
  FileNameOutput = vtkFileName
  @show GridLengthMin,GridLengthMax
  @show nAdveVel
  @show dtau
  @show nprint
else 
    print("Error")
end

Examples.InitialProfile!(backend,FTB,Model,Problem,Param,Phys)

#Output
vtkSkeletonMesh = Outputs.vtkStruct{Float64}(backend,Grid,Grid.NumFaces,Flat;Refine=ref)

#Quadrature rules
if Grid.Type == Grids.Quad()
  nQuad = 2
  nQuadM = 3
  nQuadS = 3
elseif Grid.Type == Grids.Tri()
  nQuad = 4
  nQuadM = 6
  nQuadS = 4
end  

#Finite elements
DG = FEMSei.DGStruct{FTB}(backend,k,Grid.Type,Grid)
VecDG = FEMSei.VecDGStruct{FTB}(backend,k+1,Grid.Type,Grid)
RT = FEMSei.RTStruct{FTB}(backend,k,Grid.Type,Grid)
ND = FEMSei.NDStruct{FTB}(backend,k,Grid.Type,Grid)
CG = FEMSei.CGStruct{FTB}(backend,k+1,Grid.Type,Grid)

#massmatrix und LU-decomposition
DG.M = FEMSei.MassMatrix(backend,FTB,DG,Grid,nQuadM,FEMSei.Jacobi!)
DG.LUM = lu(DG.M)
CG.M = FEMSei.MassMatrix(backend,FTB,CG,Grid,nQuadM,FEMSei.Jacobi!)
CG.LUM = lu(CG.M)
VecDG.M = FEMSei.MassMatrix(backend,FTB,VecDG,Grid,nQuadM,FEMSei.Jacobi!)
VecDG.LUM = lu(VecDG.M)
RT.M = FEMSei.MassMatrix(backend,FTB,RT,Grid,nQuadM,FEMSei.Jacobi!)
RT.LUM = lu(RT.M)
ND.M = FEMSei.MassMatrix(backend,FTB,ND,Grid,nQuadM,FEMSei.Jacobi!)
ND.LUM = lu(ND.M)
Curl = FEMSei.CurlMatrix(backend,FTB,ND,DG,Grid,nQuad,FEMSei.Jacobi!)

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
cName = ["h";"Vort";"uS";"vS"]

FEMSei.InterpolatehRT!(Uhu,RT,FEMSei.Jacobi!,Grid,Grid.Type,nQuad,Model.InitialProfile)
FEMSei.InterpolateDG!(Uh,DG,FEMSei.Jacobi!,Grid,Grid.Type,Model.InitialProfile)
# Output of the initial values
FileNumber=0
NumRefine = size(vtkSkeletonMesh.RefineMidPoints,1)
VelSp = zeros(Grid.NumFaces*NumRefine,2)
hout = zeros(Grid.NumFaces*NumRefine)
Vort = zeros(Grid.NumFaces*NumRefine)
FEMSei.ConvertScalar!(backend,FTB,hout,Uh,DG,Grid,FEMSei.Jacobi!,vtkSkeletonMesh.RefineMidPoints)
FEMSei.ConvertScalarVelocitySp!(backend,FTB,VelSp,Uhu,RT,Uh,DG,Grid,FEMSei.Jacobi!,vtkSkeletonMesh.RefineMidPoints)
FEMSei.Vorticity!(backend,FTB,Vort,CG,Uhu,RT,Uh,DG,Grid,Grid.Type,nQuad,FEMSei.Jacobi!,
  vtkSkeletonMesh.RefineMidPoints)
Outputs.vtkSkeleton!(vtkSkeletonMesh, FileNameOutput, Proc, ProcNumber, [hout Vort VelSp] ,FileNumber,cName)

for i = 1 : nAdveVel
  @show i,(i-1)*dtau/3600 
  @. F = 0  
  # Tendency h
  FEMSei.DivRhs!(backend,FTB,Fh,DG,Uhu,RT,Grid,Grid.Type,nQuad,FEMSei.Jacobi!)
  ldiv!(DG.LUM,Fh)
  # Tendency hu
  FEMSei.InterpolateScalarHDivVecDG!(backend,FTB,uRec,VecDG,Uh,DG,Uhu,RT,Grid,
    Grid.Type,nQuad,FEMSei.Jacobi!)

  FEMSei.DivMomentumVector!(backend,FTB,Fhu,RT,Uhu,RT,uRec,VecDG,Grid,Grid.Type,nQuad,FEMSei.Jacobi!)
  FEMSei.CrossRhs!(backend,FTB,Fhu,RT,Uhu,RT,Grid,Grid.Type,nQuad,FEMSei.Jacobi!)
  FEMSei.GradHeightSquared!(backend,FTB,Fhu,RT,Uh,DG,Grid,Grid.Type,nQuad,FEMSei.Jacobi!)
  ldiv!(RT.LUM,Fhu)
  
  @. UNew = U + 1 / 3 * dtau * F

  @. F = 0  
  # Tendency h
  FEMSei.DivRhs!(backend,FTB,Fh,DG,UNewhu,RT,Grid,Grid.Type,nQuad,FEMSei.Jacobi!)
  ldiv!(DG.LUM,Fh)
  # Tendency hu
  FEMSei.InterpolateScalarHDivVecDG!(backend,FTB,uRec,VecDG,UNewh,DG,UNewhu,RT,Grid,
    Grid.Type,nQuad,FEMSei.Jacobi!)
  FEMSei.DivMomentumVector!(backend,FTB,Fhu,RT,UNewhu,RT,uRec,VecDG,Grid,Grid.Type,nQuad,FEMSei.Jacobi!)
  FEMSei.CrossRhs!(backend,FTB,Fhu,RT,UNewhu,RT,Grid,Grid.Type,nQuad,FEMSei.Jacobi!)
  FEMSei.GradHeightSquared!(backend,FTB,Fhu,RT,UNewh,DG,Grid,Grid.Type,nQuad,FEMSei.Jacobi!)
  ldiv!(RT.LUM,Fhu)
  @. UNew = U + 0.5 * dtau * F

  @. F = 0  
  # Tendency h
  FEMSei.DivRhs!(backend,FTB,Fh,DG,UNewhu,RT,Grid,Grid.Type,nQuad,FEMSei.Jacobi!)
  ldiv!(DG.LUM,Fh)
  # Tendency hu
  FEMSei.InterpolateScalarHDivVecDG!(backend,FTB,uRec,VecDG,UNewh,DG,UNewhu,RT,Grid,
    Grid.Type,nQuad,FEMSei.Jacobi!)
  FEMSei.DivMomentumVector!(backend,FTB,Fhu,RT,UNewhu,RT,uRec,VecDG,Grid,Grid.Type,nQuad,FEMSei.Jacobi!)
  FEMSei.CrossRhs!(backend,FTB,Fhu,RT,UNewhu,RT,Grid,Grid.Type,nQuad,FEMSei.Jacobi!)
  FEMSei.GradHeightSquared!(backend,FTB,Fhu,RT,UNewh,DG,Grid,Grid.Type,nQuad,FEMSei.Jacobi!)
  ldiv!(RT.LUM,Fhu)
  @. U = U + dtau * F
    
  # Output
  if mod(i,nprint) == 0 
    global FileNumber += 1
    FEMSei.ConvertScalar!(backend,FTB,hout,Uh,DG,Grid,FEMSei.Jacobi!,vtkSkeletonMesh.RefineMidPoints)
    FEMSei.ConvertScalarVelocitySp!(backend,FTB,VelSp,Uhu,RT,Uh,DG,Grid,FEMSei.Jacobi!,vtkSkeletonMesh.RefineMidPoints)
    FEMSei.Vorticity!(backend,FTB,Vort,CG,Uhu,RT,Uh,DG,Grid,Grid.Type,nQuad,FEMSei.Jacobi!,
      vtkSkeletonMesh.RefineMidPoints)
    Outputs.vtkSkeleton!(vtkSkeletonMesh, FileNameOutput, Proc, ProcNumber, [hout Vort VelSp] ,FileNumber,cName)
  end
end
@show "finished"

