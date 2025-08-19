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

Problem = "GalewskySphere"
RadEarth = Phys.RadEarth
dtau = 6
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
Examples.InitialProfile!(Model,Problem,Param,Phys)

RefineLevel = 5
nz = 1
nPanel = 60
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
  Model,ParallelCom;order=false)

#=
for iE = 1 : Grid.NumEdges
  iF1 = Grid.Edges[iE].F[1]
  for j = 1 : length(Grid.Faces[iF1].E)
    jE = Grid.Faces[iF1].E[j]  
    if jE == iE
      if Grid.Faces[iF1].OrientE[j] == -1
        Grid.Edges[iE].F[1] = Grid.Edges[iE].F[2]
        Grid.Edges[iE].F[2] = iF1
      end  
    end  
  end  
end  

@show Grid.NumFaces
=#


Grids.TestGrid(Grid)
vtkSkeletonMesh = Outputs.vtkStruct{Float64}(backend,Grid,Grid.NumFaces,Flat)
vtkSkeletonMeshGhost = Outputs.vtkStruct{Float64}(backend,Grid,Grid.NumFaces+Grid.NumFacesG,Flat)
#KiteFaces = FiniteVolumes.MatrixTangential(FiniteVolumes.JacobiSphere,Grid)

MetricFV = FiniteVolumes.MetricFiniteVolume(backend,FTB,Grid)

#Div = FiniteVolumes.Divergence(backend,FTB,MetricFV,Grid)
#GradFV = FiniteVolumes.Gradient(backend,FTB,MetricFV,Grid)
Grad,Inter = FiniteVolumes.GradMPFATri(backend,FTB,Grid)
GradQuad,InterQuad = FiniteVolumes.GradMPFA(backend,FTB,Grid)
Div = FiniteVolumes.DivMPFA(backend,FTB,MetricFV,Grid)


pPosS = 1
pPosEI = Grid.NumFaces
pPosE = Grid.NumFaces + Grid.NumFacesG
uPosS = pPosE + 1
uPosE = pPosE + Grid.NumEdges
U = zeros(FTB,Grid.NumFaces+Grid.NumFacesG+Grid.NumEdges)
r = similar(U)
UNew = similar(U)
@views rp = r[pPosS:pPosE]
@views rpI = r[pPosS:pPosEI]
@views ru = r[uPosS:uPosE]
@views Up = U[pPosS:pPosE]
@views UpI = U[pPosS:pPosEI]
@views Uu = U[uPosS:uPosE]
@views UNewp = UNew[pPosS:pPosE]
@views UNewpI = UNew[pPosS:pPosEI]
@views UNewu = UNew[uPosS:uPosE]
FiniteVolumes.ProjectFace!(backend,FTB,UpI,Grid,Model.InitialProfile)
FiniteVolumes.ProjectEdge!(backend,FTB,Uu,Grid,Model.InitialProfile)
FileNumber = 0
VelCa = zeros(Grid.NumFaces,Grid.Dim)
VelSp = zeros(Grid.NumFaces,2)
FiniteVolumes.ConvertVelocityCart!(backend,FTB,VelCa,Uu,Grid)
FiniteVolumes.ConvertVelocitySp!(backend,FTB,VelSp,Uu,Grid)
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType*"FV", Proc, ProcNumber, [UpI VelCa VelSp], FileNumber)
uCurlN = zeros(FTB,Grid.NumNodes)
uTan = zeros(FTB,Grid.NumEdges)
K = zeros(FTB,Grid.NumFaces)
Curl = FiniteVolumes.CurlNodeMatrix(MetricFV,Grid)
Tang = FiniteVolumes.TagentialVelocityMatrix(MetricFV,Grid)
mul!(uCurlN,Curl,Uu)
mul!(uTan,Tang,Uu)
FiniteVolumes.TagentialVelocity2(uTan,Uu,MetricFV,Grid)
FiniteVolumes.KineticEnergy(K,Uu,MetricFV,Grid)
FiniteVolumes.KineticEnergy(K,Uu,uTan,MetricFV,Grid)
FiniteVolumes.ConvertVelocityTCart!(backend,FTB,VelCa,uTan,Grid)
FiniteVolumes.ConvertVelocityTSp!(backend,FTB,VelSp,uTan,Grid)

uCurl = zeros(FTB,Grid.NumFaces)
DualVolume = zeros(FTB,Grid.NumFaces)
for iF = 1 : Grid.NumFaces
   for iN in Grid.Faces[iF].N
     uCurl[iF] += uCurlN[iN] * MetricFV.DualVolume[iN]
     DualVolume[iF] += MetricFV.DualVolume[iN]
   end
end   
@. uCurl /= DualVolume

@show GridType*"FVCurl"*"$RefineLevel"
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType*"FVCurl", Proc, ProcNumber, [uCurl K VelSp], 0)

hE = zeros(FTB,Grid.NumEdges)
mul!(hE,Grad,Up)
FiniteVolumes.ConvertVelocitySp!(backend,FTB,VelSp,hE,Grid)
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType*"FVCurl", Proc, ProcNumber, [uCurl K VelSp], 1)


hEExact = zeros(FTB,Grid.NumEdges)
hR = zeros(FTB,Grid.NumFaces)
mul!(hE,Inter,Up)
FiniteVolumes.ProjectEdgeScalar!(backend,FTB,hEExact,Grid,Model.InitialProfile)
for iF = 1 : Grid.NumFaces
   for iE in Grid.Faces[iF].E
     if Grid.Edges[iE].F[1] == iF  
       hR[iF] += (hE[iE] - hEExact[iE]) * MetricFV.DualEdgeVolume[1,iE] / MetricFV.PrimalVolume[iF]
     else  
       hR[iF] += (hE[iE] - hEExact[iE]) * MetricFV.DualEdgeVolume[2,iE] / MetricFV.PrimalVolume[iF]
     end  
   end
end
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType*"FVCurl", Proc, ProcNumber, [hR K VelSp], 2)
stop
#@. Uu = Uu * hE
mul!(hR,Div,Uu)
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType*"FVCurl", Proc, ProcNumber, [hR K VelSp], 3)
stop
pNeu = zeros(FTB,Grid.NumFaces)
stop

time = 0.0
Parallels.ExchangeDataFSendGPU(reshape(Up,1,length(Up),1),Exchange)
Parallels.ExchangeDataFRecvGPU!(reshape(Up,1,length(Up),1),Exchange)
Outputs.vtkSkeleton!(vtkSkeletonMeshGhost, GridType*"FVG", Proc, ProcNumber, [Up Up] , FileNumber)


@show dtau,nAdveVel
for i = 1 : nAdveVel
  if Proc == 1  
    @show i  
  end  
  Parallels.ExchangeDataFSendGPU(reshape(Up,1,length(Up),1),Exchange)
  Parallels.ExchangeDataFRecvGPU!(reshape(Up,1,length(Up),1),Exchange)
  mul!(rpI,Div,Uu)
  mul!(ru,Grad,Up)
  @. ru = -ru
  @. UNew = U + 0.5 * dtau * r  

  Parallels.ExchangeDataFSendGPU(reshape(UNewp,1,length(Up),1),Exchange)
  Parallels.ExchangeDataFRecvGPU!(reshape(UNewp,1,length(Up),1),Exchange)
  mul!(rpI,Div,UNewu)
  mul!(ru,Grad,UNewp)
  @. ru = -ru
  @. U = U + dtau * r  
  
  if mod(i,PrintStp) == 0
    global FileNumber += 1
    FiniteVolumes.ConvertVelocityCart!(backend,FTB,VelCa,Uu,Grid)
    FiniteVolumes.ConvertVelocitySp!(backend,FTB,VelSp,Uu,Grid)
    Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType*"FV", Proc, ProcNumber, [UpI VelCa VelSp], FileNumber)
    MPI.Barrier(comm)
  end  

  #time = time + dtau
end


FileNumber += 1
FiniteVolumes.ConvertVelocityCart!(backend,FTB,VelCa,Uu,Grid)
FiniteVolumes.ConvertVelocitySp!(backend,FTB,VelSp,Uu,Grid)
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType*"FV", Proc, ProcNumber, [UpI VelCa VelSp], FileNumber)


#= ToDO
Gleichungen testen mit @show
mehr Quadraturregeln für Dreiecke höhere ordnungen
RT1, DG1
Wiki von github
=#
