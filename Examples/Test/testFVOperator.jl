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

function TestOperators()
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
RadEarth = 300.0
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

RefineLevel = 6
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

Grids.TestGrid(Grid)
vtkSkeletonMesh = Outputs.vtkStruct{Float64}(backend,Grid,Grid.NumFaces,Flat)
vtkSkeletonMeshGhost = Outputs.vtkStruct{Float64}(backend,Grid,Grid.NumFaces+Grid.NumFacesG,Flat)

MetricFV = FiniteVolumes.MetricFiniteVolume(backend,FTB,Grid)

Grad,Inter = FiniteVolumes.GradMPFATri(backend,FTB,Grid)
DivMatrix = FiniteVolumes.DivMPFA(backend,FTB,MetricFV,Grid)
CurlMatrix = FiniteVolumes.CurlNodeMatrix(MetricFV,Grid)
Tang = FiniteVolumes.TagentialVelocityMatrix(MetricFV,Grid)


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

h = zeros(FTB,Grid.NumFaces)
hEx = zeros(FTB,Grid.NumFaces)
hEdge = zeros(FTB,Grid.NumEdges)
hEdgeEx = zeros(FTB,Grid.NumEdges)
Div = zeros(FTB,Grid.NumFaces)
DivEx = zeros(FTB,Grid.NumFaces)
CurlN = zeros(FTB,Grid.NumNodes)
CurlNEx = zeros(FTB,Grid.NumNodes)
Curl = zeros(FTB,Grid.NumFaces)
CurlEx = zeros(FTB,Grid.NumFaces)
uN = zeros(FTB,Grid.NumEdges)
uT = zeros(FTB,Grid.NumEdges)
uTEx = zeros(FTB,Grid.NumEdges)
hGrad = zeros(FTB,Grid.NumEdges)
hGradEx = zeros(FTB,Grid.NumEdges)
VelSp = zeros(Grid.NumFaces,2)
VelSpE = zeros(Grid.NumFaces,2)
VelCart = zeros(3)
VelSphere = zeros(3)

# Test gradient
for iF = 1 : Grid.NumFaces
  x = Grid.Faces[iF].Mid.x  
  y = Grid.Faces[iF].Mid.y  
  z = Grid.Faces[iF].Mid.z  
  h[iF] = 3*x^4 + 4*y^3 + 5 * z^2  
end  
for iE = 1 : Grid.NumEdges
  x = Grid.Edges[iE].Mid.x  
  y = Grid.Edges[iE].Mid.y  
  z = Grid.Edges[iE].Mid.z  
# h[iF] = 3*x^4 + 4*y^3 + 5 * z^2  
  hGradx = 12 * x^3
  hGrady = 12 * y^2
  hGradz = 10 * z
  hGradEx[iE] = hGradx * Grid.Edges[iE].n.x +
               hGrady * Grid.Edges[iE].n.y
               hGradz * Grid.Edges[iE].n.z
end               
mul!(hGrad,Grad,h)
FiniteVolumes.ConvertVelocitySp!(backend,FTB,VelSp,hGrad,Grid)
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType*"FVOperator", Proc, ProcNumber, [h VelSp], 0)
FiniteVolumes.ConvertVelocitySp!(backend,FTB,VelSpE,hGradEx,Grid)
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType*"FVOperator", Proc, ProcNumber, [h VelSpE], 1)
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType*"FVOperator", Proc, ProcNumber, [h VelSp-VelSpE], 2)

# Test tangential
for iE = 1 : Grid.NumEdges
  x = Grid.Edges[iE].Mid.x
  y = Grid.Edges[iE].Mid.y
  z = Grid.Edges[iE].Mid.z
  lon,lat,_= Grids.cart2sphere(x,y,z)
  VelSphere[1] = sin(lon)^4 * cos(lat)^3 
  VelSphere[2] = sin(lon)^3 * cos(lat)^2
  VelCart = FiniteVolumes.VelSphere2Cart(VelSphere,lon,lat)
  n1 = Grid.Edges[iE].n.x
  n2 = Grid.Edges[iE].n.y
  n3 = Grid.Edges[iE].n.z
  uN[iE] = n1 * VelCart[1] + n2 * VelCart[2] + n3 * VelCart[3]
  t1 = Grid.Edges[iE].t.x
  t2 = Grid.Edges[iE].t.y
  t3 = Grid.Edges[iE].t.z
  uTEx[iE] = t1 * VelCart[1] + t2 * VelCart[2] + t3 * VelCart[3]
end
mul!(uT,Tang,uN)
FiniteVolumes.ConvertVelocityTSp!(backend,FTB,VelSpE,uTEx,Grid)
FiniteVolumes.ConvertVelocityTSp!(backend,FTB,VelSp,uT,Grid)
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType*"FVOperator", Proc, ProcNumber, [h VelSp], 3)
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType*"FVOperator", Proc, ProcNumber, [h VelSpE], 4)
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType*"FVOperator", Proc, ProcNumber, [h VelSp-VelSpE], 5)

# Test Divergence 
for iF = 1 : Grid.NumFaces
  x = Grid.Faces[iF].Mid.x
  y = Grid.Faces[iF].Mid.y
  z = Grid.Faces[iF].Mid.z
  lon,lat,_= Grids.cart2sphere(x,y,z)
  if abs(cos(lat)) == 0
    if lat > 0
      lat -= 1.e-4
    else  
      lat += 1.e-4
    end
  end  
  uS = sin(lon)^4 * cos(lat)^3 
  vS = sin(lon)^3 * cos(lat)^2
  duSdlon = 4 * sin(lon)^3 * cos(lon) * cos(lat)^3
  duSdlat = -3 * sin(lon)^4 * cos(lat)^2 * sin(lat)
  dvSdlon = 3 * sin(lon)^2 * cos(lon) * cos(lat)^2
  dvSdlat = -2 * sin(lon)^3 * cos(lat) * sin(lat)
  DivEx[iF] = -1/cos(lat) * (-sin(lat)*vS + cos(lat) * dvSdlat + duSdlon) / RadEarth
end  
mul!(Div,DivMatrix,uN)
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType*"FVOperator", Proc, ProcNumber, [DivEx VelSp], 6)
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType*"FVOperator", Proc, ProcNumber, [Div VelSp], 7)
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType*"FVOperator", Proc, ProcNumber, [DivEx-Div VelSp], 8)

# Test Rotation
for iN = 1 : Grid.NumNodes
  x = Grid.Nodes[iN].P.x
  y = Grid.Nodes[iN].P.y
  z = Grid.Nodes[iN].P.z
  lon,lat,_= Grids.cart2sphere(x,y,z)
  if abs(cos(lat)) == 0
    if lat > 0
      lat -= 1.e-4
    else
      lat += 1.e-4
    end
  end
  uS = sin(lon)^4 * cos(lat)^3
  vS = sin(lon)^3 * cos(lat)^2
  duSdlon = 4 * sin(lon)^3 * cos(lon) * cos(lat)^3
  duSdlat = -3 * sin(lon)^4 * cos(lat)^2 * sin(lat)
  dvSdlon = 3 * sin(lon)^2 * cos(lon) * cos(lat)^2
  dvSdlat = -2 * sin(lon)^3 * cos(lat) * sin(lat)
  CurlNEx[iN] = -1/cos(lat) * (-sin(lat) * uS + cos(lat) * duSdlat - dvSdlon) / RadEarth
end
mul!(CurlN,CurlMatrix,uN)
DualVolume = zeros(FTB,Grid.NumFaces)
for iF = 1 : Grid.NumFaces
   for iN in Grid.Faces[iF].N
     Curl[iF] += CurlN[iN] * MetricFV.DualVolume[iN]
     CurlEx[iF] += CurlNEx[iN] * MetricFV.DualVolume[iN]
     DualVolume[iF] += MetricFV.DualVolume[iN]
   end
end
@. Curl /= DualVolume
@. CurlEx /= DualVolume
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType*"FVOperator", Proc, ProcNumber, [CurlEx VelSp], 9)
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType*"FVOperator", Proc, ProcNumber, [Curl VelSp], 10)
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType*"FVOperator", Proc, ProcNumber, [CurlEx-Curl VelSp], 11)

# Test interpolation
for iF = 1 : Grid.NumFaces
  x = Grid.Faces[iF].Mid.x
  y = Grid.Faces[iF].Mid.y
  z = Grid.Faces[iF].Mid.z
  hEx[iF] = 3*x^4 + 4*y^3 + 5 * z^2
end
for iE = 1 : Grid.NumEdges
  x = Grid.Edges[iE].Mid.x
  y = Grid.Edges[iE].Mid.y
  z = Grid.Edges[iE].Mid.z
  hEdgeEx[iE] = 3*x^4 + 4*y^3 + 5 * z^2
end
mul!(hEdge,Inter,hEx)
@. h = 0
for iE = 1 : Grid.NumEdges
  iF1 = Grid.Edges[iE].F[1]  
  iF2 = Grid.Edges[iE].F[2]  
  h[iF1] += hEdge[iE] * MetricFV.DualEdgeVolume[1,iE]
  h[iF2] += hEdge[iE] * MetricFV.DualEdgeVolume[2,iE]
end  
@. h /= MetricFV.PrimalVolume


Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType*"FVOperator", Proc, ProcNumber, [hEx VelSp], 12)
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType*"FVOperator", Proc, ProcNumber, [h VelSp], 13)
Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType*"FVOperator", Proc, ProcNumber, [h - hEx VelSp], 14)
return 
end



# theta lat 
# phi   lon
# Div = 1/cos(theta)*d/dtheta(cos(theta)*u_theta) + 1/cos(theta)*d/dphi(u_phi)
# Curl = 1/cos(theta)*(d/dtheta(cos(theta)*u_phi) - d/dphi(u_theta))



