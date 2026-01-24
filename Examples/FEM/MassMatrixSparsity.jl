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
using SparseArrays

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

RefineLevel = 5
nz = 1
nPanel = 50
nQuad = 3
Decomp = ""
Decomp = "EqualArea"
Problem = "GalewskySphere"
RadEarth = Phys.RadEarth
Problem = "LinearBlob"
Param = Examples.Parameters(FTB,Problem)
Examples.InitialProfile!(backend,FTB,Model,Problem,Param,Phys)
GridType = "TriangularSphere"

#TRI
Grid, Exchange = Grids.InitGridSphere(backend,FTB,OrdPoly,nz,nPanel,RefineLevel,ns,
  nLat,nLon,LatB,GridType,Decomp,RadEarth,Model,ParallelCom;ChangeOrient=2)

KiteGrid = Grids.Grid2KiteGrid(backend,FTB,Grid,Grids.OrientFaceSphere)


Jacobi = FEM.Jacobi!  

# Quadrature rules
if Grid.Type == Grids.Quad()
  nQuad = 3
  nQuadM = 3
  nQuadS = 3
elseif Grid.Type == Grids.Tri()
  nQuad = 4
  nQuadM = 4
  nQuadS = 4
end

# Finite elements
k = 1
@show k
RT = FEM.RTKiteDualHDiv{FTB}(k,Grids.Quad(),backend,KiteGrid)
DG = FEM.CGKitePrimalStruct{FTB}(k+1,Grids.Quad(),backend,KiteGrid)
CG = FEM.CGStruct{FTB}(backend,k+1,Grids.Quad(),KiteGrid)
ND = FEM.CG1KiteDualHCurl{FTB}(Grids.Quad(),backend,KiteGrid)


ModelFEM = FEM.ModelFEMVecI(backend,FTB,ND,RT,CG,DG,KiteGrid,nQuadM,nQuadS,Jacobi)

global MDG = Array{AbstractSparseMatrix,1}(undef,KiteGrid.NumFaces)
global iM = 1
for iN = 1 : KiteGrid.NumNodes
  Index = []  
  if KiteGrid.Nodes[iN].Type == 'F'
    for iF in KiteGrid.Nodes[iN].F
      for iDoF = 1 : DG.DoF  
        push!(Index,DG.Glob[iDoF,iF])  
      end  
    end
    unique!(Index)
    global MDG[iM] = DG.M[Index,Index]
    global iM += 1
  end    
end

global NumNodesN = 0
for iN = 1 : KiteGrid.NumNodes
  if KiteGrid.Nodes[iN].Type == 'N'
    global NumNodesN += 1
  end
end  

global MRT = Array{AbstractSparseMatrix,1}(undef, NumNodesN)
global IndexGlob = Array{Array{Int, 1}, 1}(undef, NumNodesN)
global iM = 1
for iN = 1 : KiteGrid.NumNodes
  Index = []
  if KiteGrid.Nodes[iN].Type == 'N'
    for iF in KiteGrid.Nodes[iN].F
      for iDoF = 1 : RT.DoF
        push!(Index,RT.Glob[iDoF,iF])
      end
    end
    unique!(Index)
    global MRT[iM] = RT.M[Index,Index]
    global IndexGlob[iM] = Index
    global iM += 1
  end
end

nothing









