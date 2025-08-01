import CGDycore:
  Examples, Parallels, Models, Grids, Outputs, Integration,  GPUS, DyCore, Surfaces
using MPI
using Base
using CUDA
using AMDGPU
using Metal
using KernelAbstractions
using StaticArrays
using ArgParse

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

# Physical parameters
Phys = DyCore.PhysParameters{FTB}()

z0M = 0.05
z0H = 0.05
RiB = 0.4
LandClass = 1
TSurf = 280.0
T1 = 280.0
VT = 2.0
uf = Surfaces.Businger{FTB}()
theta = FTB(300)
thetaS = FTB(301.0)
RhoPos = 1
uPos = 2
vPos = 3
wPos = 4
ThPos = 5
NumG = 10
nz = 2
U = zeros(FTB,nz,NumG,5)
p = zeros(FTB,nz,NumG)
@. U[:,:,1] = 1.0
@. U[:,:,2] = 2.0
@. U[:,:,3] = 2.0
@. U[:,:,5] = theta
dXdxI = zeros(FTB,3)
dXdxI[3] = 1
nS = zeros(FTB,3,NumG)
dz = zeros(nz,NumG)
@. dz = 10.0
@. nS[3,:] = 1
@. p = 1
#thetaS = FTB(295)

SurfaceData = Surfaces.SurfaceData{FTB}(backend,Surfaces.LenSurfaceData,NumG)
LandUseData = Surfaces.LandUseData{FTB}(backend,NumG)
@. LandUseData.LandClass[1:5] = 1
@. LandUseData.LandClass[6:end] = 2

for iG = 1 : NumG
  a = rand(1)  
  SurfaceData.Data[Surfaces.TSurfPos,iG] = 300.0 + a[1]
end  

SurfaceFluxValues = Surfaces.MOSurfaceFlux()(Surfaces.Businger(),Phys,RhoPos,uPos,
      vPos,wPos,ThPos,LandUseData)


LandClass = LandUseData.LandClass
NumG = size(U,2)
groupS = (1)
ndrangeS = (NumG)
KSurfaceFluxDataKernel! = Surfaces.SurfaceFluxDataKernel!(backend,groupS)
KSurfaceFluxDataKernel!(SurfaceFluxValues,SurfaceData.Data,U,p,dz,nS,
    LandClass,ndrange=ndrangeS)
@. U[:,:,5] = theta + 0.05
@show "Second "
KSurfaceFluxDataKernel!(SurfaceFluxValues,SurfaceData.Data,U,p,dz,nS,
    LandClass,ndrange=ndrangeS)
stop

#SurfaceData = Surfaces.MOSurface()(uf,Phys,RhoPos,uPos,vPos,wPos,ThPos)
#uStar, CT, CH = SurfaceData(dz,U,p,dXdxI,nS,landuse)
#@show uStar, CT, CH
#stop    

@show dz, VT, theta, thetaS
CM, CT, uStar = Surfaces.MOSTIteration(uf,z0M,z0H,dz,VT,theta,thetaS,LandClass,Phys)
@show CM, CT, uStar
VT = 20.0
theta = FTB(300)
thetaS = FTB(300)
z0M = 0.00001
z0H = 0.00001
@show dz, VT, theta, thetaS
CM, CT, uStar = Surfaces.MOSTSeaIteration(uf,z0M,z0H,dz,VT,theta,thetaS,LandClass,Phys)
@show CM, CT, uStar
theta = FTB(301)
thetaS = FTB(300)
@show dz, VT, theta, thetaS
CM, CT, uStar = Surfaces.MOSTSeaIteration(uf,z0M,z0H,dz,VT,theta,thetaS,LandClass,Phys)
@show CM, CT, uStar
theta = FTB(300)
thetaS = FTB(301)
@show dz, VT, theta, thetaS
CM, CT, uStar = Surfaces.MOSTSeaIteration(uf,z0M,z0H,dz,VT,theta,thetaS,LandClass,Phys)
@show CM, CT, uStar
stop


RiB  = Phys.Grav/TSurf*((T1-TSurf)*(dz-z0M))/(VT*VT)
@show RiB
Cm, Ch = Surfaces.RichardIteration(z0M,z0h,dz,RiB,LandClass)
@show Cm, Ch




RiB  = Phys.Grav/TSurf*((T1-TSurf)*(dz-z0M))/(VT*VT)
Cm, Ch = Surfaces.RichardIteration(z0M,z0h,dz,RiB,LandClass)
@show Cm, Ch


