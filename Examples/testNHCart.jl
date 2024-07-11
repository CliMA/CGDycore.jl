import CGDycore:
  Thermodynamics, Examples, Parallels, Models, Grids, Surfaces,  Outputs, Integration,  GPU, DyCore
using MPI
using Base
using CUDA
using AMDGPU
using Metal
using KernelAbstractions
using StaticArrays
using ArgParse
using MPI


# Model
parsed_args = DyCore.parse_commandline()
Problem = parsed_args["Problem"]
ProfRho = parsed_args["ProfRho"]
ProfTheta = parsed_args["ProfTheta"]
PertTh = parsed_args["PertTh"]
ProfVel = parsed_args["ProfVel"]
ProfpBGrd = parsed_args["ProfpBGrd"]
ProfRhoBGrd = parsed_args["ProfRhoBGrd"]
RhoVPos = parsed_args["RhoVPos"]
RhoCPos = parsed_args["RhoCPos"]
RhoIPos = parsed_args["RhoIPos"]
RhoRPos = parsed_args["RhoRPos"]
TkePos = parsed_args["TkePos"]
HorLimit = parsed_args["HorLimit"]
Upwind = parsed_args["Upwind"]
Damping = parsed_args["Damping"]
Relax = parsed_args["Relax"]
StrideDamp = parsed_args["StrideDamp"]
NumV = parsed_args["NumV"]
NumTr = parsed_args["NumTr"]
BoundaryWE = parsed_args["BoundaryWE"]
BoundarySN = parsed_args["BoundarySN"]
BoundaryBT = parsed_args["BoundaryBT"]
Thermo = parsed_args["Thermo"]
State = parsed_args["State"]
RefProfile = parsed_args["RefProfile"]
Profile = parsed_args["Profile"]
Curl = parsed_args["Curl"]
ModelType = parsed_args["ModelType"]
Equation = parsed_args["Equation"]
Microphysics = parsed_args["Microphysics"]
TypeMicrophysics = parsed_args["TypeMicrophysics"]
RelCloud = parsed_args["RelCloud"]
Rain = parsed_args["Rain"]
#Orography
TopoS = parsed_args["TopoS"]
P1 = parsed_args["P1"]
P2 = parsed_args["P2"]
P3 = parsed_args["P3"]
P4 = parsed_args["P4"]

# Parallel
Decomp = parsed_args["Decomp"]
SimDays = parsed_args["SimDays"]
SimHours = parsed_args["SimHours"]
SimMinutes = parsed_args["SimMinutes"]
SimSeconds = parsed_args["SimSeconds"]
SimTime = parsed_args["SimTime"]
dtau = parsed_args["dtau"]
IntMethod = parsed_args["IntMethod"]
Table = parsed_args["Table"]
GridType = parsed_args["GridType"]
AdaptGridType = parsed_args["AdaptGridType"]
Coriolis = parsed_args["Coriolis"]
CoriolisType = parsed_args["CoriolisType"]
Buoyancy = parsed_args["Buoyancy"]
Turbulence = parsed_args["Turbulence"]
Source = parsed_args["Source"]
VerticalDiffusion = parsed_args["VerticalDiffusion"]
JacVerticalDiffusion = parsed_args["JacVerticalDiffusion"]
JacVerticalAdvection = parsed_args["JacVerticalAdvection"]
VerticalDiffusionMom = parsed_args["VerticalDiffusionMom"]
SurfaceFlux = parsed_args["SurfaceFlux"]
SurfaceFluxMom = parsed_args["SurfaceFluxMom"]
SurfaceScheme = parsed_args["SurfaceScheme"]
# Grid
nx = parsed_args["nx"]
ny = parsed_args["ny"]
nz = parsed_args["nz"]
H = parsed_args["H"]
Stretch = parsed_args["Stretch"]
StretchType = parsed_args["StretchType"]
OrdPoly = parsed_args["OrdPoly"]
Lx = parsed_args["Lx"]
Ly = parsed_args["Ly"]
x0 = parsed_args["x0"]
y0 = parsed_args["y0"]
# Viscosity
HyperVisc = parsed_args["HyperVisc"]
HyperDCurl = parsed_args["HyperDCurl"]
HyperDGrad = parsed_args["HyperDGrad"]
HyperDDiv = parsed_args["HyperDDiv"]
# Output
OrdPrint = parsed_args["OrdPrint"]
PrintDays = parsed_args["PrintDays"]
PrintHours = parsed_args["PrintHours"]
PrintMinutes = parsed_args["PrintMinutes"]
PrintSeconds = parsed_args["PrintSeconds"]
PrintTime = parsed_args["PrintTime"]
PrintStartTime = parsed_args["PrintStartTime"]
# Device
Device = parsed_args["Device"]
GPUType = parsed_args["GPUType"]
FloatTypeBackend = parsed_args["FloatTypeBackend"]
NumberThreadGPU = parsed_args["NumberThreadGPU"]

MPI.Init()
comm = MPI.COMM_WORLD
Proc = MPI.Comm_rank(comm) + 1
ProcNumber = MPI.Comm_size(comm)
ParallelCom = DyCore.ParallelComStruct()
ParallelCom.Proc = Proc
ParallelCom.ProcNumber  = ProcNumber

JuliaDevice = get(ENV, "JuliaDevice", "CPU")
JuliaGPU = get(ENV, "JuliaGPU", "CUDA")

if JuliaDevice == "CPU" 
  backend = CPU()
elseif JuliaDevice == "GPU" 
  if JuliaGPU == "CUDA"
    backend = CUDABackend()
    CUDA.allowscalar(false)
    CUDA.device!(MPI.Comm_rank(MPI.COMM_WORLD))
  elseif JuliaGPU == "AMD"
    backend = ROCBackend()
    AMDGPU.allowscalar(false)
  elseif JuliaGPU == "Metal"
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
Param = Examples.Parameters(FTB,Problem)

KernelAbstractions.synchronize(backend)

OrdPolyZ=1
Parallel = true

# Physical parameters
Phys = DyCore.PhysParameters{FTB}()

#ModelParameters
Model = DyCore.ModelStruct{FTB}()

# Initial conditions
Model.NumV=NumV
Model.NumTr=NumTr
Model.Problem=Problem
if ProfRho == ""
  Model.ProfRho = Problem
else
  Model.ProfRho = ProfRho  
end  
if ProfTheta == ""
  Model.ProfTheta = Problem
else
  Model.ProfTheta = ProfTheta  
end  
if ProfVel == ""
  Model.ProfVel = Problem
else
  Model.ProfVel = ProfVel  
end  
Model.ProfpBGrd = ProfpBGrd
Model.ProfRhoBGrd = ProfRhoBGrd
Model.RefProfile = RefProfile
Model.Profile = Profile
Model.RhoPos=1
Model.uPos=2
Model.vPos=3
Model.wPos=4
Model.ThPos=5
Model.RhoVPos  = RhoVPos
Model.RhoCPos  = RhoCPos
Model.RhoIPos  = RhoIPos
Model.RhoRPos  = RhoRPos
Model.TkePos  = TkePos
Model.HorLimit = HorLimit
Model.Upwind = Upwind
Model.Damping = Damping
Model.StrideDamp = StrideDamp
Model.Relax = Relax
Model.Coriolis = Coriolis
Model.CoriolisType = CoriolisType
Model.Buoyancy = Buoyancy
Model.Turbulence = Turbulence
Model.VerticalDiffusion = VerticalDiffusion
Model.JacVerticalDiffusion = JacVerticalDiffusion
Model.JacVerticalAdvection = JacVerticalAdvection
Model.VerticalDiffusionMom = VerticalDiffusionMom
Model.Source = Source
Model.Microphysics = Microphysics
Model.RelCloud = RelCloud
Model.Rain = Rain
Model.Source = Source
Model.SurfaceFlux = SurfaceFlux
Model.SurfaceFluxMom = SurfaceFluxMom
Model.Thermo = Thermo
Model.Curl = Curl
Model.Stretch = Stretch
Model.StretchType = StretchType
Model.ModelType = ModelType
Model.HyperVisc = HyperVisc
Model.HyperDCurl = HyperDCurl # =7.e15
Model.HyperDGrad = HyperDGrad # =7.e15
Model.HyperDDiv = HyperDDiv # =7.e15




Boundary = Grids.Boundary()
Boundary.WE = BoundaryWE
Boundary.SN = BoundarySN
Boundary.BT = BoundaryBT
Topography=(TopoS=TopoS,
            H=H,
            P1=P1,
            P2=P2,
            P3=P3,
            P4=P4,
            )

# Equation
if Equation == "CompressibleShallow"
  Model.Equation = Models.CompressibleShallow()
elseif Equation == "CompressibleDeep"
  Model.Equation = Models.CompressibleDeep()
end
Grid, Exchange = Grids.InitGridCart(backend,FTB,OrdPoly,nx,ny,Lx,Ly,x0,y0,Boundary,nz,Model,ParallelCom)


#Topography
if TopoS == "AgnesiHill"
  TopoProfile = Examples.AgnesiHill()()
elseif TopoS == "SchaerHill"
  TopoProfile = Examples.SchaerHill()()
else
  TopoProfile = Examples.Flat()()  
end  

Grid.AdaptGrid = Grids.AdaptGrid(FTB,AdaptGridType,H)

@show "InitCart"
(CG, Metric, Global) = DyCore.InitCart(backend,FTB,OrdPoly,OrdPolyZ,H,Topography,Model,
  Phys,TopoProfile,Exchange,Grid,ParallelCom)

# Initial values
#if Problem == "Stratified" || Problem == "HillAgnesiXCart"
#  Profile = Examples.StratifiedExample()(Param,Phys)
#elseif Problem == "WarmBubble2DXCart"
#  Profile = Examples.WarmBubbleCartExample()(Param,Phys)
#elseif Problem == "BryanFritschCart"
#  ProfileBF = Models.TestRes(Phys)
#  Profile = Examples.BryanFritsch(ProfileBF)(Param,Phys)
#end


# Initial values
Examples.InitialProfile!(Model,Problem,Param,Phys)
U = GPU.InitialConditions(backend,FTB,CG,Metric,Phys,Global,Model.InitialProfile,Param)


#Coriolis
if Coriolis
  if Equation == "CompressibleShallow"
    CoriolisFun = GPU.CoriolisShallow()(Phys)
    Model.CoriolisFun = CoriolisFun
  elseif Equation == "CompressibleDeep"
    CoriolisFun = GPU.CoriolisDeep()(Phys)
    Model.CoriolisFun = CoriolisFun
  else
    CoriolisFun = GPU.CoriolisNo()()
    Model.CoriolisFun = CoriolisFun
  end
else
  CoriolisFun = GPU.CoriolisNo()()
  Model.CoriolisFun = CoriolisFun
end

#Buoyancy
if Buoyancy
  if Equation == "CompressibleShallow"
    GravitationFun = GPU.GravitationShallow()(Phys)
    Model.GravitationFun = GravitationFun
  elseif Equation == "CompressibleDeep"
    GravitationFun = GPU.GravitationDeep()(Phys)
    Model.GravitationFun = GravitationFun
  else
    GravitationFun = GPU.GravitationNo()
    Model.GravitationFun = GravitationFun
  end
else
  GravitationFun = GPU.GravitationNo()
  Model.GravitationFun = GravitationFun
end

# Pressure
if State == "Dry"
  Pressure = Models.Dry()(Phys)
  Model.Pressure = Pressure
elseif State == "Moist"
  Pressure = Models.Moist()(Phys,Model.RhoPos,Model.ThPos,
    Model.RhoVPos+NumV,Model.RhoCPos+NumV)
  Model.Pressure = Pressure
end

#Microphysics
if Microphysics
  if TypeMicrophysics == "SimpleMicrophysics"
    MicrophysicsSource  = Models.SimpleMicrophysics()(Phys,Model.RhoPos,Model.ThPos,
      Model.RhoVPos+NumV,Model.RhoCPos+NumV,Model.RelCloud,Model.Rain)
    Model.MicrophysicsSource = MicrophysicsSource
  else
    @show "False Type Microphysics"
  end
end

# Surface flux
Global.SurfaceData = Surfaces.SurfaceData{FTB}(backend,CG.NumG)
Global.LandUseData = Surfaces.LandUseData{FTB}(backend,CG.NumG)
@. Global.LandUseData.z0M = 0.01
@. Global.LandUseData.z0H = 0.01
@. Global.LandUseData.LandClass = 5
if Model.SurfaceFlux || Model.VerticalDiffusion || Model.SurfaceFluxMom || Model.VerticalDiffusionMom
  if SurfaceScheme == ""
    @show "Warning: No surface scheme"  
  elseif SurfaceScheme == "MOST"
    @show "SurfaceScheme MOST"
    SurfaceValues, SurfaceFluxValues = Surfaces.MOSurface()(Surfaces.Businger(),Phys,Model.RhoPos,Model.uPos,
      Model.vPos,Model.wPos,Model.ThPos)
    Model.SurfaceValues = SurfaceValues
    Model.SurfaceFluxValues = SurfaceFluxValues
  end
end
if Model.SurfaceFlux
  if RhoVPos > 0  
    Model.SurfaceFluxRhs = Surfaces.SurfaceFlux(Phys,Param,Model.ThPos,Model.RhoPos,Model.RhoVPos)
  else  
    Model.SurfaceFluxRhs = Surfaces.SurfaceFlux(Phys,Param,Model.ThPos,Model.RhoPos)
  end  
end

#Vertical Diffusion
if Model.VerticalDiffusion || Model.VerticalDiffusionMom
  if Model.Turbulence
    Model.Eddy = Examples.TkeKoefficient()(Param,Phys,TkePos,Model.RhoPos)
  else
    Model.Eddy = Examples.SimpleKoefficient()(Param,Phys)
  end
end

#Turbulence
if Model.Turbulence
  Model.TurbulenceSource = Examples.TKEModel()(Param,Phys,Model.RhoPos,Model.uPos,
    Model.vPos,Model.ThPos,Model.TkePos)
end

# Forcing
Force =  Examples.NoForcing()(Param,Phys)
# Damping
if Damping
  Damp = GPU.DampingW()(FTB(H),FTB(StrideDamp),FTB(Relax),Model.wPos)
  Model.Damp = Damp
end

# Output
Global.Output.vtkFileName=string(Problem*"_")
Global.Output.vtk=0
  Global.Output.Flat=true
  Global.Output.H=H
  if ModelType == "VectorInvariant" || ModelType == "Advection"
    if State == "Dry"
      Global.Output.cNames = [
        "Rho",
        "u",
        "v",
        "wB",
        "Th",
        "Pres",
        ]
      if TkePos > 0
        push!(Global.Output.cNames,"Tke")
      end
      if VerticalDiffusion
        push!(Global.Output.cNames,"DiffKoeff")
      end
    elseif State == "Moist"
      Global.Output.cNames = [
        "Rho",
        "u",
        "v",
        "wB",
        "Th",
        "ThE",
        "Pres",
        "Tr1",
        "Tr2",
        ]
    end    
  elseif ModelType == "Conservative"
    Global.Output.cNames = [
      "Rho",
      "Rhou",
      "Rhov",
      "w",
      "Th",
      "Vort",
      ]
  end

  Global.Output.PrintDays = PrintDays
  Global.Output.PrintHours = PrintHours
  Global.Output.PrintMinutes = PrintMinutes
  Global.Output.PrintSeconds = PrintSeconds
  Global.Output.PrintTime = PrintTime
  Global.Output.PrintStartTime = PrintStartTime
  if OrdPrint == 0
    Global.Output.OrdPrint = CG.OrdPoly
  else  
    Global.Output.OrdPrint = OrdPrint
  end  
  Global.vtkCache = Outputs.vtkStruct{FTB}(backend,Global.Output.OrdPrint,Grids.TransCartX!,CG,Metric,Global)

  # TimeStepper
  time=[0.0]
  Global.TimeStepper.IntMethod = IntMethod
  Global.TimeStepper.Table = Table
  Global.TimeStepper.dtau = dtau
  Global.TimeStepper.SimDays = SimDays
  Global.TimeStepper.SimHours = SimHours
  Global.TimeStepper.SimMinutes = SimMinutes
  Global.TimeStepper.SimSeconds = SimSeconds
  Global.TimeStepper.SimTime = SimTime
  if ModelType == "VectorInvariant" || ModelType == "Advection"
    DiscType = Val(:VectorInvariant)
  elseif ModelType == "Conservative"
    DiscType = Val(:Conservative)
  end

  if Device == "CPU"  || Device == "GPU"
  Global.ParallelCom.NumberThreadGPU = NumberThreadGPU
  nT = max(7 + NumTr, NumV + NumTr)
  Parallels.InitExchangeData3D(backend,FTB,nz,nT,Exchange)
  @show "vor Timestepper"
  Integration.TimeStepper!(U,GPU.FcnGPU!,GPU.FcnPrepareGPU!,DyCore.JacSchurGPU!,
    Grids.TransCartX,CG,Metric,Phys,Exchange,Global,Param,Model.Equation)
else
  nT = max(7 + NumTr, NumV + NumTr)
  Parallels.InitExchangeData3D(backend,FTB,nz,nT,Exchange)
  @show "vor CPU Timestepper"
  Integration.TimeStepper!(U,DyCore.Fcn!,DyCore.FcnPrepare!,DyCore.JacSchurGPU!,
    Grids.TransCartX,CG,Metric,Phys,Exchange,Global,Param,Model.Equation)
end

