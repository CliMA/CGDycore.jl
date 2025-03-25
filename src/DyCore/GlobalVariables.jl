mutable struct TimeStepperStruct{FT<:AbstractFloat}
  IntMethod::String
  Table::String
  dtau::Float64
  dtauStage::Float64
  SimDays::Int
  SimHours::Int
  SimMinutes::Int
  SimSeconds::Int
  SimTime::Float64
  ROS::Integration.RosenbrockStruct{FT}
  LinIMEX::Integration.LinIMEXStruct
  IMEX::Integration.IMEXStruct
  MIS::Integration.MISStruct
  RK::Integration.RungeKuttaStruct
  SSP::Integration.SSPRungeKuttaStruct
end
function TimeStepperStruct{FT}(backend) where FT<:AbstractFloat
  IntMethod = ""
  Table = ""
  dtau  = 0.0
  dtauStage  = 0.0
  SimDays = 0
  SimHours = 0
  SimMinutes = 0
  SimSeconds = 0
  SimTime = 0.0
  ROS=Integration.RosenbrockStruct{FT}()
  LinIMEX=Integration.LinIMEXMethod()
  IMEX=Integration.IMEXMethod()
  MIS=Integration.MISMethod()
  RK=Integration.RungeKuttaMethod()
  SSP=Integration.SSPRungeKuttaMethod()
  return TimeStepperStruct(
    IntMethod,
    Table,
    dtau,
    dtauStage,
    SimDays,
    SimHours,
    SimMinutes,
    SimSeconds,
    SimTime,
    ROS,
    LinIMEX,
    IMEX,
    MIS,
    RK,
    SSP,
  )  
end

mutable struct OutputStruct
  vtk::Int
  vtkFileName::String
  Flat::Bool
  cNames::Array{String, 1}
  nPanel::Int
  nIter::Int
  PrintDays::Int
  PrintHours::Int
  PrintMinutes::Int
  PrintSeconds::Int
  PrintTime::Float64
  PrintStartTime::Float64
  StartAverageDays::Int
  PrintInt::Int
  PrintStartInt::Int
  RadPrint::Float64
  H::Float64
  OrdPrint::Int
  OrdPrintZ::Int
  dTol::Float64
end
function OutputStruct()
  vtk=0
  vtkFileName=""
  Flat=false
  cNames=[]
  nPanel=1
  nIter = 0
  PrintDays = 0
  PrintHours = 0
  PrintMinutes = 0
  PrintSeconds = 0
  PrintTime = 0
  PrintStartTime = 0
  StartAverageDays = -1
  PrintInt = 0
  PrintStartInt = 0
  RadPrint=1000.0
  H=1000.0
  OrdPrint=1
  OrdPrintZ=1
  dTol=1.e-8
  return OutputStruct(
  vtk,
  vtkFileName,
  Flat,
  cNames,
  nPanel,
  nIter,
  PrintDays,
  PrintHours,
  PrintMinutes,
  PrintSeconds,
  PrintTime,
  PrintStartTime,
  StartAverageDays,
  PrintInt,
  PrintStartInt,
  RadPrint,
  H,
  OrdPrint,
  OrdPrintZ,
  dTol,
  )
end  

mutable struct MetricStruct{FT<:AbstractFloat,
                            AT2<:AbstractArray,
                            AT3<:AbstractArray,
                            AT4<:AbstractArray,
                            AT5<:AbstractArray,
                            AT6<:AbstractArray}
  J::AT4
  X::AT5
  dXdxI::AT6
  dXdx::AT6
  Rotate::AT6
  nSS::AT2
  nS::AT3
  FS::AT2
  dz::AT2
  zP::AT2
  JC::AT3
  JCW::AT3
  xS::AT2
  VolSurfH::AT4
  NH::AT5
  T1H::AT5
  T2H::AT5
  VolSurfV::AT3
  NV::AT4
  T1V::AT4
  T2V::AT4
end
function MetricStruct{FT}(backend,nQuad,OPZ,NF,nz,NumG) where FT<:AbstractFloat
    J      = KernelAbstractions.zeros(backend,FT,nQuad,OPZ,nz,NF)
    X      = KernelAbstractions.zeros(backend,FT,nQuad,OPZ,3,nz,NF)
    dXdxI  = KernelAbstractions.zeros(backend,FT,3,3,OPZ,nQuad,nz,NF)
    dXdx   = KernelAbstractions.zeros(backend,FT,3,3,OPZ,nQuad,nz,NF)
    Rotate  = KernelAbstractions.zeros(backend,FT,3,3,OPZ,nQuad,nz,NF)
    nSS  = KernelAbstractions.zeros(backend,FT,3,NumG)
    nS = KernelAbstractions.zeros(backend,FT,nQuad,3,NF)
    FS = KernelAbstractions.zeros(backend,FT,nQuad,NF)
    dz = KernelAbstractions.zeros(backend,FT,0,0)
    zP = KernelAbstractions.zeros(backend,FT,0,0)
    JC     = KernelAbstractions.zeros(backend,FT,0,0,0)
    JCW    = KernelAbstractions.zeros(backend,FT,0,0,0)
    xS    = KernelAbstractions.zeros(backend,FT,2,NumG)
    VolSurfH = KernelAbstractions.zeros(backend,FT,0,0,0,0)
    NH = KernelAbstractions.zeros(backend,FT,0,0,0,0,0)
    T1H = KernelAbstractions.zeros(backend,FT,0,0,0,0,0)
    T2H = KernelAbstractions.zeros(backend,FT,0,0,0,0,0)
    VolSurfV = KernelAbstractions.zeros(backend,FT,0,0,0)
    NV = KernelAbstractions.zeros(backend,FT,0,0,0,0)
    T1V = KernelAbstractions.zeros(backend,FT,0,0,0,0)
    T2V = KernelAbstractions.zeros(backend,FT,0,0,0,0)
    return MetricStruct{FT,
                        typeof(zP),
                        typeof(nS),
                        typeof(J),
                        typeof(X),
                        typeof(dXdxI)}(
        J,
        X,
        dXdxI,
        dXdx,
        Rotate,
        nSS,
        nS, 
        FS, 
        dz,
        zP,
        JC,
        JCW,
        xS,
        VolSurfH,
        NH,
        T1H,
        T2H,
        VolSurfV,
        NV,
        T1V,
        T2V,
    )
end

struct PhysParameters{FT<:AbstractFloat}
  RadEarth::FT 
  Grav::FT 
  Cpd::FT
  Cvd::FT
  Cpv::FT
  Cvv::FT
  Cpl::FT
  Cpi::FT
  Rd::FT
  Rv::FT
  L0V::FT
  L0S::FT
  L0F::FT
  p0::FT
  Rho0::FT
  Gamma::FT
  kappa::FT
  Omega::FT
  T0::FT
  T00::FT
  Cd::FT # Dissipation Koefficient
  PrTke::FT # Prandtl number for Tke
  Karm::FT
  cS::FT
end
function PhysParameters{FT}() where FT<:AbstractFloat
  RadEarth::FT = 6.37122e+6
  Grav::FT =  9.80616
  Cpd::FT = 1004.0
  Cvd::FT = 717.0
  Cpv::FT = 1885.0
  Cvv::FT = 1424.0
  Cpl::FT = 4186.0
  Cpi::FT = 2110.0
  Rd::FT = Cpd - Cvd
  Rv::FT = Cpv - Cvv
# L00 = 2.5000e6 + (Cpl - Cpv) * 273.15
  L0V::FT =  2.5000e6 # 2500800 
  L0S::FT =  2.834e6
  L0F::FT =  L0S - L0V
  p0::FT = 1.0e5
  Rho0::FT = 1.41e0
  Gamma::FT = Cpd / Cvd
  kappa::FT = Rd / Cpd
  Omega::FT = 2 * pi / 24.0 / 3600.0
  T0::FT = 273.15
  T00::FT = 273.15 -35.0
  Cd::FT = 0.125
  PrTke::FT = 1.0
  Karm::FT = 0.4
  cS::FT = 360
 return PhysParameters{FT}(
  RadEarth,
  Grav,
  Cpd,
  Cvd,
  Cpv,
  Cvv,
  Cpl,
  Cpi,
  Rd,
  Rv,
  L0V,
  L0S,
  L0F,
  p0,
  Rho0,
  Gamma,
  kappa,
  Omega,
  T0,
  T00,
  Cd,
  PrTke,
  Karm,
  cS,
  )
end 

mutable struct ParallelComStruct
  Proc::Int
  ProcNumber::Int
  NumberThreadGPU::Int
  NumberThreadTriGPU::Int
end  
function ParallelComStruct()
  Proc = 1
  ProcNumber = 1
  NumberThreadGPU = 256
  NumberThreadTriGPU = 32
  return ParallelComStruct(
    Proc,
    ProcNumber,
    NumberThreadGPU,
    NumberThreadTriGPU,
  )
end  

Base.@kwdef mutable struct ModelStruct{FT}
  Problem::String
  Discretization::String
  Profile::Bool
  ProfRho::String
  ProfTheta::String
  ProfTr::String
  ProfVel::String
  ProfVelGeo::String
  ProfVelW::String
  ProfpBGrd::String
  ProfRhoBGrd::String
  ProfTest::String
  RhoPos::Int
  uPos::Int
  vPos::Int
  wPos::Int
  ThPos::Int
  RhoTPos::Int
  RhoVPos::Int
  RhoCPos::Int
  RhoIPos::Int
  RhoRPos::Int
  TkePos::Int
  TrPos::Int
  EDMFPos::Int
  RhoEDMFPos::Int
  wEDMFPos::Int
  ThEDMFPos::Int
  NumV::Int
  NumTr::Int
  NumThermo::Int
  Equation::Models.EquationType
  Thermo::String
  ModelType::String
  Source::Bool
  Forcing::Bool
  Damping::Bool
  Geos::Bool
  Relax::FT
  StrideDamp::FT
  Coriolis::Bool
  CoriolisType::String
  Buoyancy::Bool
  Turbulence::Bool
  EDMF::Bool
  NDEDMF::Int
  RefProfile::Bool
  HyperVisc::Bool
  HyperDCurl::FT
  HyperDGrad::FT
  HyperDRhoDiv::FT
  HyperDDiv::FT
  HyperDDivW::FT
  Upwind::Bool
  HorLimit::Bool
  Microphysics::Bool
  Sedimentation::Bool
  TypeMicrophysics::String
  RelCloud::FT
  Rain::FT
  JacVerticalDiffusion::Bool
  JacVerticalAdvection::Bool
  VerticalDiffusion::Bool
  VerticalDiffusionMom::Bool
  SurfaceFlux::Bool
  SurfaceFluxMom::Bool
  Deep::Bool
  Curl::Bool
  Stretch::Bool
  StretchType::String
  State::String
  InitialProfile::Any
  Eddy::Any
  Force::Any
  Damp::Any
  Pressure::Any
  RiemannSolver::Any
  FluxAverage::Any
  dPresdRhoTh::Any
  dPresdRho::Any
  CoriolisFun::Any
  GravitationFun::Any
  MicrophysicsSource::Any
  SedimentationSource::Any
  TurbulenceSource::Any
  SurfaceFluxRhs::Any
  SurfaceValues::Any
  SurfaceFluxValues::Any
end

function ModelStruct{FT}() where FT <:AbstractFloat
  Problem = ""
  Discretization = "CG"
  Profile = false
  ProfRho = ""
  ProfTheta = ""
  ProfTr = ""
  ProfVel = ""
  ProfVelGeo = ""
  ProfVelW = ""
  ProfpBGrd = ""
  ProfRhoBGrd = ""
  ProfTest = ""
  RhoPos = 0
  uPos = 0
  vPos = 0
  wPos = 0
  ThPos = 0
  RhoTPos = 0
  RhoVPos = 0
  RhoCPos = 0
  RhoIPos = 0
  RhoRPos = 0
  TkePos = 0
  TrPos = 0
  EDMFPos = 0
  RhoEDMFPos = 0
  wEDMFPos = 0
  ThEDMFPos = 0
  NumV = 0
  NumTr = 0
  NumThermo = 0
  Equation = Models.CompressibleShallow()
  Thermo = ""
  ModelType = "VectorInvariant"
  Source = false
  Forcing = false
  Damping = false
  Geos = false
  Relax = 0.0
  StrideDamp = 0.0
  Coriolis = false
  CoriolisType = ""
  Buoyancy = true
  Turbulence = false
  EDMF = false
  NDEDMF = 0
  RefProfile = false
  HyperVisc = false
  HyperDCurl = 0.0
  HyperDGrad = 0.0
  HyperDRhoDiv = 0.0
  HyperDDiv = 0.0
  HyperDDivW = 0.0
  Upwind = false
  HorLimit = false
  Microphysics = false
  Sedimentation = false
  TypeMicrophysics = ""
  RelCloud = 0.0
  Rain = 0.0
  JacVerticalDiffusion = false
  JacVerticalAdvection = false
  VerticalDiffusion = false
  VerticalDiffusionMom = false
  SurfaceFlux = false
  SurfaceFluxMom = false
  Deep = false
  Curl = true
  Stretch = false
  StretchType = ""
  State = ""
  InitialProfile = ""
  Eddy = ""
  Force = ""
  Damp = ""
  Pressure = ""
  RiemannSolver = ""
  FluxAverage = ""
  dPresdRhoTh = ""
  dPresdRho = ""
  CoriolisFun = ""
  GravitationFun = ""
  MicrophysicsSource = ""
  SedimentationSource = ""
  TurbulenceSource = ""
  SurfaceFluxRhs = ""
  SurfaceValues = ""
  SurfaceFluxValues = ""
  return ModelStruct{FT}(
   Problem,
   Discretization,
   Profile,
   ProfRho,
   ProfTheta,
   ProfTr,
   ProfVel,
   ProfVelGeo,
   ProfVelW,
   ProfpBGrd,
   ProfRhoBGrd,
   ProfTest,
   RhoPos,
   uPos,
   vPos,
   wPos,
   ThPos,
   RhoTPos,
   RhoVPos,
   RhoCPos,
   RhoIPos,
   RhoRPos,
   TkePos,
   TrPos,
   EDMFPos,
   RhoEDMFPos,
   wEDMFPos,
   ThEDMFPos,
   NumV,
   NumTr,
   NumThermo,
   Equation,
   Thermo,
   ModelType,
   Source,
   Forcing,
   Damping,
   Geos,
   Relax,
   StrideDamp,
   Coriolis,
   CoriolisType,
   Buoyancy,
   Turbulence,
   EDMF,
   NDEDMF,
   RefProfile,
   HyperVisc,
   HyperDCurl,
   HyperDGrad,
   HyperDRhoDiv,
   HyperDDiv,
   HyperDDivW,
   Upwind,
   HorLimit,
   Microphysics,
   Sedimentation,
   TypeMicrophysics,
   RelCloud,
   Rain,
   JacVerticalDiffusion,
   JacVerticalAdvection,
   VerticalDiffusion,
   VerticalDiffusionMom,
   SurfaceFlux,
   SurfaceFluxMom,
   Deep,
   Curl,
   Stretch,
   StretchType,
   State,
   InitialProfile,
   Eddy,
   Force,
   Damp,
   Pressure,
   RiemannSolver,
   FluxAverage,
   dPresdRhoTh,
   dPresdRho,
   CoriolisFun,
   GravitationFun,
   MicrophysicsSource,
   SedimentationSource,
   TurbulenceSource,
   SurfaceFluxRhs,
   SurfaceValues,
   SurfaceFluxValues,
   )
end  

mutable struct GlobalStruct{FT<:AbstractFloat,
                            AT2<:AbstractArray}
  Grid::Grids.GridStruct
  Model::ModelStruct{FT}
  ParallelCom::ParallelComStruct
  TimeStepper::TimeStepperStruct
  Output::OutputStruct
  vtkCache::Outputs.vtkStruct{FT}
  SurfaceData::Surfaces.SurfaceData{FT}
  LandUseData::Surfaces.LandUseData{FT}
  latN::Array{Float64, 1}
  ThetaBGrd::Array{Float64, 2}
  TBGrd::Array{Float64, 2}
  pBGrd::Array{Float64, 2}
  RhoBGrd::Array{Float64, 2}
  UGeo::AT2
  VGeo::AT2
end
function GlobalStruct{FT}(backend,Grid::Grids.GridStruct,
                Model::ModelStruct,
                TimeStepper::TimeStepperStruct,
                ParallelCom::ParallelComStruct,
                Output::OutputStruct,
                DoF,nz,NumV,NumTr) where FT<:AbstractFloat
  vtkCache = Outputs.vtkStruct{FT}(backend)
  latN=zeros(0)
  ThetaBGrd = zeros(0,0)
  TBGrd = zeros(0,0)
  pBGrd = zeros(0,0)
  RhoBGrd = zeros(0,0)
  UGeo = KernelAbstractions.zeros(backend,FT,0,0)
  VGeo = KernelAbstractions.zeros(backend,FT,0,0)
  SurfaceData = Surfaces.SurfaceData{FT}(backend,0,0)
  LandUseData = Surfaces.LandUseData{FT}(backend,0)
  return GlobalStruct{FT,
                      typeof(UGeo)}(
    Grid,
    Model,
    ParallelCom,
    TimeStepper,
    Output,
    vtkCache,
    SurfaceData,
    LandUseData,
    latN,
    ThetaBGrd,
    TBGrd,
    pBGrd,
    RhoBGrd,
    UGeo,
    VGeo,
    )
end  
