mutable struct CacheStruct
CacheE1::Array{Float64, 2}
CacheE2::Array{Float64, 2}
CacheE3::Array{Float64, 2}
CacheE4::Array{Float64, 2}
CacheE5::Array{Float64, 2}
CacheF1::Array{Float64, 3}
CacheF2::Array{Float64, 3}
CacheF3::Array{Float64, 3}
CacheF4::Array{Float64, 3}
CacheF5::Array{Float64, 3}
CacheF6::Array{Float64, 3}
CacheC1::Array{Float64, 3}
CacheC2::Array{Float64, 3}
CacheC3::Array{Float64, 3}
CacheC4::Array{Float64, 3}
CacheC5::Array{Float64, 3}
CacheC6::Array{Float64, 3}
Cache1::Array{Float64, 2}
Cache2::Array{Float64, 2}
Cache3::Array{Float64, 2}
Cache4::Array{Float64, 2}
Pres::Array{Float64, 4}
PresG::Array{Float64, 2}
Temp::Array{Float64, 4}
KE::Array{Float64, 3}
uStar::Array{Float64, 3}
cTrS::Array{Float64, 4}
TSurf::Array{Float64, 3}
FCG::Array{Float64, 5}
FCC::Array{Float64, 4}
FwCC::Array{Float64, 3}
Vn::Array{Float64, 3}
RhoCG::Array{Float64, 3}
v1CG::Array{Float64, 3}
v2CG::Array{Float64, 3}
wCG::Array{Float64, 3}
Omega::Array{Float64, 3}
wCCG::Array{Float64, 3}
ThCG::Array{Float64, 3}
TrCG::Array{Float64, 4}
Rot1CG::Array{Float64, 3}
Rot2CG::Array{Float64, 3}
Grad1CG::Array{Float64, 3}
Grad2CG::Array{Float64, 3}
DivCG::Array{Float64, 3}
zPG::Array{Float64, 3}
pBGrdCG::Array{Float64, 3}
RhoBGrdCG::Array{Float64, 3}
Rot1C::Array{Float64, 3}
Rot2C::Array{Float64, 3}
Grad1C::Array{Float64, 3}
Grad2C::Array{Float64, 3}
DivC::Array{Float64, 3}
KV::Array{Float64, 3}
Temp1::Array{Float64, 3}
k::Array{Float64, 4}
Ymyn::Array{Float64, 4}
Y::Array{Float64, 4}
Z::Array{Float64, 4}
fV::Array{Float64, 3}
R::Array{Float64, 3}
dZ::Array{Float64, 3}
fS::Array{Float64, 4}
fRhoS::Array{Float64, 3}
VS::Array{Float64, 4}
RhoS::Array{Float64, 3}
f::Array{Float64, 4}
qMin::Array{Float64, 3}
qMax::Array{Float64, 3}
end
function CacheStruct()
CacheE1=zeros(0,0);
CacheE2=zeros(0,0);
CacheE3=zeros(0,0);
CacheE4=zeros(0,0);
CacheE5=zeros(0,0);
CacheF1=zeros(0,0,0);
CacheF2=zeros(0,0,0);
CacheF3=zeros(0,0,0);
CacheF4=zeros(0,0,0);
CacheF5=zeros(0,0,0);
CacheF6=zeros(0,0,0);
CacheC1 = view(CacheF1,:,:,:)
CacheC2 = view(CacheF2,:,:,:)
CacheC3 = view(CacheF3,:,:,:)
CacheC4 = view(CacheF4,:,:,:)
CacheC5 = view(CacheF5,:,:,:)
CacheC6 = view(CacheF6,:,:,:)
Cache1=zeros(0,0)
Cache2=zeros(0,0)
Cache3=zeros(0,0)
Cache4=zeros(0,0)
Pres=zeros(0,0,0,0)
PresG=zeros(0,0)
Temp=zeros(0,0,0,0)
KE=zeros(0,0,0)
uStar=zeros(0,0,0)
cTrS=zeros(0,0,0,0)
TSurf=zeros(0,0,0)
FCG=zeros(0,0,0,0,0)
FCC=zeros(0,0,0,0)
FwCC=zeros(0,0,0)
Vn=zeros(0,0,0)
RhoCG=zeros(0,0,0)
v1CG=zeros(0,0,0)
v2CG=zeros(0,0,0)
wCG=zeros(0,0,0)
Omega=zeros(0,0,0)
wCCG=zeros(0,0,0)
ThCG=zeros(0,0,0)
TrCG=zeros(0,0,0,0)
Rot1CG=zeros(0,0,0)
Rot2CG=zeros(0,0,0)
Grad1CG=zeros(0,0,0)
Grad2CG=zeros(0,0,0)
DivCG=zeros(0,0,0)
zPG=zeros(0,0,0)
pBGrdCG=zeros(0,0,0)
RhoBGrdCG=zeros(0,0,0)
Rot1C=zeros(0,0,0)
Rot2C=zeros(0,0,0)
Grad1C=zeros(0,0,0)
Grad2C=zeros(0,0,0)
DivC=zeros(0,0,0)
KV=zeros(0,0,0)
Temp1=zeros(0,0,0)
k=zeros(0,0,0,0)
Ymyn=zeros(0,0,0,0)
Y=zeros(0,0,0,0)
Z=zeros(0,0,0,0)
fV=zeros(0,0,0)
R=zeros(0,0,0)
dZ=zeros(0,0,0)
fS=zeros(0,0,0,0)
fRhoS=zeros(0,0,0)
VS=zeros(0,0,0,0)
RhoS=zeros(0,0,0)
f=zeros(0,0,0,0)
qMin=zeros(0,0,0)
qMax=zeros(0,0,0)
return CacheStruct(
  CacheE1,
  CacheE2,
  CacheE3,
  CacheE4,
  CacheE5,
  CacheF1,
  CacheF2,
  CacheF3,
  CacheF4,
  CacheF5,
  CacheF6,
  CacheC1,
  CacheC2,
  CacheC3,
  CacheC4,
  CacheC5,
  CacheC6,
  Cache1,
  Cache2,
  Cache3,
  Cache4,
  Pres,
  PresG,
  Temp,
  KE,
  uStar,
  cTrS,
  TSurf,
  FCG,
  FCC,
  FwCC,
  Vn,
  RhoCG,
  v1CG,
  v2CG,
  wCG,
  Omega,
  wCCG,
  ThCG,
  TrCG,
  Rot1CG,
  Rot2CG,
  Grad1CG,
  Grad2CG,
  DivCG,
  zPG,
  pBGrdCG,
  RhoBGrdCG,
  Rot1C,
  Rot2C,
  Grad1C,
  Grad2C,
  DivC,
  KV,
  Temp1,
  k,
  Ymyn,
  Y,
  Z,
  fV,
  R,
  dZ,
  fS,
  fRhoS,
  VS,
  RhoS,
  f,
  qMin,
  qMax
)
end

function CacheCreate(OP,NF,NGF,NumG,nz,NumV,NumTr)
CacheE1=zeros(OP,OP);
CacheE2=zeros(OP,OP);
CacheE3=zeros(OP,OP);
CacheE4=zeros(OP,OP);
CacheE5=zeros(OP,OP);
CacheF1=zeros(OP,OP,nz+1);
CacheF2=zeros(OP,OP,nz+1);
CacheF3=zeros(OP,OP,nz+1);
CacheF4=zeros(OP,OP,nz+1);
CacheF5=zeros(OP,OP,nz+1);
CacheF6=zeros(OP,OP,nz+1);
CacheC1 = view(CacheF1,:,:,1:nz)
CacheC2 = view(CacheF2,:,:,1:nz)
CacheC3 = view(CacheF3,:,:,1:nz)
CacheC4 = view(CacheF4,:,:,1:nz)
CacheC5 = view(CacheF5,:,:,1:nz)
CacheC6 = view(CacheF6,:,:,1:nz)
Cache1=zeros(nz,NumG)
Cache2=zeros(nz,NumG)
Cache3=zeros(nz,NumG)
Cache4=zeros(nz,NumG)
Pres=zeros(OP,OP,nz,NF)
PresG=zeros(nz,NumG)
Temp=zeros(OP,OP,nz,NF)
KE=zeros(OP,OP,nz)
uStar=zeros(OP,OP,NF)
cTrS=zeros(OP,OP,NF,NumTr)
TSurf=zeros(0,0,0)
FCG=zeros(OP,OP,nz,NF,NumV+NumTr)
FCC=zeros(OP,OP,nz,NumV+NumTr)
FwCC=zeros(OP,OP,nz+1)
Vn=zeros(nz,NumG,NumV+NumTr)
RhoCG=zeros(OP,OP,nz)
v1CG=zeros(OP,OP,nz)
v2CG=zeros(OP,OP,nz)
wCG=zeros(OP,OP,nz+1)
Omega=zeros(OP,OP,nz+1)
wCCG=zeros(OP,OP,nz)
ThCG=zeros(OP,OP,nz)
TrCG=zeros(OP,OP,nz,NumTr)
Rot1CG=zeros(OP,OP,nz)
Rot2CG=zeros(OP,OP,nz)
Grad1CG=zeros(OP,OP,nz)
Grad2CG=zeros(OP,OP,nz)
DivCG=zeros(OP,OP,nz)
zPG=zeros(OP,OP,nz)
pBGrdCG=zeros(OP,OP,nz)
RhoBGrdCG=zeros(OP,OP,nz)
Rot1C=zeros(OP,OP,nz)
Rot2C=zeros(OP,OP,nz)
Grad1C=zeros(OP,OP,nz)
Grad2C=zeros(OP,OP,nz)
DivC=zeros(OP,OP,nz)
KV=zeros(OP,OP,nz)
Temp1=zeros(nz,NumG,NumV+NumTr+3)
k=zeros(0,0,0,0)
Ymyn=zeros(0,0,0,0)
Y=zeros(0,0,0,0)
Z=zeros(0,0,0,0)
fV=zeros(0,0,0)
R=zeros(0,0,0)
dZ=zeros(0,0,0)
fS=zeros(0,0,0,0)
fRhoS=zeros(0,0,0)
VS=zeros(0,0,0,0)
RhoS=zeros(0,0,0)
f=zeros(0,0,0,0)
qMin=zeros(nz,NF+NGF,NumTr)
qMax=zeros(nz,NF+NGF,NumTr)
return CacheStruct(
  CacheE1,
  CacheE2,
  CacheE3,
  CacheE4,
  CacheE5,
  CacheF1,
  CacheF2,
  CacheF3,
  CacheF4,
  CacheF5,
  CacheF6,
  CacheC1,
  CacheC2,
  CacheC3,
  CacheC4,
  CacheC5,
  CacheC6,
  Cache1,
  Cache2,
  Cache3,
  Cache4,
  Pres,
  PresG,
  Temp,
  KE,
  uStar,
  cTrS,
  TSurf,
  FCG,
  FCC,
  FwCC,
  Vn,
  RhoCG,
  v1CG,
  v2CG,
  wCG,
  Omega,
  wCCG,
  ThCG,
  TrCG,
  Rot1CG,
  Rot2CG,
  Grad1CG,
  Grad2CG,
  DivCG,
  zPG,
  pBGrdCG,
  RhoBGrdCG,
  Rot1C,
  Rot2C,
  Grad1C,
  Grad2C,
  DivC,
  KV,
  Temp1,
  k,
  Ymyn,
  Y,
  Z,
  fV,
  R,
  dZ,
  fS,
  fRhoS,
  VS,
  RhoS,
  f,
  qMin,
  qMax,
)
end
mutable struct TimeStepperStruct
  IntMethod::String
  Table::String
  dtau::Float64
  dtauStage::Float64
  SimDays::Int
  SimHours::Int
  SimMinutes::Int
  SimSeconds::Int
  SimTime::Float64
  ROS::RosenbrockStruct
  LinIMEX::LinIMEXStruct
  IMEX::IMEXStruct
  MIS::MISStruct
  RK::RungeKuttaStruct
  SSP::SSPRungeKuttaStruct
end
function TimeStepper()
  IntMethod = ""
  Table = ""
  dtau  = 0.0
  dtauStage  = 0.0
  SimDays = 0
  SimHours = 0
  SimMinutes = 0
  SimSeconds = 0
  SimTime = 0.0
  ROS=RosenbrockMethod()
  LinIMEX=LinIMEXMethod()
  IMEX=IMEXMethod()
  MIS=MISMethod()
  RK=RungeKuttaMethod()
  SSP=SSPRungeKuttaMethod()
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
  PrintSeconds::Int
  PrintTime::Float64
  PrintStartDays::Int
  PrintInt::Int
  PrintStartInt::Int
  RadPrint::Float64
  H::Float64
  OrdPrint::Int
  Topography::NamedTuple
end
function Output(Topography::NamedTuple)
  vtk=0
  vtkFileName=""
  Flat=false
  cNames=[]
  nPanel=1
  nIter = 0
  PrintDays = 0
  PrintHours = 0
  PrintSeconds = 0
  PrintTime = 0
  PrintStartDays = 0
  PrintInt = 0
  PrintStartInt = 0
  RadPrint=1000.0
  H=1000.0
  OrdPrint=1
  return OutputStruct(
  vtk,
  vtkFileName,
  Flat,
  cNames,
  nPanel,
  nIter,
  PrintDays,
  PrintHours,
  PrintSeconds,
  PrintTime,
  PrintStartDays,
  PrintInt,
  PrintStartInt,
  RadPrint,
  H,
  OrdPrint,
  Topography,
  )
end  

mutable struct MetricStruct
  lat::Array{Float64, 3}
  JC::Array{Float64, 4}
  JF::Array{Float64, 4}
  J::Array{Float64, 5}
  X::Array{Float64, 6}
  dXdxIF::Array{Float64, 6}
  dXdxI::Array{Float64, 7}
  dXdxIC::Array{Float64, 6}
  nS::Array{Float64, 4}
  dz::Array{Float64, 2}
  zP::Array{Float64, 2}
end
function MetricStruct()
    lat    = zeros(0,0,0)
    JC     = zeros(0,0,0,0)
    JF     = zeros(0,0,0,0)
    J      = zeros(0,0,0,0,0)
    X      = zeros(0,0,0,0,0,0)
    dXdxIF = zeros(0,0,0,0,0,0)
    dXdxI  = zeros(0,0,0,0,0,0,0)
    dXdxIC = zeros(0,0,0,0,0,0)
    nS = zeros(0,0,0,0)
    dz = zeros(0,0)
    zP = zeros(0,0)
    return MetricStruct(
        lat,
        JC,
        JF,
        J,
        X,
        dXdxIF,
        dXdxI,
        dXdxIC,
        nS,
        dz,
        zP,
    )
end
function Metric(OP,OPZ,NF,nz)
    lat    = zeros(OP,OP,NF)
    JC     = zeros(OP,OP,nz,NF)
    JF     = zeros(OP,OP,nz+1,NF)
    J      = zeros(OP,OP,OPZ,nz,NF)
    X      = zeros(OP,OP,OPZ,3,nz,NF)
    dXdxIF = zeros(OP,OP,nz+1,3,3,NF)
    dXdxI  = zeros(OP,OP,OPZ,nz,3,3,NF)
    dXdxIC = zeros(OP,OP,nz,3,3,NF)
    nS = zeros(OP,OP,3,NF)
    dz = zeros(0,0)
    zP = zeros(0,0)
    return MetricStruct(
        lat,
        JC,
        JF,
        J,
        X,
        dXdxIF,
        dXdxI,
        dXdxIC,
        nS, 
        dz,
        zP,
    )
end

mutable struct PhysParameters
  RadEarth::Float64 
  Grav::Float64 
  Cpd::Float64
  Cvd::Float64
  Cpv::Float64
  Cvv::Float64
  Cpl::Float64
  Rd::Float64
  Rv::Float64
  L00::Float64
  p0::Float64
  Gamma::Float64
  kappa::Float64
  Omega::Float64
end
function PhysParameters()
  RadEarth = 6.37122e+6
  Grav = 9.81e0
  Cpd=1004.0e0
  Cvd=717.0e0
  Cpv=1885.0e0
  Cvv=1424.0e0
  Cpl=4186.0e0
  Rd=Cpd-Cvd
  Rv=Cpv-Cvv
  L00 = 2.5000e6 + (Cpl - Cpv) * 273.15
  p0=1.0e5
  Gamma=Cpd/Cvd
  kappa=Rd/Cpd
  Omega=2*pi/24.0/3600.0
 return PhysParameters(
  RadEarth,
  Grav,
  Cpd,
  Cvd,
  Cpv,
  Cvv,
  Cpl,
  Rd,
  Rv,
  L00,
  p0,
  Gamma,
  kappa,
  Omega,
  )
end 
mutable struct ParallelComStruct
  Proc::Int
  ProcNumber::Int
end  

mutable struct ModelStruct
  Problem::String
  ProfRho::String
  ProfTheta::String
  ProfTr::String
  ProfVel::String
  ProfVelW::String
  ProfpBGrd::String
  ProfRhoBGrd::String
  RhoPos::Int
  uPos::Int
  vPos::Int
  wPos::Int
  ThPos::Int
  RhoVPos::Int
  RhoCPos::Int
  NumV::Int
  NumTr::Int
  Equation::String
  Thermo::String
  ModelType::String
  Source::Bool
  Damping::Bool
  Relax::Float64
  StrideDamp::Float64
  Coriolis::Bool
  CoriolisType::String
  Buoyancy::Bool
  RefProfile::Bool
  HyperVisc::Bool
  HyperDCurl::Float64
  HyperDGrad::Float64
  HyperDDiv::Float64
  Upwind::Bool
  HorLimit::Bool
  Microphysics::Bool
  RelCloud::Float64
  Rain::Float64
  VerticalDiffusion::Bool
  SurfaceFlux::Bool
  Deep::Bool
end
function ParallelCom()
  Proc = 1
  ProcNumber = 1
  return ParallelComStruct(
    Proc,
    ProcNumber,
  )
end  

function Model()
  Problem=""
  ProfRho=""
  ProfTheta=""
  ProfTr=""
  ProfVel=""
  ProfVelW=""
  ProfpBGrd=""
  ProfRhoBGrd=""
  RhoPos = 0
  uPos = 0
  vPos = 0
  wPos = 0
  ThPos = 0
  RhoVPos = 0
  RhoCPos = 0
  NumV = 0
  NumTr = 0
  Equation="Compressible"
  Thermo=""
  ModelType="Curl"
  Source=false
  Damping=false
  Relax=0.0
  StrideDamp=0.0
  Coriolis=false
  CoriolisType=""
  Buoyancy=true
  RefProfile=false
  HyperVisc=false
  HyperDCurl=0.0
  HyperDGrad=0.0
  HyperDDiv=0.0
  Upwind=false
  HorLimit=false
  Microphysics=false
  RelCloud=0.0
  Rain=1.0
  VerticalDiffusion=false
  SurfaceFlux=false
  Deep=false
  return ModelStruct(
   Problem,
   ProfRho,
   ProfTheta,
   ProfTr,
   ProfVel,
   ProfVelW,
   ProfpBGrd,
   ProfRhoBGrd,
   RhoPos,
   uPos,
   vPos,
   wPos,
   ThPos,
   RhoVPos,
   RhoCPos,
   NumV,
   NumTr,
   Equation,
   Thermo,
   ModelType,
   Source,
   Damping,
   Relax,
   StrideDamp,
   Coriolis,
   CoriolisType,
   Buoyancy,
   RefProfile,
   HyperVisc,
   HyperDCurl,
   HyperDGrad,
   HyperDDiv,
   Upwind,
   HorLimit,
   Microphysics,
   RelCloud,
   Rain,
   VerticalDiffusion,
   SurfaceFlux,
   Deep,

   )
end  

mutable struct GlobalStruct{TCache}
  Metric::MetricStruct
  Grid::GridStruct
  Model::ModelStruct
  ParallelCom::ParallelComStruct
  TimeStepper::TimeStepperStruct
  Phys::PhysParameters
  Output::OutputStruct
  Exchange::ExchangeStruct
  vtkCache::vtkStruct
  Cache::CacheStruct
  J::JStruct
  latN::Array{Float64, 1}
  ThreadCache::TCache
  ThetaBGrd::Array{Float64, 2}
  TBGrd::Array{Float64, 2}
  pBGrd::Array{Float64, 2}
  RhoBGrd::Array{Float64, 2}
end
function Global(Grid::GridStruct,
                Model::ModelStruct,
                TimeStepper::TimeStepperStruct,
                ParallelCom::ParallelComStruct,
                Phys::PhysParameters,
                Output::OutputStruct,
                Exchange::ExchangeStruct,
                OP,nz,NumV,NumTr,init_tcache=NamedTuple())
  Metric=MetricStruct()
  Cache=CacheStruct()
  vtkCache = vtkStruct()
  J=JStruct()
  latN=zeros(0)
  tcache=(;CreateCache(OP,nz,NumV,NumTr)...,init_tcache)
  ThetaBGrd = zeros(0,0)
  TBGrd = zeros(0,0)
  pBGrd = zeros(0,0)
  RhoBGrd = zeros(0,0)
  return GlobalStruct{typeof(tcache)}(
    Metric,
    Grid,
    Model,
    ParallelCom,
    TimeStepper,
    Phys,
    Output,
    Exchange,
    vtkCache,
    Cache,
    J,
    latN,
    tcache,
    ThetaBGrd,
    TBGrd,
    pBGrd,
    RhoBGrd,
    )
end  
