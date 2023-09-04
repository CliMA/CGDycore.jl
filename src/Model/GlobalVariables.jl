using KernelAbstractions
mutable struct CacheStruct{FT<:Real,
                           AT3<:AbstractArray,
                           AT4<:AbstractArray}
CacheE1::Array{FT, 2}
CacheE2::Array{FT, 2}
CacheE3::Array{FT, 2}
CacheE4::Array{FT, 2}
CacheE5::Array{FT, 2}
CacheF1::Array{FT, 3}
CacheF2::Array{FT, 3}
CacheF3::Array{FT, 3}
CacheF4::Array{FT, 3}
CacheF5::Array{FT, 3}
CacheF6::Array{FT, 3}
CacheC1::SubArray{FT, 3}
CacheC2::SubArray{FT, 3}
CacheC3::SubArray{FT, 3}
CacheC4::SubArray{FT, 3}
CacheC5::SubArray{FT, 3}
CacheC6::SubArray{FT, 3}
Cache1::Array{FT, 2}
Cache2::Array{FT, 2}
Cache3::Array{FT, 2}
Cache4::Array{FT, 2}
PresCG::Array{FT, 3}
AuxG::Array{FT, 3}
Aux2DG::Array{FT, 3}
Temp::Array{FT, 4}
KE::Array{FT, 3}
uStar::Array{FT, 3}
cTrS::Array{FT, 4}
TSurf::Array{FT, 3}
FCG::Array{FT, 5}
FCC::Array{FT, 4}
FwCC::Array{FT, 3}
Vn::Array{FT, 3}
RhoCG::Array{FT, 3}
v1CG::Array{FT, 3}
v2CG::Array{FT, 3}
wCG::Array{FT, 3}
Omega::Array{FT, 3}
wCCG::Array{FT, 3}
ThCG::Array{FT, 3}
TrCG::Array{FT, 4}
Rot1CG::Array{FT, 3}
Rot2CG::Array{FT, 3}
Grad1CG::Array{FT, 3}
Grad2CG::Array{FT, 3}
DivCG::Array{FT, 3}
DivThCG::Array{FT, 3}
DivwCG::Array{FT, 3}
zPG::Array{FT, 3}
pBGrdCG::Array{FT, 3}
RhoBGrdCG::Array{FT, 3}
Rot1C::Array{FT, 3}
Rot2C::Array{FT, 3}
Grad1C::Array{FT, 3}
Grad2C::Array{FT, 3}
DivC::Array{FT, 3}
DivThC::Array{FT, 3}
DivwC::Array{FT, 3}
KVCG::Array{FT, 3}
Temp1::AT3
#Temp1::Array{FT, 3}
k::Array{FT, 4}
Ymyn::Array{FT, 4}
Y::Array{FT, 4}
Z::Array{FT, 4}
fV::AT3
R::Array{FT, 3}
dZ::Array{FT, 3}
fS::AT4
fRhoS::AT3
VS::AT4
RhoS::AT3
f::Array{FT, 4}
qMin::Array{FT, 3}
qMax::Array{FT, 3}
end

function CacheStruct{FT}(backend) where FT<:Real
CacheE1=zeros(FT,0,0);
CacheE2=zeros(FT,0,0);
CacheE3=zeros(FT,0,0);
CacheE4=zeros(FT,0,0);
CacheE5=zeros(FT,0,0);
CacheF1=zeros(FT,0,0,0);
CacheF2=zeros(FT,0,0,0);
CacheF3=zeros(FT,0,0,0);
CacheF4=zeros(FT,0,0,0);
CacheF5=zeros(FT,0,0,0);
CacheF6=zeros(FT,0,0,0);
CacheC1 = view(CacheF1,:,:,:)
CacheC2 = view(CacheF2,:,:,:)
CacheC3 = view(CacheF3,:,:,:)
CacheC4 = view(CacheF4,:,:,:)
CacheC5 = view(CacheF5,:,:,:)
CacheC6 = view(CacheF6,:,:,:)
Cache1=zeros(FT,0,0)
Cache2=zeros(FT,0,0)
Cache3=zeros(FT,0,0)
Cache4=zeros(FT,0,0)
PresCG=zeros(FT,0,0,0)
AuxG=zeros(FT,0,0,0)
Aux2DG=zeros(FT,0,0,0)
Temp=zeros(FT,0,0,0,0)
KE=zeros(FT,0,0,0)
uStar=zeros(FT,0,0,0)
cTrS=zeros(FT,0,0,0,0)
TSurf=zeros(FT,0,0,0)
FCG=zeros(FT,0,0,0,0,0)
FCC=zeros(FT,0,0,0,0)
FwCC=zeros(FT,0,0,0)
Vn=zeros(FT,0,0,0)
RhoCG=zeros(FT,0,0,0)
v1CG=zeros(FT,0,0,0)
v2CG=zeros(FT,0,0,0)
wCG=zeros(FT,0,0,0)
Omega=zeros(FT,0,0,0)
wCCG=zeros(FT,0,0,0)
ThCG=zeros(FT,0,0,0)
TrCG=zeros(FT,0,0,0,0)
Rot1CG=zeros(FT,0,0,0)
Rot2CG=zeros(FT,0,0,0)
Grad1CG=zeros(FT,0,0,0)
Grad2CG=zeros(FT,0,0,0)
DivCG=zeros(FT,0,0,0)
DivThCG=zeros(FT,0,0,0)
DivwCG=zeros(FT,0,0,0)
zPG=zeros(FT,0,0,0)
pBGrdCG=zeros(FT,0,0,0)
RhoBGrdCG=zeros(FT,0,0,0)
Rot1C=zeros(FT,0,0,0)
Rot2C=zeros(FT,0,0,0)
Grad1C=zeros(FT,0,0,0)
Grad2C=zeros(FT,0,0,0)
DivC=zeros(FT,0,0,0)
DivThC=zeros(FT,0,0,0)
DivwC=zeros(FT,0,0,0)
KVCG=zeros(FT,0,0,0)
Temp1=KernelAbstractions.zeros(backend,FT,0,0,0)
#Temp1=zeros(FT,0,0,0)
k=zeros(FT,0,0,0,0)
Ymyn=zeros(FT,0,0,0,0)
Y=zeros(FT,0,0,0,0)
Z=zeros(FT,0,0,0,0)
fV=KernelAbstractions.zeros(backend,FT,0,0,0)
R=zeros(FT,0,0,0)
dZ=zeros(FT,0,0,0)
fS=KernelAbstractions.zeros(backend,FT,0,0,0,0)
fRhoS=KernelAbstractions.zeros(backend,FT,0,0,0)
VS=KernelAbstractions.zeros(backend,FT,0,0,0,0)
RhoS=KernelAbstractions.zeros(backend,FT,0,0,0)
f=zeros(FT,0,0,0,0)
qMin=zeros(FT,0,0,0)
qMax=zeros(FT,0,0,0)
return CacheStruct{FT,
                   typeof(Temp1),
                   typeof(VS)}(
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
  PresCG,
  AuxG,
  Aux2DG,
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
  DivThCG,
  DivwCG,
  zPG,
  pBGrdCG,
  RhoBGrdCG,
  Rot1C,
  Rot2C,
  Grad1C,
  Grad2C,
  DivC,
  DivThC,
  DivwC,
  KVCG,
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

function CacheStruct{FT}(backend,OP,NF,NGF,NumG,nz,NumV,NumTr) where FT<:Real
CacheE1=zeros(FT,OP,OP);
CacheE2=zeros(FT,OP,OP);
CacheE3=zeros(FT,OP,OP);
CacheE4=zeros(FT,OP,OP);
CacheE5=zeros(FT,OP,OP);
CacheF1=zeros(FT,OP,OP,nz+1);
CacheF2=zeros(FT,OP,OP,nz+1);
CacheF3=zeros(FT,OP,OP,nz+1);
CacheF4=zeros(FT,OP,OP,nz+1);
CacheF5=zeros(FT,OP,OP,nz+1);
CacheF6=zeros(FT,OP,OP,nz+1);
CacheC1 = view(CacheF1,:,:,1:nz)
CacheC2 = view(CacheF2,:,:,1:nz)
CacheC3 = view(CacheF3,:,:,1:nz)
CacheC4 = view(CacheF4,:,:,1:nz)
CacheC5 = view(CacheF5,:,:,1:nz)
CacheC6 = view(CacheF6,:,:,1:nz)
Cache1=zeros(FT,nz,NumG)
Cache2=zeros(FT,nz,NumG)
Cache3=zeros(FT,nz,NumG)
Cache4=zeros(FT,nz,NumG)
PresCG=zeros(FT,OP,OP,nz)
AuxG=zeros(FT,nz,NumG,4)
Aux2DG=zeros(FT,1,NumG,NumTr+1)
Temp=zeros(FT,OP,OP,nz,NF)
KE=zeros(FT,OP,OP,nz)
uStar=zeros(FT,OP,OP,NF)
cTrS=zeros(FT,OP,OP,NF,NumTr)
TSurf=zeros(FT,0,0,0)
FCG=zeros(FT,OP,OP,nz,NF,NumV+NumTr)
FCC=zeros(FT,OP,OP,nz,NumV+NumTr)
FwCC=zeros(FT,OP,OP,nz+1)
Vn=zeros(FT,nz,NumG,NumV+NumTr)
RhoCG=zeros(FT,OP,OP,nz)
v1CG=zeros(FT,OP,OP,nz)
v2CG=zeros(FT,OP,OP,nz)
wCG=zeros(FT,OP,OP,nz+1)
Omega=zeros(FT,OP,OP,nz+1)
wCCG=zeros(FT,OP,OP,nz)
ThCG=zeros(FT,OP,OP,nz)
TrCG=zeros(FT,OP,OP,nz,NumTr)
Rot1CG=zeros(FT,OP,OP,nz)
Rot2CG=zeros(FT,OP,OP,nz)
Grad1CG=zeros(FT,OP,OP,nz)
Grad2CG=zeros(FT,OP,OP,nz)
DivCG=zeros(FT,OP,OP,nz)
DivThCG=zeros(FT,OP,OP,nz)
DivwCG=zeros(FT,OP,OP,nz+1)
zPG=zeros(FT,OP,OP,nz)
pBGrdCG=zeros(FT,OP,OP,nz)
RhoBGrdCG=zeros(FT,OP,OP,nz)
Rot1C=zeros(FT,OP,OP,nz)
Rot2C=zeros(FT,OP,OP,nz)
Grad1C=zeros(OP,OP,nz)
Grad2C=zeros(OP,OP,nz)
DivC=zeros(FT,OP,OP,nz)
DivThC=zeros(FT,OP,OP,nz)
DivwC=zeros(FT,OP,OP,nz+1)
KVCG=zeros(FT,OP,OP,nz)
Temp1=KernelAbstractions.zeros(backend,FT,nz,NumG,max(NumV+NumTr,9+NumTr))
#Temp1=zeros(FT,nz,NumG,max(NumV+NumTr,9+NumTr))
k=zeros(FT,0,0,0,0)
Ymyn=zeros(FT,0,0,0,0)
Y=zeros(FT,0,0,0,0)
Z=zeros(FT,0,0,0,0)
fV=KernelAbstractions.zeros(backend,FT,0,0,0)
R=zeros(FT,0,0,0)
dZ=zeros(FT,0,0,0)
fS=KernelAbstractions.zeros(backend,FT,0,0,0,0)
fRhoS=KernelAbstractions.zeros(backend,FT,0,0,0)
VS=KernelAbstractions.zeros(backend,FT,0,0,0,0)
RhoS=KernelAbstractions.zeros(backend,FT,0,0,0)
f=zeros(FT,0,0,0,0)
qMin=zeros(FT,nz,NF+NGF,NumTr+1)
qMax=zeros(FT,nz,NF+NGF,NumTr+1)
return CacheStruct{FT,
                   typeof(Temp1),
                   typeof(VS)}(
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
  PresCG,
  AuxG,
  Aux2DG,
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
  DivThCG,
  DivwCG,
  zPG,
  pBGrdCG,
  RhoBGrdCG,
  Rot1C,
  Rot2C,
  Grad1C,
  Grad2C,
  DivC,
  DivThC,
  DivwC,
  KVCG,
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
function TimeStepperStruct()
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
  Topography::NamedTuple
end
function OutputStruct(Topography::NamedTuple)
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
  Topography,
  )
end  

mutable struct MetricStruct{FT<:Real,
                            AT2<:AbstractArray,
                            AT3<:AbstractArray,
                            AT4<:AbstractArray,
                            AT5<:AbstractArray,
                            AT6<:AbstractArray,
                            AT7<:AbstractArray}
  lat::AT3
  J::AT5
  X::AT6
  dXdxI::AT7
  nS::AT4
  FS::AT3
  dz::AT2
  zP::AT2
end
function MetricStruct{FT}(backend) where FT<:Real
  lat    = KernelAbstractions.zeros(backend,FT,0,0,0)
  J      = KernelAbstractions.zeros(backend,FT,0,0,0,0,0)
  X      = KernelAbstractions.zeros(backend,FT,0,0,0,0,0,0)
  dXdxI  = KernelAbstractions.zeros(backend,FT,0,0,0,0,0,0,0)
  nS = KernelAbstractions.zeros(backend,FT,0,0,0,0)
  FS = KernelAbstractions.zeros(backend,FT,0,0,0)
  dz = KernelAbstractions.zeros(backend,FT,0,0)
  zP = KernelAbstractions.zeros(backend,FT,0,0)
    return MetricStruct{FT,
                        typeof(zP),
                        typeof(FS),
                        typeof(nS),
                        typeof(J),
                        typeof(X),
                        typeof(dXdxI)}(
        lat,
        J,
        X,
        dXdxI,
        nS,
        FS,
        dz,
        zP,
    )
end
function MetricStruct{FT}(backend,OP,OPZ,NF,nz) where FT<:Real
    lat    = KernelAbstractions.zeros(backend,FT,OP,OP,NF)
    J      = KernelAbstractions.zeros(backend,FT,OP,OP,OPZ,nz,NF)
    X      = KernelAbstractions.zeros(backend,FT,OP,OP,OPZ,3,nz,NF)
#   dXdxI  = KernelAbstractions.zeros(backend,FT,OP,OP,OPZ,nz,3,3,NF)
    dXdxI  = KernelAbstractions.zeros(backend,FT,3,3,OPZ,OP,OP,nz,NF)
    nS = KernelAbstractions.zeros(backend,FT,OP,OP,3,NF)
    FS = KernelAbstractions.zeros(backend,FT,OP,OP,NF)
    dz = KernelAbstractions.zeros(backend,FT,0,0)
    zP = KernelAbstractions.zeros(backend,FT,0,0)
    return MetricStruct{FT,
                        typeof(zP),
                        typeof(FS),
                        typeof(nS),
                        typeof(J),
                        typeof(X),
                        typeof(dXdxI)}(
        lat,
        J,
        X,
        dXdxI,
        nS, 
        FS, 
        dz,
        zP,
    )
end

struct PhysParameters{FT<:Real}
  RadEarth::FT 
  Grav::FT 
  Cpd::FT
  Cvd::FT
  Cpv::FT
  Cvv::FT
  Cpl::FT
  Rd::FT
  Rv::FT
  L00::FT
  p0::FT
  Gamma::FT
  kappa::FT
  Omega::FT
  T0::FT
end
function PhysParameters{FT}() where FT<:Real
  RadEarth = 6.37122e+6
  Grav =  9.81
  Cpd = 1004.0
  Cvd = 717.0
  Cpv = 1885.0
  Cvv = 1424.0
  Cpl = 4186.0
  Rd = Cpd - Cvd
  Rv = Cpv - Cvv
# L00 = 2.5000e6 + (Cpl - Cpv) * 273.15
  L00 =  2.5000e6 
  p0 = 1.0e5
  Gamma = Cpd / Cvd
  kappa = Rd / Cpd
  Omega = 2 * pi / 24.0 / 3600.0
  T0 = 273.15
 return PhysParameters{FT}(
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
  T0,
  )
end 

mutable struct ParallelComStruct
  Proc::Int
  ProcNumber::Int
end  
function ParallelComStruct()
  Proc = 1
  ProcNumber = 1
  return ParallelComStruct(
    Proc,
    ProcNumber,
  )
end  

mutable struct ModelStruct
  Problem::String
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
  PertTh::Bool
  RhoVPos::Int
  RhoCPos::Int
  NumV::Int
  NumTr::Int
  Equation::String
  Thermo::String
  ModelType::String
  Source::Bool
  Damping::Bool
  Geos::Bool
  Relax::Float64
  StrideDamp::Float64
  Coriolis::Bool
  CoriolisType::String
  Buoyancy::Bool
  RefProfile::Bool
  HyperVisc::Bool
  HyperDCurl::Float64
  HyperDGrad::Float64
  HyperDRhoDiv::Float64
  HyperDDiv::Float64
  Upwind::Bool
  HorLimit::Bool
  Microphysics::Bool
  RelCloud::Float64
  Rain::Float64
  VerticalDiffusion::Bool
  JacVerticalDiffusion::Bool
  JacVerticalAdvection::Bool
  VerticalDiffusionMom::Bool
  SurfaceFlux::Bool
  SurfaceFluxMom::Bool
  Deep::Bool
  Curl::Bool
  Stretch::Bool
  StretchType::String
end

function Model()
  Problem = ""
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
  PertTh = false
  RhoVPos = 0
  RhoCPos = 0
  NumV = 0
  NumTr = 0
  Equation = "Compressible"
  Thermo = ""
  ModelType = "VectorInvariant"
  Source = false
  Damping = false
  Geos = false
  Relax = 0.0
  StrideDamp = 0.0
  Coriolis = false
  CoriolisType = ""
  Buoyancy = true
  RefProfile = false
  HyperVisc = false
  HyperDCurl = 0.0
  HyperDGrad = 0.0
  HyperDRhoDiv = 0.0
  HyperDDiv = 0.0
  Upwind = false
  HorLimit = false
  Microphysics = false
  RelCloud = 0.0
  Rain = 0.0
  VerticalDiffusion = false
  JacVerticalDiffusion = false
  JacVerticalAdvection = false
  VerticalDiffusionMom = false
  SurfaceFlux = false
  SurfaceFluxMom = false
  Deep = false
  Curl = true
  Stretch = false
  StretchType = ""
  return ModelStruct(
   Problem,
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
   PertTh,
   RhoVPos,
   RhoCPos,
   NumV,
   NumTr,
   Equation,
   Thermo,
   ModelType,
   Source,
   Damping,
   Geos,
   Relax,
   StrideDamp,
   Coriolis,
   CoriolisType,
   Buoyancy,
   RefProfile,
   HyperVisc,
   HyperDCurl,
   HyperDGrad,
   HyperDRhoDiv,
   HyperDDiv,
   Upwind,
   HorLimit,
   Microphysics,
   RelCloud,
   Rain,
   VerticalDiffusion,
   JacVerticalDiffusion,
   JacVerticalAdvection,
   VerticalDiffusionMom,
   SurfaceFlux,
   SurfaceFluxMom,
   Deep,
   Curl,
   Stretch,
   StretchType,
   )
end  

mutable struct GlobalStruct{FT<:Real,
                            TCache}
# Metric::MetricStruct{FT}
  Grid::GridStruct
  Model::ModelStruct
  ParallelCom::ParallelComStruct
  TimeStepper::TimeStepperStruct
  Phys::PhysParameters
  Output::OutputStruct
  Exchange::ExchangeStruct
  vtkCache::vtkStruct{FT}
# Cache::CacheStruct{FT}
  J::JStruct
  latN::Array{Float64, 1}
  ThreadCache::TCache
  ThetaBGrd::Array{Float64, 2}
  TBGrd::Array{Float64, 2}
  pBGrd::Array{Float64, 2}
  RhoBGrd::Array{Float64, 2}
  UGeo::Array{Float64, 2}
  VGeo::Array{Float64, 2}
end
function GlobalStruct{FT}(backend,Grid::GridStruct,
                Model::ModelStruct,
                TimeStepper::TimeStepperStruct,
                ParallelCom::ParallelComStruct,
                Phys::PhysParameters,
                Output::OutputStruct,
                Exchange::ExchangeStruct,
                OP,nz,NumV,NumTr,init_tcache=NamedTuple()) where FT<:Real
# Metric=MetricStruct{FT}(backend)
# Cache=CacheStruct{FT}(backend)
  vtkCache = vtkStruct{FT}(backend)
  J=JStruct()
  latN=zeros(0)
  tcache=(;CreateCache(FT,OP,nz,NumV,NumTr)...,init_tcache)
  ThetaBGrd = zeros(0,0)
  TBGrd = zeros(0,0)
  pBGrd = zeros(0,0)
  RhoBGrd = zeros(0,0)
  UGeo = zeros(0,0)
  VGeo = zeros(0,0)
  return GlobalStruct{FT,typeof(tcache)}(
#   Metric,
    Grid,
    Model,
    ParallelCom,
    TimeStepper,
    Phys,
    Output,
    Exchange,
    vtkCache,
#   Cache,
    J,
    latN,
    tcache,
    ThetaBGrd,
    TBGrd,
    pBGrd,
    RhoBGrd,
    UGeo,
    VGeo,
    )
end  
