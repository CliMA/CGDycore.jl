mutable struct CacheStruct{FT<:AbstractFloat,
                           AT2<:AbstractArray,
                           AT3<:AbstractArray,
                           AT4<:AbstractArray}
PresCG::Array{FT, 2}
AuxG::AT3
KV::AT2
Aux2DG::Array{FT, 3}
Temp::Array{FT, 3}
KE::AT2
cTrS::Array{FT, 3}
FCG::Array{FT, 4}
FCC::Array{FT, 3}
FwCC::Array{FT, 2}
Vn::AT3
RhoCG::Array{FT, 2}
v1CG::Array{FT, 2}
v2CG::Array{FT, 2}
wCG::Array{FT, 2}
Omega::Array{FT, 2}
wCCG::Array{FT, 2}
ThCG::Array{FT, 2}
TrCG::Array{FT, 3}
DivThCG::Array{FT, 2}
DivwCG::Array{FT, 2}
zPG::Array{FT, 2}
pBGrdCG::Array{FT, 2}
RhoBGrdCG::Array{FT, 2}
DivThC::Array{FT, 2}
DivwC::Array{FT, 2}
KVCG::Array{FT, 2}
Temp1::AT3
k::AT4
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
f::AT4
qMin::AT3
qMax::AT3
end

function CacheStruct{FT}(backend) where FT<:AbstractFloat
PresCG=zeros(FT,0,0)
AuxG=KernelAbstractions.zeros(backend,FT,0,0,0)
KV=KernelAbstractions.zeros(backend,FT,0,0)
Aux2DG=zeros(FT,0,0,0)
Temp=zeros(FT,0,0,0)
KE=KernelAbstractions.zeros(backend,FT,0,0)
cTrS=zeros(FT,0,0,0)
FCG=zeros(FT,0,0,0,0)
FCC=zeros(FT,0,0,0)
FwCC=zeros(FT,0,0)
Vn=KernelAbstractions.zeros(backend,FT,0,0,0)
RhoCG=zeros(FT,0,0)
v1CG=zeros(FT,0,0)
v2CG=zeros(FT,0,0)
wCG=zeros(FT,0,0)
Omega=zeros(FT,0,0)
wCCG=zeros(FT,0,0)
ThCG=zeros(FT,0,0)
TrCG=zeros(FT,0,0,0)
DivThCG=zeros(FT,0,0)
DivwCG=zeros(FT,0,0)
zPG=zeros(FT,0,0)
pBGrdCG=zeros(FT,0,0)
RhoBGrdCG=zeros(FT,0,0)
DivThC=zeros(FT,0,0)
DivwC=zeros(FT,0,0)
KVCG=zeros(FT,0,0)
Temp1=KernelAbstractions.zeros(backend,FT,0,0,0)
k=KernelAbstractions.zeros(backend,FT,0,0,0,0)
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
f=KernelAbstractions.zeros(backend,FT,0,0,0,0)
qMin=KernelAbstractions.zeros(backend,FT,0,0,0)
qMax=KernelAbstractions.zeros(backend,FT,0,0,0)
return CacheStruct{FT,
                   typeof(KE),
                   typeof(RhoS),
                   typeof(VS)}(
  PresCG,
  AuxG,
  KV,
  Aux2DG,
  Temp,
  KE,
  cTrS,
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
  DivThCG,
  DivwCG,
  zPG,
  pBGrdCG,
  RhoBGrdCG,
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

function CacheStruct{FT}(backend,DoF,NF,NGF,NumG,nz,NumV,NumTr) where FT<:AbstractFloat
PresCG=zeros(FT,DoF,nz)
AuxG=KernelAbstractions.zeros(backend,FT,nz,NumG,4)
KV=KernelAbstractions.zeros(backend,FT,nz,NumG)
Aux2DG=zeros(FT,1,NumG,NumTr+1)
Temp=zeros(FT,DoF,nz,NF)
KE=KernelAbstractions.zeros(backend,FT,DoF,nz)
cTrS=zeros(FT,DoF,NF,NumTr)
RhoVSurf=KernelAbstractions.zeros(backend,FT,DoF,NF)
CT=KernelAbstractions.zeros(backend,FT,DoF,NF)
CH=KernelAbstractions.zeros(backend,FT,DoF,NF)
FCG=zeros(FT,DoF,nz,NF,NumV+NumTr)
FCC=zeros(FT,DoF,nz,NumV+NumTr)
FwCC=zeros(FT,DoF,nz+1)
Vn=KernelAbstractions.zeros(backend,FT,nz,NumG,NumV+NumTr)
RhoCG=zeros(FT,DoF,nz)
v1CG=zeros(FT,DoF,nz)
v2CG=zeros(FT,DoF,nz)
wCG=zeros(FT,DoF,nz+1)
Omega=zeros(FT,DoF,nz+1)
wCCG=zeros(FT,DoF,nz)
ThCG=zeros(FT,DoF,nz)
TrCG=zeros(FT,DoF,nz,NumTr)
DivThCG=zeros(FT,DoF,nz)
DivwCG=zeros(FT,DoF,nz+1)
zPG=zeros(FT,DoF,nz)
pBGrdCG=zeros(FT,DoF,nz)
RhoBGrdCG=zeros(FT,DoF,nz)
DivThC=zeros(FT,DoF,nz)
DivwC=zeros(FT,DoF,nz+1)
KVCG=zeros(FT,DoF,nz)
Temp1=KernelAbstractions.zeros(backend,FT,nz,NumG,max(NumV+NumTr,7+NumTr))
k=KernelAbstractions.zeros(backend,FT,0,0,0,0)
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
f=KernelAbstractions.zeros(backend,FT,0,0,0,0)
qMin=KernelAbstractions.zeros(backend,FT,nz,NF+NGF,NumTr+1)
qMax=KernelAbstractions.zeros(backend,FT,nz,NF+NGF,NumTr+1)
return CacheStruct{FT,
                   typeof(KE),
                   typeof(RhoS),
                   typeof(VS)}(
  PresCG,
  AuxG,
  KV,
  Aux2DG,
  Temp,
  KE,
  cTrS,
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
  DivThCG,
  DivwCG,
  zPG,
  pBGrdCG,
  RhoBGrdCG,
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
