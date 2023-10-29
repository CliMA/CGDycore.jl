using KernelAbstractions
mutable struct CacheStruct{FT<:AbstractFloat,
                           AT1<:AbstractArray,
                           AT2<:AbstractArray,
                           AT3<:AbstractArray,
                           AT4<:AbstractArray}
CacheE1::Array{FT, 1}
CacheE2::Array{FT, 1}
CacheE3::Array{FT, 1}
CacheE4::Array{FT, 1}
CacheE5::Array{FT, 1}
CacheF1::Array{FT, 2}
CacheF2::Array{FT, 2}
CacheF3::Array{FT, 2}
CacheF4::Array{FT, 2}
CacheF5::Array{FT, 2}
CacheF6::Array{FT, 2}
CacheC1::SubArray{FT, 2}
CacheC2::SubArray{FT, 2}
CacheC3::SubArray{FT, 2}
CacheC4::SubArray{FT, 2}
CacheC5::SubArray{FT, 2}
CacheC6::SubArray{FT, 2}
Cache1::Array{FT, 2}
Cache2::Array{FT, 2}
Cache3::Array{FT, 2}
Cache4::Array{FT, 2}
PresCG::Array{FT, 2}
AuxG::AT3
Aux2DG::Array{FT, 3}
Temp::Array{FT, 3}
KE::AT2
uStar::AT1
cTrS::Array{FT, 3}
TSurf::Array{FT, 3}
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
Rot1CG::Array{FT, 2}
Rot2CG::Array{FT, 2}
Grad1CG::Array{FT, 2}
Grad2CG::Array{FT, 2}
DivCG::Array{FT, 2}
DivThCG::Array{FT, 2}
DivwCG::Array{FT, 2}
zPG::Array{FT, 2}
pBGrdCG::Array{FT, 2}
RhoBGrdCG::Array{FT, 2}
Rot1C::Array{FT, 2}
Rot2C::Array{FT, 2}
Grad1C::Array{FT, 2}
Grad2C::Array{FT, 2}
DivC::Array{FT, 2}
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
qMin::Array{FT, 3}
qMax::Array{FT, 3}
end

function CacheStruct{FT}(backend) where FT<:AbstractFloat
CacheE1=zeros(FT,0);
CacheE2=zeros(FT,0);
CacheE3=zeros(FT,0);
CacheE4=zeros(FT,0);
CacheE5=zeros(FT,0);
CacheF1=zeros(FT,0,0);
CacheF2=zeros(FT,0,0);
CacheF3=zeros(FT,0,0);
CacheF4=zeros(FT,0,0);
CacheF5=zeros(FT,0,0);
CacheF6=zeros(FT,0,0);
CacheC1 = view(CacheF1,:,:)
CacheC2 = view(CacheF2,:,:)
CacheC3 = view(CacheF3,:,:)
CacheC4 = view(CacheF4,:,:)
CacheC5 = view(CacheF5,:,:)
CacheC6 = view(CacheF6,:,:)
Cache1=zeros(FT,0,0)
Cache2=zeros(FT,0,0)
Cache3=zeros(FT,0,0)
Cache4=zeros(FT,0,0)
PresCG=zeros(FT,0,0)
AuxG=KernelAbstractions.zeros(backend,FT,0,0,0)
Aux2DG=zeros(FT,0,0,0)
Temp=zeros(FT,0,0,0)
KE=KernelAbstractions.zeros(backend,FT,0,0)
uStar=KernelAbstractions.zeros(backend,FT,0,0)
cTrS=zeros(FT,0,0,0)
TSurf=zeros(FT,0,0,0)
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
Rot1CG=zeros(FT,0,0)
Rot2CG=zeros(FT,0,0)
Grad1CG=zeros(FT,0,0)
Grad2CG=zeros(FT,0,0)
DivCG=zeros(FT,0,0)
DivThCG=zeros(FT,0,0)
DivwCG=zeros(FT,0,0)
zPG=zeros(FT,0,0)
pBGrdCG=zeros(FT,0,0)
RhoBGrdCG=zeros(FT,0,0)
Rot1C=zeros(FT,0,0)
Rot2C=zeros(FT,0,0)
Grad1C=zeros(FT,0,0)
Grad2C=zeros(FT,0,0)
DivC=zeros(FT,0,0)
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
qMin=zeros(FT,0,0,0)
qMax=zeros(FT,0,0,0)
return CacheStruct{FT,
                   typeof(uStar),
                   typeof(KE),
                   typeof(RhoS),
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

function CacheStruct{FT}(backend,DoF,NF,NGF,NumG,nz,NumV,NumTr) where FT<:AbstractFloat
CacheE1=zeros(FT,DoF);
CacheE2=zeros(FT,DoF);
CacheE3=zeros(FT,DoF);
CacheE4=zeros(FT,DoF);
CacheE5=zeros(FT,DoF);
CacheF1=zeros(FT,DoF,nz+1);
CacheF2=zeros(FT,DoF,nz+1);
CacheF3=zeros(FT,DoF,nz+1);
CacheF4=zeros(FT,DoF,nz+1);
CacheF5=zeros(FT,DoF,nz+1);
CacheF6=zeros(FT,DoF,nz+1);
CacheC1 = view(CacheF1,:,1:nz)
CacheC2 = view(CacheF2,:,1:nz)
CacheC3 = view(CacheF3,:,1:nz)
CacheC4 = view(CacheF4,:,1:nz)
CacheC5 = view(CacheF5,:,1:nz)
CacheC6 = view(CacheF6,:,1:nz)
Cache1=zeros(FT,nz,NumG)
Cache2=zeros(FT,nz,NumG)
Cache3=zeros(FT,nz,NumG)
Cache4=zeros(FT,nz,NumG)
PresCG=zeros(FT,DoF,nz)
AuxG=KernelAbstractions.zeros(backend,FT,nz,NumG,4)
Aux2DG=zeros(FT,1,NumG,NumTr+1)
Temp=zeros(FT,DoF,nz,NF)
KE=KernelAbstractions.zeros(backend,FT,DoF,nz)
uStar=KernelAbstractions.zeros(backend,FT,DoF,NF)
cTrS=zeros(FT,DoF,NF,NumTr)
TSurf=zeros(FT,0,0,0)
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
Rot1CG=zeros(FT,DoF,nz)
Rot2CG=zeros(FT,DoF,nz)
Grad1CG=zeros(FT,DoF,nz)
Grad2CG=zeros(FT,DoF,nz)
DivCG=zeros(FT,DoF,nz)
DivThCG=zeros(FT,DoF,nz)
DivwCG=zeros(FT,DoF,nz+1)
zPG=zeros(FT,DoF,nz)
pBGrdCG=zeros(FT,DoF,nz)
RhoBGrdCG=zeros(FT,DoF,nz)
Rot1C=zeros(FT,DoF,nz)
Rot2C=zeros(FT,DoF,nz)
Grad1C=zeros(DoF,nz)
Grad2C=zeros(DoF,nz)
DivC=zeros(FT,DoF,nz)
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
qMin=zeros(FT,nz,NF+NGF,NumTr+1)
qMax=zeros(FT,nz,NF+NGF,NumTr+1)
return CacheStruct{FT,
                   typeof(uStar),
                   typeof(KE),
                   typeof(RhoS),
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
