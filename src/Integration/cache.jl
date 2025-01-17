mutable struct CacheStruct{FT<:AbstractFloat,
                           AT2<:AbstractArray,
                           AT3<:AbstractArray,
                           AT4<:AbstractArray}
AuxG::AT3
KV::AT2
Aux2DG::Array{FT, 3}
Temp::Array{FT, 3}
KE::AT2
cTrS::Array{FT, 3}
Vn::AT3
zPG::Array{FT, 2}
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
q::AT3
#EDMF variables
RhoEDMF::AT3
end

function CacheStruct{FT}(backend) where FT<:AbstractFloat
AuxG=KernelAbstractions.zeros(backend,FT,0,0,0)
KV=KernelAbstractions.zeros(backend,FT,0,0)
Aux2DG=zeros(FT,0,0,0)
Temp=zeros(FT,0,0,0)
KE=KernelAbstractions.zeros(backend,FT,0,0)
cTrS=zeros(FT,0,0,0)
Vn=KernelAbstractions.zeros(backend,FT,0,0,0)
zPG=zeros(FT,0,0)
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
q=KernelAbstractions.zeros(backend,FT,0,0,0)
RhoEDMF=KernelAbstractions.zeros(backend,FT,0,0,0)
return CacheStruct{FT,
                   typeof(KE),
                   typeof(RhoS),
                   typeof(VS)}(
  AuxG,
  KV,
  Aux2DG,
  Temp,
  KE,
  cTrS,
  Vn,
  zPG,
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
  q,
  RhoEDMF
)
end

function CacheStruct{FT}(backend,DoF,NF,NGF,NumG,nz,NumV,NumTr,ND) where FT<:AbstractFloat
AuxG=KernelAbstractions.zeros(backend,FT,nz,NumG,4)
@. AuxG[:,:,2] = FT(250.0)
KV=KernelAbstractions.zeros(backend,FT,nz,NumG)
Aux2DG=zeros(FT,1,NumG,NumTr+1)
Temp=zeros(FT,DoF,nz,NF)
KE=KernelAbstractions.zeros(backend,FT,DoF,nz)
cTrS=zeros(FT,DoF,NF,NumTr)
RhoVSurf=KernelAbstractions.zeros(backend,FT,DoF,NF)
CT=KernelAbstractions.zeros(backend,FT,DoF,NF)
CH=KernelAbstractions.zeros(backend,FT,DoF,NF)
Vn=KernelAbstractions.zeros(backend,FT,nz,NumG,NumV+NumTr)
zPG=zeros(FT,DoF,nz)
lenghthTemp1 = 5 + 1 + 1 + 1 + NumTr + ND*(1 + 1 +  1 + NumTr)
Temp1=KernelAbstractions.zeros(backend,FT,nz,NumG,lenghthTemp1)
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
q=KernelAbstractions.zeros(backend,FT,nz,NF+NGF,2*(NumTr+1))
if ND > 0
  RhoEDMF=KernelAbstractions.zeros(backend,FT,nz,NumG,ND)
else  
  RhoEDMF=KernelAbstractions.zeros(backend,FT,nz,NumG,ND)
end  
return CacheStruct{FT,
                   typeof(KE),
                   typeof(RhoS),
                   typeof(VS)}(
  AuxG,
  KV,
  Aux2DG,
  Temp,
  KE,
  cTrS,
  Vn,
  zPG,
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
  q,
  RhoEDMF
)
end
