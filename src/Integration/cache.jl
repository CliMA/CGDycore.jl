mutable struct CacheStruct{FT<:AbstractFloat,
                           AT3<:AbstractArray,
                           AT4<:AbstractArray,
                           AT5<:AbstractArray}
  KV::AT3
  KE::AT3
  Thermo::AT4
  Vn::AT4
  fV::AT4
  Temp1::AT4
  fRhoS::AT4
  RhoS::AT4
  q::AT4
  RhoEDMF::AT4
  k::AT5
  fS::AT5
  VS::AT5
  f::AT5
  zPG::Array{FT, 2}
  Aux2DG::Array{FT, 3}
  Temp::Array{FT, 3}
  cTrS::Array{FT, 3}
  R::Array{FT, 3}
  dZ::Array{FT, 3}
  Ymyn::Array{FT, 4}
  Y::Array{FT, 4}
  Z::Array{FT, 4}
end

function CacheStruct{FT}(backend) where FT<:AbstractFloat
  #AT3
  KV = KernelAbstractions.zeros(backend,FT,0,0,0)
  KE = KernelAbstractions.zeros(backend,FT,0,0,0)
  #AT4
  Thermo = KernelAbstractions.zeros(backend,FT,0,0,0,0)
  Vn = KernelAbstractions.zeros(backend,FT,0,0,0,0)
  fV = KernelAbstractions.zeros(backend,FT,0,0,0,0)
  Temp1 = KernelAbstractions.zeros(backend,FT,0,0,0,0)
  fRhoS = KernelAbstractions.zeros(backend,FT,0,0,0,0)
  RhoS = KernelAbstractions.zeros(backend,FT,0,0,0,0)
  q = KernelAbstractions.zeros(backend,FT,0,0,0,0)
  RhoEDMF = KernelAbstractions.zeros(backend,FT,0,0,0,0)
  #AT4
  k = KernelAbstractions.zeros(backend,FT,0,0,0,0,0)
  fS = KernelAbstractions.zeros(backend,FT,0,0,0,0,0)
  VS = KernelAbstractions.zeros(backend,FT,0,0,0,0,0)
  f = KernelAbstractions.zeros(backend,FT,0,0,0,0,0)
  # CPU Arrays
  # Dim2  
  zPG = zeros(FT,0,0)
  # Dim3  
  Aux2DG = zeros(FT,0,0,0)
  Temp=zeros(FT,0,0,0)
  cTrS=zeros(FT,0,0,0)
  R=zeros(FT,0,0,0)
  dZ=zeros(FT,0,0,0)
  # Dim4  
  Ymyn=zeros(FT,0,0,0,0)
  Y=zeros(FT,0,0,0,0)
  Z=zeros(FT,0,0,0,0)

  return CacheStruct{FT,
                   typeof(KV),
                   typeof(Thermo),
                   typeof(k)}(
  # AT3                   
    KV,
    KE,
  # AT4                   
    Thermo,
    Vn,
    fV,
    Temp1,
    fRhoS,
    RhoS,
    q,
    RhoEDMF,
  # AT5                   
    k,
    fS,
    VS,
    f,
  # CPU Arrays
  # Dim2  
    zPG,
  # Dim3  
    Aux2DG,
    Temp,
    cTrS,
    R,
    dZ,
  # Dim4  
    Ymyn,
    Y,
    Z,

  )
end

function CacheStruct{FT}(backend,DoF,NF,NGF,NumG,M,nz,NumV,NumTr,ND,NumThermo) where FT<:AbstractFloat
  #AT3
  KV = KernelAbstractions.zeros(backend,FT,M,nz,NumG)
  KE = KernelAbstractions.zeros(backend,FT,DoF,M,nz)
  #AT4
  Thermo = KernelAbstractions.zeros(backend,FT,M,nz,NumG,NumThermo)
  @views @. Thermo[:,:,:,2] = FT(250.0)
  Vn = KernelAbstractions.zeros(backend,FT,0,0,0,0)
  fV = KernelAbstractions.zeros(backend,FT,0,0,0,0)
  lenghthTemp1 = 5 + 1 + 1 + 1 + NumTr + ND*(1 + 1 +  1 + NumTr)
  Temp1 = KernelAbstractions.zeros(backend,FT,M,nz,NumG,lenghthTemp1)
  fRhoS = KernelAbstractions.zeros(backend,FT,0,0,0,0)
  RhoS = KernelAbstractions.zeros(backend,FT,0,0,0,0)
  q = KernelAbstractions.zeros(backend,FT,M,nz,NF+NGF,2*(NumTr+1))
  RhoEDMF = KernelAbstractions.zeros(backend,FT,0,0,0,0)
  if ND > 0
    RhoEDMF = KernelAbstractions.zeros(backend,FT,M,nz,NumG,ND)
  else
    RhoEDMF = KernelAbstractions.zeros(backend,FT,0,0,0,0)
  end
  #AT5
  k = KernelAbstractions.zeros(backend,FT,0,0,0,0,0)
  fS = KernelAbstractions.zeros(backend,FT,0,0,0,0,0)
  VS = KernelAbstractions.zeros(backend,FT,0,0,0,0,0)
  f = KernelAbstractions.zeros(backend,FT,0,0,0,0,0)
  # CPU Arrays
  # Dim2  
  zPG = zeros(FT,DoF,nz)
  # Dim3  
  Aux2DG = zeros(FT,1,NumG,NumTr+1)
  Temp = zeros(FT,DoF,nz,NF)
  cTrS = zeros(FT,DoF,NF,NumTr)
  R=zeros(FT,0,0,0)
  dZ=zeros(FT,0,0,0)
  # Dim4  
  Ymyn=zeros(FT,0,0,0,0)
  Y=zeros(FT,0,0,0,0)
  Z=zeros(FT,0,0,0,0)

  return CacheStruct{FT,
                   typeof(KV),
                   typeof(Thermo),
                   typeof(k)}(
  # AT3
    KV,
    KE,
  # AT4
    Thermo,
    Vn,
    fV,
    Temp1,
    fRhoS,
    RhoS,
    q,
    RhoEDMF,
  # AT5
    k,
    fS,
    VS,
    f,
  # CPU Arrays
  # Dim2
    zPG,
  # Dim3
    Aux2DG,
    Temp,
    cTrS,
    R,
    dZ,
  # Dim4
    Ymyn,
    Y,
    Z,
  )
end

mutable struct CacheDGStruct{FT<:AbstractFloat,
                           AT4<:AbstractArray,
                           AT5<:AbstractArray}
  Vn::AT4
  k::AT5
  fV::AT4
end

function CacheDGStruct{FT}(backend) where FT<:AbstractFloat
  Vn=KernelAbstractions.zeros(backend,FT,0,0,0,0)
  k=KernelAbstractions.zeros(backend,FT,0,0,0,0,0)
  fV=KernelAbstractions.zeros(backend,FT,0,0,0,0)
return CacheDGStruct{FT,
                   typeof(Vn),
                   typeof(k)}(
  Vn,
  k,
  fV,
)
end

