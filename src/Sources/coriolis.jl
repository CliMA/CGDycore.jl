abstract type CoriolisType end

Base.@kwdef struct FPlane <: CoriolisType end

function (CoriolisFun::FPlane)(Param,Phys)
  @inline function Coriolis(x,y,z,u,v,w1,w2)
    FT = eltype(x)
    W = Param.f0 + Param.beta0 * (y - Param.y0)
    Fu = v * W
    Fv = - u * W
    return Fu,Fv,FT(0)
  end
  return Coriolis
end  

Base.@kwdef struct CoriolisShallow <: CoriolisType end

function (CoriolisFun::CoriolisShallow)(Phys)
  @inline function Coriolis(x,y,z,u,v,w1,w2)
    FT = eltype(x)
    sinlat = z / sqrt(x^2 + y^2 + z^2)
    W = FT(2) * Phys.Omega * sinlat 
    Fu = v * W
    Fv = -u * W
    return Fu,Fv,FT(0)
  end
  return Coriolis
end

Base.@kwdef struct CoriolisShallowDG <: CoriolisType end

function (CoriolisFun::CoriolisShallowDG)(Phys)
  @inline function Coriolis(x,y,z,u,v,w)
    FT = eltype(x)
    sinlat = z / sqrt(x^2 + y^2 + z^2)
    W = FT(2) * Phys.Omega * sinlat
    Fu = v * W
    Fv = -u * W
    return Fu,Fv,FT(0)
  end
  return Coriolis
end

Base.@kwdef struct CoriolisDeep <: CoriolisType end

function (CoriolisFun::CoriolisDeep)(uPos,vPos,wPos,::Examples.VelocityS)
  @inline function Coriolis(F,U,lon,lat)
    FT = eltype(F)
    sinlat = sin(lat)
    coslat = cos(lat)
    Ws = FT(2) * P.Omega * sinlat
    Wc = FT(2) * P.Omega * coslat
    F[uPos] += U[vPos] * Ws - FT(0.5) * Wc * U[wPos]
    F[vPos] += -U[uPos] * Ws
    F[wPos] += Wc * U[uPos]
  end
  return Coriolis
end

function (CoriolisFun::CoriolisDeep)(uPos,vPos,wPos,::Examples.VelocityC)
  @inline function Coriolis(F,U,lon,lat)
    FT = eltype(F)
    F[uPos] += FT(2) * P.Omega * U[vPos]
    F[vPos] += -FT(2) * P.Omega * U[uPos]
  end
  return Coriolis
end

Base.@kwdef struct CoriolisNo <: CoriolisType end

function (CoriolisFun::CoriolisNo)()
  @inline function Coriolis(x,y,z,u,v,w1,w2)
    FT = eltype(x)
    return FT(0),FT(0),FT(0)
  end
  return Coriolis
end  

function Coriolis!(Cor,F,U,Aux,FE,Metric,NumberThreadGPU)
  backend = get_backend(F)
  M = size(F,1)
  Nz = size(F,2)
  NumG = size(F,3)
  NumV = size(F,4)
  NumGG = min(div(NumberThreadGPU,Nz*M),NumG)
  group = (M,Nz,NumGG)
  ndrange = (M,Nz,NumG)
  KCoriolisKernel! = CoriolisKernel!(backend,group)
  KCoriolisKernel!(Cor,F,U,Metric.xS,Val(NumV);ndrange=ndrange)
end

@kernel inbounds = true function CoriolisKernel!(Cor,F,@Const(U),@Const(xS),
  ::Val{NUMV}) where {NUMV}

  _,_,iD  = @index(Local, NTuple)
  K,Iz,ID = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  FLoc = @private eltype(F) (NUMV,)
  ULoc = @private eltype(F) (NUMV,)
  xSLoc = @private eltype(F) (2,)

  if ID <= ND
    xSLoc[1] = xS[1,ID]
    xSLoc[2] = xS[2,ID]
    @unroll for iv = 1 : NUMV
      ULoc[iv] = U[K,Iz,ID,iv]
      FLoc[iv] = eltype(F)(0)
    end  
    Cor(FLoc,ULoc,xSLoc[1],xSLoc[2])
    @unroll for iv = 1 : NUMV
      F[K,Iz,ID,iv] += FLoc[iv]
    end  
  end
end

