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

function (CoriolisFun::CoriolisDeep)(Phys)
  @inline function Coriolis(x,y,z,u,v,w1,w2)
    FT = eltype(x)
    sinlat = z / sqrt(x^2 + y^2 + z^2)
    coslat = sqrt(FT(1) - sinlat^2)
    Ws = FT(2) * Phys.Omega * sinlat
    Wc = FT(2) * Phys.Omega * coslat
    Fu = v * Ws - FT(0.5) * Wc * (w1 + w2) 
    Fv = -u * Ws
    Fw = Wc * u
    return Fu,Fv,Fw
  end
  return Coriolis
end

Base.@kwdef struct CoriolisDeepDG <: CoriolisType end

function (CoriolisFun::CoriolisDeepDG)(uPos,vPos,wPos)
  @inline function Coriolis(F,U,X)
    FT = eltype(X)
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



