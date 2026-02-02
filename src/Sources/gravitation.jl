abstract type GravitationType end

Base.@kwdef struct GravitationShallow <: GravitationType end

function (GravitationFun::GravitationShallow)(Phys)
  @inline function Gravitation(x,y,z)
    FT = eltype(x)
    return Phys.Grav
  end
  return Gravitation
end

Base.@kwdef struct GravitationDeep <: GravitationType end

function (GravitationFun::GravitationDeep)(Phys)
  @inline function Gravitation(x,y,z)
    FT = eltype(x)
    r = sqrt(x^2 + y^2 + z^2)
    return Phys.Grav * (Phys.RadEarth / r)^2
  end
  return Gravitation
end

Base.@kwdef struct GravitationNo <: GravitationType end

function (GravitationFun::GravitationNo)()
  @inline function Gravitation(x,y,z)
    FT = eltype(x)
    return FT(0)
  end
  return Gravitation
end

abstract type GeoPotentialType end

Base.@kwdef struct GeoPotentialShallow <: GeoPotentialType end

function (GeoPotentialFun::GeoPotentialShallow)(Phys)
  @inline function GeoPotential(x,y,z)
    FT = eltype(x)
    r = sqrt(x^2 + y^2 + z^2)
    return Phys.Grav * r 
  end
  return GeoPotential
end


Base.@kwdef struct GeoPotentialDeep <: GeoPotentialType end

function (GeoPotentialFun::GeoPotentialDeep)(Phys,::Grids.SphericalGrid)
  @inline function GeoPotential(x,y,z)
    FT = eltype(x)
    r = sqrt(x^2 + y^2 + z^2)
    return Phys.Grav * (Phys.RadEarth - Phys.RadEarth^2 / r)
#   return Phys.Grav * max(r- Phys.RadEarth,FT(0))
  end
  return GeoPotential
end

function (GeoPotentialFun::GeoPotentialDeep)(Phys,::Grids.CartesianGrid)
  @inline function GeoPotential(x,y,z)
    return Phys.Grav * z 
  end
  return GeoPotential
end

Base.@kwdef struct GeoPotentialNo <: GeoPotentialType end

function (GeoPotentialFun::GeoPotentialNo)()
  @inline function GeoPotential(x,y,z)
    FT = eltype(x)
    return FT(0)
  end
  return GeoPotential
end

abstract type BuoyancyType end

Base.@kwdef struct BuoyancyBoussinesq <: BuoyancyType end

function (::BuoyancyBoussinesq)(Param,wPos,bPos)
  @inline function BuoyancyFun(F,U,x)
    FT = eltype(F)
    F[wPos] += U[bPos]
    F[bPos] += -Param.N^2 * U[wPos]
  end
  return BuoyancyFun
end

