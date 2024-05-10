abstract type SurfaceStress end

Base.@kwdef struct DivergentSphereExample <: SurfaceStress end

function (profile::DivergentSphereExample)(Param,Phys)
  function local_profile(x,time)
    FT = eltype(x)
    Rho = FT(1)
    Lon,Lat,R = Gri

