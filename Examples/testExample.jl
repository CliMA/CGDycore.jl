Base.@kwdef struct ParamGalewskySphere
  H0G = 10000.0
  hH = 120.0
  alphaG = 1.0/3.0
  betaG = 1.0/15.0
  lat0G = pi/7.0
  lat1G = pi/2.0-lat0G
  eN = exp(-4.0/(lat1G-lat0G)^2.0)
  uM = 80.0
  Omega = 2*pi/24.0/3600.0
end

Param = ParamGalewskySphere()

abstract type Example end

Base.@kwdef struct Example1 <: Example end

mutable struct Model
  Eddy::Any
  Prof::Any
end



function (profile::Example1)(Param)
  function local_profile(x)
    FT = eltype(x)
    A = 2 * x
    return A
  end
  function eddy(x)
    FT = eltype(x)
    A = 3 * x
    return A
  end
  return local_profile,eddy
end

Prof, Eddy = Example1()(Param)
@show eltype(Prof)
@show eltype(Eddy)
Model1 = Model(Prof,Eddy)

Prof(5)
Eddy(7)

