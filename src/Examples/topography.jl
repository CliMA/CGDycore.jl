abstract type Topography end

Base.@kwdef struct AgnesiHill{T} <: Topography 
  x0C::T = 0
  aC::T = 1000
  hC::T = 400
end

function (profile::AgnesiHill)()
  (; x0C, aC, hC) = profile
  function local_profile(x)
    FT = eltype(x)
    h = hC / (((x[1]  -x0C) / aC)^2 + FT(1));
    return h
  end  
  return local_profile
end

Base.@kwdef struct Flat <: Topography end

function (profile::Flat)()
  function local_profile(x)
    FT = eltype(x)
    h = FT(0)
  end
  return local_profile
end


