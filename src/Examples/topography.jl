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

Base.@kwdef struct SchaerHill{T} <: Topography 
  d0::T = 5000.0
  ksi0::T = 4000.0
  h0::T = 250.0
end

function (profile::SchaerHill)()
  (; d0, ksi0, h0) = profile
  function local_profile(x)
    FT = eltype(x)
    h = h0 * exp( -(x[1]  / d0)^2) * cos(pi * x[1] / ksi0)^2
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

Base.@kwdef struct BaroWaveHill{T} <: Topography 
  h0::T = 2.e3
  lon1::T = (72 + 10) * pi / 180
  lon2::T = (140 + 10) * pi / 180
  lonBar::T = 7 * pi / 180
  c::T = 0.5 * lonBar*(-log(0.1))^(-1/2)
  lat1::T = pi / 4 
  lat2::T = pi / 4
  latBar::T = 40 * pi / 180
  d::T = 0.5 * latBar*(-log(0.1))^(-1/6)
end  

function (profile::BaroWaveHill)()
  (; h0, lon1, lon2, lonBar, c, lat1, lat2, latBar, d) = profile
  function local_profile(x)
    FT = eltype(x)
    (lon,lat,r) = Grids.cart2sphere(x[1],x[2],x[3])
    d1 = lon - lon1
    d2 = lon - lon2
    l1 = min(d1, FT(2*pi) - d1)
    l2 = min(d2, FT(2*pi) - d2)
    h = h0 * (exp(-((lat-lat1) / d )^6 - (l1 / c)^2) + 
      exp(-((lat-lat2) / d)^6 - (l2 / c)^2)) 
    return h
  end  
  return local_profile
end

Base.@kwdef struct SchaerSphereCircle{T} <: Topography
  PertLat::T = 0.0
  PertLon::T = pi / 4
  d0::T = 5000.0
  ksi0::T = 4000.0
  h0::T = 250.0
end

function (profile::SchaerSphereCircle)(Param,Phys)
  (; PertLat, PertLon, d0, ksi0, h0) = profile
  function local_profile(x)
    FT = eltype(x)
    (lon,lat,r) = Grids.cart2sphere(x[1],x[2],x[3])
    rad = Phys.RadEarth / Param.X
    z = max(FT(0), r - rad)
    GreatCircleR = rad * acos(sin(PertLat) * sin(lat) +
        cos(PertLat) * cos(lat) * cos(lon - PertLon))
    h = h0 * exp( -(GreatCircleR  / d0)^2) * cos(pi * GreatCircleR / ksi0)^2
    return h
  end  
  return local_profile
end



