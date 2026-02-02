struct EarthParameters{FT<:AbstractFloat}
  RadEarth::FT
  Grav::FT
  Omega::FT
end
function EarthParameters{FT}(ScaleFactor) where FT<:AbstractFloat
  RadEarth::FT = 6.37122e+6
  Grav::FT =  9.80616
  Omega::FT = 2 * pi / 24.0 / 3600.0 * ScaleFactor

  return EarthParameters{FT}(
    RadEarth,
    Grav,
    Omega,
  )  
end  
