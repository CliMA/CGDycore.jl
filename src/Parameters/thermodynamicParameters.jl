struct ThermodynamicParameters{FT<:AbstractFloat}
  Cpd::FT
  Cvd::FT
  Cpv::FT
  Cvv::FT
  Cpl::FT
  Cpi::FT
  Rd::FT
  Rv::FT
  L0V::FT
  L0S::FT
  L0F::FT
  p0::FT
  Rho0::FT
  Gamma::FT
  kappa::FT
  T0::FT
  T00::FT
  cS::FT
end
function ThermodynamicParameters{FT}() where FT<:AbstractFloat
  Cpd::FT = 1004.0
  Cvd::FT = 717.0
  Cpv::FT = 1885.0
  Cvv::FT = 1424.0
  Cpl::FT = 4186.0
  Cpi::FT = 2110.0
  Rd::FT = Cpd - Cvd
  Rv::FT = Cpv - Cvv
# L00 = 2.5000e6 + (Cpl - Cpv) * 273.15
  L0V::FT =  2.5000e6 # 2500800 
  L0S::FT =  2.834e6
  L0F::FT =  L0S - L0V
  p0::FT = 1.0e5
  Rho0::FT = 1.41e0
  Gamma::FT = Cpd / Cvd
  kappa::FT = Rd / Cpd
  T0::FT = 273.15
  T00::FT = 273.15 -35.0
  cS::FT = 360

  return ThermodynamicParameters{FT}(
    Cpd,
    Cvd,
    Cpv,
    Cvv,
    Cpl,
    Cpi,
    Rd,
    Rv,
    L0V,
    L0S,
    L0F,
    p0,
    Rho0,
    Gamma,
    kappa,
    T0,
    T00,
    cS,
  )  
end  
