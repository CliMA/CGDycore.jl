mutable struct PhysParameters
  RadEarth::Float64 
  Grav::Float64 
  Cpd::Float64
  Cvd::Float64
  Cpv::Float64
  Cvv::Float64
  Cpl::Float64
  Rd::Float64
  Rv::Float64
  L00::Float64
  p0::Float64
  Gamma::Float64
  kappa::Float64
  Omega::Float64
end
function PhysParameters()
  RadEarth = 6.37122e+6
  Grav = 9.81e0
  Cpd=1004.0e0
  Cvd=717.0e0
  Cpv=1885.0e0
  Cvv=1424.0e0
  Cpl=4186.0e0
  Rd=Cpd-Cvd
  Rv=Cpv-Cvv
  L00 = 2.5000e6 + (Cpl - Cpv) * 273.15
  p0=1.0e5
  Gamma=Cpd/Cvd
  kappa=Rd/Cpd
  Omega=2*pi/24.0/3600.0
 return PhysParameters(
  RadEarth,
  Grav,
  Cpd,
  Cvd,
  Cpv,
  Cvv,
  Cpl,
  Rd,
  Rv,
  L00,
  p0,
  Gamma,
  kappa,
  Omega,
  )
end 
