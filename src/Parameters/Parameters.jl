module Parameters

using ArgParse

include("parse_commandline.jl")
# Earth parameters
  const global ScaleFactor = 1.0
  const global RadEarth = 6.37122e+6 / ScaleFactor
  const global Grav =  9.80616
  const global Omega = 2 * pi / 24.0 / 3600.0 * ScaleFactor

# Thermodynamic parameters
  const global  Cpd = 1004.0
  const global Cvd = 717.0
  const global Cpv = 1885.0
  const global Cvv = 1424.0
  const global Cpl = 4186.0
  const global Cpi = 2110.0
  const global Rd = Cpd - Cvd
  const global Rv = Cpv - Cvv
# L00 = 2.5000e6 + (Cpl - Cpv) * 273.15
  const global L0V =  2.5000e6 # 2500800
  const global L0S =  2.834e6
  const global L0F =  L0S - L0V
  const global p0 = 1.0e5
  const global Rho0 = 1.41e0
  const global Gamma = Cpd / Cvd
  const global kappa = Rd / Cpd
  const global T0 = 273.15
  const global T00 = 273.15 -35.0
end 
