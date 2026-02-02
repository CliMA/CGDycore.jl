module Parameters


include("earthParameters.jl")
# Earth parameters
  RadEarth = 0
  Grav = 0
  Omega = 0

  function SetParameters(FT;ScaleFactor=FT(1))
    EP = EarthParameters{FT}(ScaleFactor)
    global RadEarth = EP.RadEarth
    global Grav = EP.Grav
    global Omega = EP.Omega
#=
    TP = ThermodynamicParameters{FT}(ScaleFactor)
    global Cpd = TP.Cpd,
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
=#  
  end  
end
