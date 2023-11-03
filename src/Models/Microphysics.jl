abstract type Microphysics end

Base.@kwdef struct SimpleMicrophysics <: Microphysics end


function (::SimpleMicrophysics)(Phys,RhoPos,ThPos,RhoVPos,RhoCPos,RelCloud,Rain)
  function Microphysics(U,p)
    FT = eltype(U)
    Rho = U[RhoPos]
    RhoTh = U[ThPos]
    RhoV = max(U[RhoVPos], FT(0))
    RhoC = max(U[RhoCPos], FT(0))
    RhoD = Rho - RhoV - RhoC
    Cpml = Phys.Cpd * RhoD + Phys.Cpv * RhoV + Phys.Cpl * RhoC
    Rm = Phys.Rd * RhoD + Phys.Rv * RhoV
    kappaM = Rm / Cpml
    T = p / Rm
    p_vs = fpvs(T,Phys.T0)
    a = p_vs / (Phys.Rv * T) - RhoV
    b = RhoC
    FPh = FT(0.5) * RelCloud * (a + b - sqrt(a * a + b * b))
    L = LatHeat(T,Phys.L00,Phys.Cpl,Phys.Cpv,Phys.T0)
    FR = -FPh * Rain
    FRhoTh = RhoTh*((-L/(Cpml*T) - log(p / Phys.p0) * (Rm / Cpml) * (Phys.Rv / Rm + 
      (Phys.Cpl - Phys.Cpv) / Cpml)  + Phys.Rv / Rm) * FPh +
      (FT(1) / Rho - log(p/Phys.p0) * (Rm / Cpml) * (Phys.Cpl / Cpml)) * FR)
    return (-FR, FRhoTh, FPh, -FPh - FR)
  end
  return Microphysics
end

