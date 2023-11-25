#=
struct Equation :<AbstractType end
struct CompressibleTheta :<Equation end
struct CompressibleMoistTheta :<Equation end
=#

@inline function fThetaV(Rho,RhoV,RhoL,T,Phys)
  Cpml = fCpml(Rho,RhoV,RhoL,Phys)
  Rm = fRm(Rho,RhoV,RhoL,Phys)
  p = Rm * T
  ThetaV = T * (Phys.p0 / p)^(Rm / Cpml) * Rm / (Phys.Rd * Rho)
  return ThetaV
end

@inline function fCpml(Rho,RhoV,RhoL,Phys)
  RhoD = Rho - RhoV - RhoL 
  Phys.Cpd * RhoD + Phys.Cpv * RhoV + Phys.Cpl * RhoL
end

@inline function fCvml(Rho,RhoV,RhoL,Phys)
  RhoD = Rho - RhoV - RhoL
  Phys.Cvd * RhoD + Phys.Cvv * RhoV + Phys.Cpl * RhoL
end

@inline function fPres(Rho,RhoV,RhoL,T,Phys)
  RhoD = Rho - RhoV - RhoL 
  (Phys.Rd * RhoD + Phys.Rv * RhoV) * T 
end

@inline function fpvs(T,T0)
  FT = eltype(T)
  # ClaudiusClapperon
  # Phys.p0 * (T / Phys.T0)^((Phys.Cpv - Phys.Cpl) / Phys.Rv) *
  #   exp((Phys.L00 / Phys.Rv) *(1.0 / Phys.T0 - 1.0 / T))
  T_C = T - T0
  FT(611.2) * exp(FT(17.62) * T_C / (FT(243.12) + T_C))
end

@inline function fpvs1(T,Phys)
  FT = eltype(T)
  # ClaudiusClapperon
  # Phys.p0 * (T / Phys.T0)^((Phys.Cpv - Phys.Cpl) / Phys.Rv) *
  #   exp((Phys.L00 / Phys.Rv) *(1.0 / Phys.T0 - 1.0 / T))
  T_C = T - Phys.T0
  FT(611.2) * exp(FT(17.62) * T_C / (FT(243.12) + T_C))
end

@inline function dfpvs1dT(T,Phys)
  FT = eltype(T)
  # ClaudiusClapperon
  # Phys.p0 * (T / Phys.T0)^((Phys.Cpv - Phys.Cpl) / Phys.Rv) *
  #   exp((Phys.L00 / Phys.Rv) *(1.0 / Phys.T0 - 1.0 / T))
  T_C = T - Phys.T0
  fpvs = FT(611.2) * exp(FT(17.62) * T_C / (FT(243.12) + T_C))
  dfpvsdT = fpvs * FT(17.62) * ((FT(243.12) + T_C) - T_C) / (FT(243.12) + T_C) / (FT(243.12) + T_C)
end


@inline function fpv(RhoV,T,Phys)
  RhoD = Rho - RhoV - RhoL 
  Phys.Rv * RhoV * T 
end

@inline function fRm(Rho,RhoV,RhoL,Phys)
  RhoD = Rho - RhoV - RhoL 
  (Phys.Rd * RhoD + Phys.Rv * RhoV) 
end

@inline function fTemp(Rho,RhoV,RhoL,RhoThetaV,Phys)
  RhoD = Rho - RhoV - RhoL 
  Cpml = fCpml(Rho,RhoV,RhoL,Phys)
  Rm = fRm(Rho,RhoV,RhoL,Phys)
  Kappa = Rm / Cpml
  (Phys.Rd * RhoThetaV / Phys.p0^Kappa)^(1.0 / (1.0 - Kappa)) / Rm
end

@inline function LatHeat(T,L00,Cpl,Cpv,T0)
  L = L00 - (Cpl - Cpv) * (T - T0)
end

@inline function LatHeat(T,Phys)
  L = Phys.L00 - (Phys.Cpl - Phys.Cpv) * (T - Phys.T0)
end

@inline function dLatHeatdT(T,Phys)
# L = Phys.L00 - (Phys.Cpl - Phys.Cpv) * (T - Phys.T0)
  dLdT = -(Phys.Cpl - Phys.Cpv)
end

@inline function fThE(Rho,RhoV,RhoL,RhoThetaV,Phys)
  RhoD = Rho - RhoV - RhoL 
  Cpml = fCpml(Rho,RhoV,RhoL,Phys)
  Rm = fRm(Rho,RhoV,RhoL,Phys)
  Kappa = Rm / Cpml
  rv = RhoV  / RhoD
  rt = (RhoV + RhoL) / RhoD
  T = (Phys.Rd * RhoThetaV / Phys.p0^Kappa)^(1.0 / (1.0 - Kappa)) / Rm
  Lv = LatHeat(T,Phys)
  ThE = T * ((RhoD * Phys.Rd * T) / Phys.p0)^(-Phys.Rd / (Phys.Cpd + Phys.Cpl * rt)) *
    exp(Lv * rv /((Phys.Cpd + Phys.Cpl * rt) * T))
end

@inline function FischerBurmeisterPres(Rho,RhoV,RhoC,p,Phys)
  RhoD = Rho - RhoV - RhoC
  Rm = Phys.Rd * RhoD + Phys.Rv * RhoV
  kappaM = Rm / Cpml
  T = p / Rm
  p_vs = fpvs(T,Phys)
  a = p_vs / (Rv * T) - RhoV
  b = RhoC
  0.5 * (a + b - sqrt(a * a + b * b))
end  

@inline function FischerBurmeisterRhoTh(Rho,RhoV,RhoC,RhoThetaV,Phys)
  RhoD = Rho - RhoV - RhoC
  Cpml = Phys.Cpd * RhoD + Phys.Cpv * RhoV + Phys.Cpl * RhoC
  Rm = (Phys.Rd * RhoD + Phys.Rv * RhoV) 
  Kappa = Rm / Cpml
  T = (Phys.Rd * RhoThetaV / Phys.p0^Kappa)^(1.0 / (1.0 - Kappa)) / Rm
  p_vs = fpvs(T,Phys)
  a = p_vs / (Phys.Rv * T) - RhoV
  b = RhoC
  0.5 * (a + b - sqrt(a * a + b * b))
end  

@inline function InternalEnergy(Rho,RhoV,RhoC,T,Phys)
  Lv = LatHeat(T,Phys)
  e = (Phys.Cvd * (Rho - RhoV -RhoC) + Phys.Cvv * RhoV + Phys.Cpl * RhoC) * T - RhoC * Lv
end

@inline function dInternalEnergydT(Rho,RhoV,RhoC,T,Phys)
  Lv = LatHeat(T,Phys)
  dLvdT = dLatHeatdT(T,Phys)
  pVS = fpvs1(T,Phys)
  dpVSdT = dfpvs1dT(T,Phys)
# RhoVS = pVS / (Phys.Rv *  T)
  dRhoVSdT = (dpVSdT * T - pVS) / (T * T * Phys.Rv)
# e = (Phys.Cvd * (Rho - RhoV -RhoC) + Phys.Cvv * RhoV + Phys.Cpl * RhoC) * T - RhoC * Lv
  dedT = (Phys.Cvd * (Rho - RhoV -RhoC) + Phys.Cvv * RhoV + Phys.Cpl * RhoC)  - RhoV * dLvdT +
    ((Phys.Cvv - Phys.Cpl) * T + Lv) * dRhoVSdT
end

