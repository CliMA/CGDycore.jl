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

@inline function fpvs(T,Phys)
  # ClaudiusClapperon
  # Phys.p0 * (T / Phys.T0)^((Phys.Cpv - Phys.Cpl) / Phys.Rv) *
  #   exp((Phys.L00 / Phys.Rv) *(1.0 / Phys.T0 - 1.0 / T))
  T_C = T - Phys.T0
  611.2 * exp(17.62 * T_C / (243.12 + T_C))
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

@inline function LatHeat(T,Phys)
  L = Phys.L00 - (Phys.Cpl - Phys.Cpv) * (T - Phys.T0)
end

function SaturationAdjustment!(RhoThetaV,Rho,RhoV,RhoC,Phys)

  # Compute T
  T = fTemp(Rho,RhoV,RhoC,RhoThetaV,Phys)

  # Determine DeltaRhoV
  # where RhoVStar = RhoV - DeltaRhoV
  # where RhoCStar = RhoV + DeltaRhoV
  # where TStar = T + Lv * DeltaRhoV
  # with
  # A := pvs(TStar) - pV(RhoVStar,TStar) >= 0 
  # and
  # B := RhoCStar >= 0
  # and
  # A * B = 0
  # (a + b - sqrt(a * a + b * b))
  A = fpvs(T) - fpv(RhoV,T)
  B = RhoC
  F = A + B - sqrt(A * A + B * B)


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
