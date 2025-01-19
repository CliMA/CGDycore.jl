module Thermodynamics


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

#@inline function fpws(T,T0)
#  FT = eltype(T)
#  # ClaudiusClapperon
#  # Phys.p0 * (T / Phys.T0)^((Phys.Cpv - Phys.Cpl) / Phys.Rv) *
#  #   exp((Phys.L0V / Phys.Rv) *(1.0 / Phys.T0 - 1.0 / T))
#  T_C = T - T0
#  FT(611.2) * exp(FT(17.62) * T_C / (FT(243.12) + T_C))
#end

#@inline function fpis(T,T0)
#  FT = eltype(T)
#  # ClaudiusClapperon
#  # Phys.p0 * (T / Phys.T0)^((Phys.Cpv - Phys.Cpl) / Phys.Rv) *
#  #   exp((Phys.L0V / Phys.Rv) *(1.0 / Phys.T0 - 1.0 / T))
#  T_C = T - T0
#  FT(611.2) * exp(FT(21.87) * T_C / (-FT(7.66) + T))
#end

@inline function fpws(T,Phys)
  FT = eltype(T)
  # ClaudiusClapperon
  # Phys.p0 * (T / Phys.T0)^((Phys.Cpv - Phys.Cpl) / Phys.Rv) *
  #   exp((Phys.L0V / Phys.Rv) *(1.0 / Phys.T0 - 1.0 / T))
  T_C = T - Phys.T0
  FT(611.2) * exp(FT(17.62) * T_C / (FT(243.12) + T_C))
end

@inline function fpis(T,Phys)
  FT = eltype(T)
  # ClaudiusClapperon
  # Phys.p0 * (T / Phys.T0)^((Phys.Cpv - Phys.Cpl) / Phys.Rv) *
  #   exp((Phys.L0V / Phys.Rv) *(1.0 / Phys.T0 - 1.0 / T))
  T_C = T - Phys.T0
  FT(611.2) * exp(FT(21.87) * T_C / (-FT(7.66) + T))
end

@inline function dfpwsdT(T,Phys)
  FT = eltype(T)
  # ClaudiusClapperon
  # Phys.p0 * (T / Phys.T0)^((Phys.Cpv - Phys.Cpl) / Phys.Rv) *
  #   exp((Phys.L0V / Phys.Rv) *(1.0 / Phys.T0 - 1.0 / T))
  T_C = T - Phys.T0
  fpws = FT(611.2) * exp(FT(17.62) * T_C / (FT(243.12) + T_C))
  dfpwsdT = fpws * FT(17.62) * ((FT(243.12) + T_C) - T_C) / (FT(243.12) + T_C) / (FT(243.12) + T_C)
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


@inline function LatHeatV(T,Phys)
  # Vaporization
  L = Phys.L0V + (Phys.Cpv - Phys.Cpl) * (T - Phys.T0)
end


@inline function LatHeatS(T,Phys)
  # Sublimation
  L = Phys.L0S + (Phys.Cpv - Phys.Cpl) * (T - Phys.T0)
end

@inline function LatHeatF(T,Phys)
  L = Phys.L0F + (Phys.Cpl - Phys.Cpi) * (T - Phys.T0)
end

@inline function SpIntEnergyDry(T,Phys)
  Phys.Cvd * (T - Phys.T0) - Phys.Rd * Phys.T0
end  

@inline function dSpIntEnergyDrydT(T,Phys)
  Phys.Cvd 
end  

@inline function SpIntEnergyVap(T,Phys)
  Phys.Cvv * (T - Phys.T0) + Phys.L0V - Phys.Rv * Phys.T0
end  

@inline function dSpIntEnergyVapdT(T,Phys)
  Phys.Cvv
end  

@inline function SpIntEnergyLiq(T,Phys)
  Phys.Cpl * (T - Phys.T0)
end  

@inline function dSpIntEnergyLiqdT(T,Phys)
  Phys.Cpl
end  

@inline function SpIntEnergyIce(T,Phys)
  Phys.Cpi * (T - Phys.T0) - Phys.L0F
end  

@inline function dSpIntEnergyIcedT(T,Phys)
  Phys.Cpi
end  


@inline function dLatHeatdT(T,Phys)
# L = Phys.L0V - (Phys.Cpl - Phys.Cpv) * (T - Phys.T0)
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
  Lv = LatHeatV(T,Phys)
  ThE = T * ((RhoD * Phys.Rd * T) / Phys.p0)^(-Phys.Rd / (Phys.Cpd + Phys.Cpl * rt)) *
    exp(Lv * rv /((Phys.Cpd + Phys.Cpl * rt) * T))
end

@inline function FischerBurmeisterPres(Rho,RhoV,RhoC,p,Phys)
  RhoD = Rho - RhoV - RhoC
  Rm = Phys.Rd * RhoD + Phys.Rv * RhoV
  kappaM = Rm / Cpml
  T = p / Rm
  p_vs = fpws(T,Phys)
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
  p_vs = fpws(T,Phys)
  a = p_vs / (Phys.Rv * T) - RhoV
  b = RhoC
  0.5 * (a + b - sqrt(a * a + b * b))
end  

@inline function InternalEnergy(Rho,RhoT,T,Phys;RhoR=0.0)
  pWS = Thermodynamics.fpws(T,Phys)
  RhoS = pWS / (Phys.Rv * T)  
  RhoV = min(RhoS,RhoT)
  RhoC = RhoT - RhoV
  e = (Rho - RhoV -RhoC) * SpIntEnergyDry(T,Phys) + 
   RhoV * SpIntEnergyVap(T,Phys) + 
   (RhoC + RhoR) * SpIntEnergyLiq(T,Phys)   
  return e, RhoV, RhoC 
end

@inline function InternalEnergyW(Rho,RhoV,RhoC,T,Phys;RhoR=0.0)
  e = (Rho - RhoV -RhoC) * SpIntEnergyDry(T,Phys) + 
   RhoV * SpIntEnergyVap(T,Phys) + 
   (RhoC + RhoR) * SpIntEnergyLiq(T,Phys)   
  return e
end

@inline function dInternalEnergyWdT(Rho,RhoT,T,Phys;RhoR=0.0)
# e = (Rho - RhoV -RhoC) * SpIntEnergyDry(T,Phys) + 
#  RhoV * SpIntEnergyVap(T,Phys) + 
#  (RhoC + RhoR) * SpIntEnergyLiq(T,Phys)   
  pWS = Thermodynamics.fpws(T,Phys)
  RhoS = pWS / (Phys.Rv * T)  
  RhoV = min(RhoS,RhoT)
  RhoC = RhoT - RhoV
  dpWSdT = Thermodynamics.dfpwsdT(T,Phys)
  dRhoSdT = (dpWSdT * Phys.Rv * T - pWS * Phys.Rv) / (Phys.Rv * T)^2 

  e = (Rho - RhoV -RhoC) * SpIntEnergyDry(T,Phys) + 
   RhoV * SpIntEnergyVap(T,Phys) + 
   (RhoC + RhoR) * SpIntEnergyLiq(T,Phys)   
  de = (Rho - RhoV -RhoC) * dSpIntEnergyDrydT(T,Phys) + 
   RhoV * dSpIntEnergyVapdT(T,Phys) + 
   dRhoSdT * SpIntEnergyVap(T,Phys) + 
   (RhoC + RhoR) * dSpIntEnergyLiqdT(T,Phys) +    
   -dRhoSdT * SpIntEnergyLiq(T,Phys)    
  return e, de, RhoV, RhoC 
end

@inline function InternalEnergyI(Rho,RhoV,RhoC,RhoI,T,Phys;RhoR=0.0,RhoS=0.0)
  e = (Rho - RhoV -RhoC - RhoI) * SpIntEnergyDry(T,Phys) + 
   RhoV * SpIntEnergyVap(T,Phys) + 
   (RhoC + RhoR) * SpIntEnergyLiq(T,Phys) +   
   (RhoI + RhoS) * SpIntEnergyIce(T,Phys)   
end

@inline function dInternalEnergydT(Rho,RhoV,RhoC,RhoI,T,Phys;RhoR=0.0,RhoS=0.0)
  e = (Rho - RhoV -RhoC - RhoI) * SpIntEnergyDry(T,Phys) + 
   RhoV * SpIntEnergyVap(T,Phys) + 
   (RhoC + RhoR) * SpIntEnergyLiq(T,Phys) +   
   (RhoI + RhoS) * SpIntEnergyIce(T,Phys)   

  de = (Rho - RhoV -RhoC) * SpIntEnergyDry(T,Phys) + 
   RhoV * SpIntEnergyVap(T,Phys) + 
   (RhoC + RhoR) * SpIntEnergyLiq(T,Phys) +   
   (RhoI + RhoS) * SpIntEnergyIce(T,Phys)   
end
end
