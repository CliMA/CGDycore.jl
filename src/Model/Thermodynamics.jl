#=
struct Equation :<AbstractType end
struct CompressibleTheta :<Equation end
struct CompressibleMoistTheta :<Equation end
=#

@inline function fThetaV(Rho,RhoV,RhoL,T,Phys)
  Cpml = fCpml(Rho,RhoV,RhoL,Phys)
  Rm = fRm(Rho,RhoV,RhoL,Phys)
  Pres = fPres(Rho,RhoV,RhoL,T,Phys)
  T * (Phys.p0 / p)^(Rm / Cpml)
end

@inline function fCpml(Rho,RhoV,RhoL,Phys)
  RhoD = Rho - RhoV - RhoL - RhoI
  Phys.Cpd * RhoD + Phys.Cpv * RhoV + Phys.Cpl * RhoL
end

@inline function fCvml(Rho,RhoV,RhoL,Phys)
  RhoD = Rho - RhoV - RhoL - RhoI
  Phys.Cvd * RhoD + Phys.Cvv * RhoV + Phys.Cpl * RhoL
end

@inline function fPres(Rho,RhoV,RhoL,T,Phys)
  RhoD = Rho - RhoV - RhoL 
  (Phys.Rd * RhoD + Phys.Rv * RhoV) * T 
end

@inline function fpvs(T,Phys)
  T_C = T - 273.15
  611.2 * exp(17.62 * T_C / (243.12 + T_C))
return

@inline function fpv(RhoV,T,Phys)
  RhoD = Rho - RhoV - RhoL 
  Phys.Rv * RhoV * T 
end

@inline function fRm(Rho,RhoV,RhoL,T,Phys)
  RhoD = Rho - RhoV - RhoL 
  (Phys.Rd * RhoD + Phys.Rv * RhoV) 
end

@inline function fTemp(Rho,RhoV,RhoL,RhoThetaV,,Phys)
  RhoD = Rho - RhoV - RhoL 
  Cpml = fCpml(Rho,RhoV,RhoL,Phys)
  Rm = fRm(Rho,RhoV,RhoL,Phys)
  Kappa = Rm / Cpml
  (Phys.Rd * RhoThetaV / Phys.p0^Kappa)^(1.0 / (1.0 - Kappa)) / Rm

end

@inline function LatHeat(T,Phys)
  L = Phys.L00 - (Phys.Cpl - Phys.Cpv) * T
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



