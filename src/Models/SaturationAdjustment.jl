@inline function SaturationAdjustmentIEW(Rho,RhoIE,RhoT,T,Phys)
#
# Find RhoV, RhoC and T  such that with
# IE = (Cvd * (Rho - RhoV -RhoC) + Cvv * RhoV + Cpl * RhoC) * T - RhoC * Lv
# RhoT = RhoV + RhoC
# RhoV <= RhoVS
# RhoC >= 0
# (RhoVS - RhoV) * RhoC = 0
  RhoV = 0.0
  RhoC = 0.0
  for i = 1 : 10
    e, de, RhoV, RhoC = Thermodynamics.dInternalEnergyWdT(T,Rho,RhoT,Phys)
    TNew = T - (e - RhoIE) / de
    if abs(T - TNew) <= 1.f-3
      break
    end  
    T = TNew
  end
  return (RhoV, RhoC, T)
end

function SaturationAdjustmentIEI(Rho,RhoIE,RhoT,T,Phys)

#
# Find RhoV, RhoC and T  such that with
# IE = (Cvd * (Rho - RhoV -RhoC) + Cvv * RhoV + Cpl * RhoC) * T - RhoC * Lv
# RhoT = RhoV + RhoC
# RhoV <= RhoVS
# RhoC >= 0
# (RhoVS - RhoV) * RhoC = 0


# Newton-Method
  RhoV = 0.0
  RhoC = 0.0
  RhoI = 0.0
  for iTer = 1 : 10
    e, RhoV, RhoC, RhoI = InternalEnergyI(T,Rho,RhoT,Phys)
    f = 1.0/e - 1.0/RhoIE

    TP = T * (1.0 + 1.e-8)
    eP, = InternalEnergyI(TP,Rho,RhoT,Phys)
    dfdT = (1.0/eP - 1.0/e) / (TP - T)

    TNew = T - f/dfdT
    @show f,T
    if abs(TNew - T) < 1.e-3
      break
    end
    T = TNew
  end
  return (RhoV, RhoC, RhoI, T)
end

@inline function Partition(T,Phys)
  T0 = Phys.T0
  T00 = Phys.T00

  lambda = max(min((T - T00) / (T0 - T00) , 1.0), 0.0)
end


@inline function InternalEnergyI(T,Rho,RhoT,Phys)
  pWS = Thermodynamics.fpws(T,Phys)
  pIS = Thermodynamics.fpis(T,Phys)
  RhoWS = pWS / (Phys.Rv * T) 
  RhoIS = pIS / (Phys.Rv * T) 
  lambda = Partition(T,Phys)
  RhoS = lambda * RhoWS + (1.0 - lambda) * RhoIS
  RhoV = min(RhoS,RhoT)
  RhoC = (RhoT - RhoV) * lambda
  RhoI = (RhoT - RhoV) * (1.0 - lambda)
  e = Thermodynamics.InternalEnergy(Rho,RhoV,RhoC,RhoI,T,Phys)
  return e, RhoV, RhoC, RhoI
end  


