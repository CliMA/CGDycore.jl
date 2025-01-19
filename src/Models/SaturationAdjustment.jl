@inline function SaturationAdjustmentIEW(Rho,RhoIE,RhoT,T,Phys)
#
# Find RhoV, RhoC and T  such that with
# IE = (Cvd * (Rho - RhoV -RhoC) + Cvv * RhoV + Cpl * RhoC) * T - RhoC * Lv
# RhoT = RhoV + RhoC
# RhoV <= RhoVS
# RhoC >= 0
# (RhoVS - RhoV) * RhoC = 0
  RhoV = eltype(Rho)(0)
  RhoC = eltype(Rho)(0)
  for i = 1 : 10
    e, de, RhoV, RhoC = Thermodynamics.dInternalEnergyWdT(T,Rho,RhoT,Phys)
    TNew = T - (e - RhoIE) / de
    if abs(T - TNew) <= eltype(Rho)(1.f-3)
      break
    end  
    T = TNew
  end
  return (RhoV, RhoC, T)
end

@inline function SaturationAdjustmentIEI(Rho,RhoIE,RhoT,T,Phys)

#
# Find RhoV, RhoC and T  such that with
# IE = (Cvd * (Rho - RhoV -RhoC) + Cvv * RhoV + Cpl * RhoC) * T - RhoC * Lv
# RhoT = RhoV + RhoC
# RhoV <= RhoVS
# RhoC >= 0
# (RhoVS - RhoV) * RhoC = 0


# Newton-Method
  RhoV = eltype(Rho)(0)
  RhoC = eltype(Rho)(0)
  RhoI = eltype(Rho)(0)
  for i = 1 : 10
    e, de, RhoV, RhoC, RhoI = Thermodynamics.dInternalEnergyIdT(T,Rho,RhoT,Phys)
    TNew = T - (e - RhoIE) / de
    if abs(T - TNew) <= eltype(Rho)(1.f-3)
      break
    end  
    T = TNew
  end
  return (RhoV, RhoC, RhoI, T)
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


