@inline function SaturationAdjustmentIEW(Rho,RhoIE,RhoT,T,Phys)

#
# Find RhoV, RhoC and T  such that with
# IE = (Cvd * (Rho - RhoV -RhoC) + Cvv * RhoV + Cpl * RhoC) * T - RhoC * Lv
# RhoT = RhoV + RhoC
# RhoV <= RhoVS
# RhoC >= 0
# (RhoVS - RhoV) * RhoC = 0

  FT = eltype(Rho)
  T_min = FT(100)
  T_max = FT(360)
  sqeps = sqrt(eps(FT))
# Newton-Method
  RhoV = FT(0)
  RhoC = FT(0)
  T0 = T
  e, RhoV, RhoC = InternalEnergyW(T,Rho,RhoT,Phys)
  f = e - RhoIE
  TR = T_max
  TL = T_min
  fL = FT(1)
  fR = FT(-1)
  right = false
  left = false
  if f > FT(0)
    fR = f
    TR = T
    right = true
  else
    fL = f
    TL = T
    left = true
  end  
  TP = T * (FT(1) + sqeps)
  eP, = InternalEnergyW(TP,Rho,RhoT,Phys)
  dfdT = (eP - e) / (TP - T)
  TNew = min(max(T - f/dfdT,T_min),T_max)
  rightF = false
  leftF = false
  for iTer = 1 : 100
    eNew, RhoVNew, RhoCNew = InternalEnergyW(TNew,Rho,RhoT,Phys)
    fNew = eNew - RhoIE
    if abs(fNew) < FT(1.e-3)
      T = TNew  
      RhoV = RhoVNew
      RhoC = RhoCNew
      break
    end
    if fNew > FT(0)
      fR = fNew  
      TR = TNew
      if rightF
        fL = FT(0.5) * fL
      end  
      right = true
      rightF = true
      leftF = false
    else  
      fL = fNew
      TL = TNew
      if leftF
        fR = FT(0.5) * fR
      end  
      left = true
      rightF = false
      leftF = true
    end  
    if left && right
      TNew = (-fR*TL+fL*TR) / (fL-fR) 
    else
      dfdT = (fNew - f) / (TNew - T)
      f = fNew
      T = TNew
      TNew = min(max(T - f/dfdT,T_min),T_max)
    end    
    if iTer > 40
      @show iTer
      @show fL,TL
      @show fR,TR
      stop
    end  
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

@inline function InternalEnergyW(T,Rho,RhoT,Phys)
  pWS = Thermodynamics.fpws(T,Phys)
  RhoS = pWS / (Phys.Rv * T) 
  RhoV = min(RhoS,RhoT)
  RhoC = RhoT - RhoV
  e = Thermodynamics.InternalEnergy(Rho,RhoV,RhoC,T,Phys)
  return e, RhoV, RhoC
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


