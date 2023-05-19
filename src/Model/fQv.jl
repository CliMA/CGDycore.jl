function fQv(x,time,Global,Param,Profile)
  Model=Global.Model
  Phys=Global.Phys
  str = lowercase(Model.Problem)
  if str == "barowavemoistsphere" || str == "heldsuarezmoistsphere" || str == "heldsuarezmoistsphereoro"
    (Lon,Lat,R)=cart2sphere(x[1],x[2],x[3])  
    q_0 = Param.q_0                # Maximum specific humidity (default: 0.018)
    q_t = Param.q_t                # Specific humidity above artificial tropopause
    lat_w = 2.0*pi / 9.0
    p_w = 34.0e3
    eta_crit = p_w / Phys.p0
    Z=max(R-Phys.RadEarth,0)
    T0=0.5*(Param.T0E+Param.T0P)
    ConstA=1.0/Param.LapseRate
    ConstB=(T0-Param.T0P)/(T0*Param.T0P)
    ConstC=0.5*(Param.K+2.0)*(Param.T0E-Param.T0P)/(Param.T0E*Param.T0P)
    ConstH=Phys.Rd*T0/Phys.Grav
    ScaledZ=Z/(Param.B*ConstH)
    Tau1=ConstA*Param.LapseRate/T0*exp(Param.LapseRate/T0*Z)+
    ConstB*(1.0-2.0*ScaledZ*ScaledZ)*exp(-ScaledZ*ScaledZ)
    Tau2=ConstC*(1.0-2.0*ScaledZ*ScaledZ)*exp(-ScaledZ*ScaledZ)
    IntTau1=ConstA*(exp(Param.LapseRate/T0*Z)-1.0)+
      ConstB*Z*exp(-ScaledZ*ScaledZ)
    IntTau2=ConstC*Z*exp(-ScaledZ*ScaledZ)
    if Model.Deep
      RRatio= R/Param.EarthRadius
    else
      RRatio = 1.0
    end
    InteriorTerm=(RRatio*cos(Lat))^Param.K -
      Param.K/(Param.K+2.0)*(RRatio*cos(Lat))^(Param.K+2.0)
    Temperature=1.0/(RRatio*RRatio)/(Tau1-Tau2*InteriorTerm)
    Pressure=Phys.p0*exp(-Phys.Grav/Phys.Rd *
      (IntTau1-IntTau2*InteriorTerm))
    eta = Pressure/Phys.p0
    if eta>eta_crit
      qv = q_0*exp(-(Lat/lat_w)^4)*exp(-((eta-1.)*Phys.p0/p_w)^2)
    else
      qv = q_t
    end
  end
  return qv
end    
