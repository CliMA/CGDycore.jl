function fpBGrd(x,time,Global,Param)
  Model=Global.Model
  Phys=Global.Phys
  str = lowercase(Model.ProfpBGrd)
  if str == "baldaufcart"
    delta=Param.Grav/(Param.Rd*Param.T0);

    p=Param.p0*exp(-delta*x[3]);
  elseif str == "barowavesdryphere"
    (Lon,Lat,R)=cart2sphere(x[1],x[2],x[3]);
    Z=max(R-Phys.RadEarth,0);
    T0=0.5*(Param.T0E+Param.T0P);
    ConstA=1.0/Param.LapseRate;
    ConstB=(T0-Param.T0P)/(T0*Param.T0P);
    ConstC=0.5*(Param.K+2.0)*(Param.T0E-Param.T0P)/(Param.T0E*Param.T0P);
    ConstH=Phys.Rd*T0/Phys.Grav;
    ScaledZ=Z/(Param.B*ConstH);
    Tau1=ConstA*Param.LapseRate/T0*exp(Param.LapseRate/T0*Z) +
      ConstB*(1.0-2.0*ScaledZ*ScaledZ)*exp(-ScaledZ*ScaledZ);
    Tau2=ConstC*(1.0-2.0*ScaledZ*ScaledZ)*exp(-ScaledZ*ScaledZ);
    IntTau1=ConstA*(exp(Param.LapseRate/T0*Z)-1.0) +
      ConstB*Z*exp(-ScaledZ*ScaledZ);
    IntTau2=ConstC*Z*exp(-ScaledZ*ScaledZ);
    if Param.Deep
      RRatio= R/Param.EarthRadius;
    else
      RRatio = 1.0;
    end
    InteriorTerm=(RRatio*cos(Lat))^Param.K -
      Param.K/(Param.K+2.0)*(RRatio*cos(Lat))^(Param.K+2.0);
    Temperature=1.0/(RRatio*RRatio)/(Tau1-Tau2*InteriorTerm);
    p=Phys.p0*exp(-Phys.Grav/Phys.Rd *
        (IntTau1-IntTau2*InteriorTerm));    
  elseif str == "isothermal"
    (Lon,Lat,R)=cart2sphere(x[1],x[2],x[3]);
    Z=max(R-Phys.RadEarth,0);
    p = Phys.p0 * exp(-Phys.Grav * Z / (Phys.Rd * Param.TEq))
  else
    p=0;
  end
return p
end
