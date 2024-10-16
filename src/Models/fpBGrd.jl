function fpBGrd(x,time,Phys,Global,Param,Profile)
  Model=Global.Model
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
  elseif str == "isothermalcart"
    p = Phys.p0 * exp(-Phys.Grav * x[3] / (Phys.Rd * Param.TEq))
  elseif str == "isothermalsphere"
    (Lon,Lat,R)=cart2sphere(x[1],x[2],x[3]);
    Z=max(R-Phys.RadEarth,0);
    p = Phys.p0 * exp(-Phys.Grav * Z / (Phys.Rd * Param.TEq))
  elseif str == "gravityhill" || str == "schaercart" || str == "agnesicart"
    z=x[3];
    NBr=Param.NBr;
    Grav=Phys.Grav;
    p0=Phys.p0;
    Cpd=Phys.Cpd;
    Rd=Phys.Rd;
    kappa=Phys.kappa;
    Th0=Param.Th0;
    S=NBr*NBr/Grav;
    ThLoc=Th0*exp(z*S);
    p =p0*(1.0-Grav/(Cpd*Th0*S)*(1.0-exp(-S*z)))^(Cpd/Rd);
  else
    p=0;
  end
return p
end
