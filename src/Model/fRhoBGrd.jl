function fRhoBGrd(x,time,Global,Param)
  Model=Global.Model
  Phys=Global.Phys
  str = lowercase(Model.ProfRhoBGrd)
  if str == "baldaufcart"
    delta=Param.Grav/(Param.Rd*Param.T0);
    p=Param.p0*exp(-delta*x[3]);
    TLoc=Param.T0;
    Rho=p/(Param.Rd*TLoc);
  elseif str == "isothermalsphere"
    (Lon,Lat,R)=cart2sphere(x[1],x[2],x[3]);
    Z=max(R-Phys.RadEarth,0);
    pLoc = Phys.p0 * exp(-Phys.Grav * Z / (Phys.Rd * Param.TEq))
    Rho = pLoc / (Phys.Rd * Param.TEq)
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
    pLoc=p0*(1-Grav/(Cpd*Th0*S)*(1-exp(-S*z))).^(Cpd/Rd);
    Rho=pLoc./((pLoc/p0).^kappa*Rd.*ThLoc);  
  else
    Rho=0;
  end
  return Rho
end
