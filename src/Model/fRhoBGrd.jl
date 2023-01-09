function fRhoBGrd(x,time,Global,Param)
  Model=Global.Model
  Phys=Global.Phys
  str = lowercase(Model.ProfRhoBGrd)
  if str == "baldaufcart"
    delta=Param.Grav/(Param.Rd*Param.T0);
    p=Param.p0*exp(-delta*x[3]);
    TLoc=Param.T0;
    Rho=p/(Param.Rd*TLoc);
  elseif str == "isothermal"
    (Lon,Lat,R)=cart2sphere(x[1],x[2],x[3]);
    Z=max(R-Phys.RadEarth,0);
    pLoc = Phys.p0 * exp(-Phys.Grav * Z / (Phys.Rd * Param.TEq))
    Rho = pLoc / (Phys.Rd * Param.TEq)
  else
    Rho=0;
  end
  return Rho
end
