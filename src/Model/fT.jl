function T = fT(x,Param)
switch lower(Param.ProfTheta)
  case 'baldaufsphere'
    [lon,lat,r]=cart2sphere(x(1),x(2),x(3));
    r=r-Param.RadEarth;
    p=Param.p0*exp(-Param.Grav*r/(Param.Rd*Param.T0));
    Rho=p/(Param.Rd*Param.T0);
    T=p/(Rho*Param.Rd);
    d=acos(sin(Param.lat0)*sin(lat)+cos(Param.lat0)*cos(lat)*cos(lon-Param.lon0));
    T=T+Param.DeltaT*exp(-Param.ExpDist*d)*sin(pi*r/Param.H);
  case 'baldaufcart'
    delta=Param.Grav/(Param.Rd*Param.T0);
    
    dT=Param.DeltaT*exp(-(x(1)-Param.xc)^2/Param.d^2)*sin(pi*x(3)/Param.H);
    
    T=Param.T0+exp(delta/2*x(3))*dT;  
end
end