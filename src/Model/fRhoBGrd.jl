function Rho = fRhoBGrd(x,Param)
switch lower(Param.ProfRhoBGrd)
  case 'baldaufcart'
    delta=Param.Grav/(Param.Rd*Param.T0);
    p=Param.p0*exp(-delta*x(3));
    TLoc=Param.T0;
    Rho=p/(Param.Rd*TLoc);
  otherwise
    Rho=0;
end
end








