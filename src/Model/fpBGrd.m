function p = fpBGrd(x,Param)
switch lower(Param.ProfpBGrd)
  case 'baldaufcart'
    delta=Param.Grav/(Param.Rd*Param.T0);
    
    p=Param.p0*exp(-delta*x(3));
  otherwise
    p=0;
end
end