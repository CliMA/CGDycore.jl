function f=hS(x,Param)
f=0;
switch Param.hS
  case 'AgnesiAnnulus'
    f=Param.h/(((x-Param.phi0)/Param.aphi)^2+1);
  case 'AgnesiCart'
    f=Param.hC/(((x-Param.x0C)/Param.aC)^2+1);
  otherwise
    f=0;
end
end