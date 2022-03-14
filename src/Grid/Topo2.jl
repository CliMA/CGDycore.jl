function [Z,dZ]=Topo(x,zeta,Param)
switch Param.TopoS
  case 'AgnesiCart'
    h=Param.hC/(((x-Param.x0C)/Param.aC)^2+1);
  otherwise
    h=0;
end
Z=zeta+(Param.H-zeta)*h/Param.H;
dZ=1-h/Param.H;
end

