function [Z,dZ]=Topo(x,y,zeta,Param)
switch Param.TopoS
  case 'AgnesiCartX'
    h=Param.hC/(((x-Param.x0C)/Param.aC)^2+1);
  case 'AgnesiCartY'
    h=Param.hC/(((y-Param.y0C)/Param.aC)^2+1);  
  otherwise
    h=0;
end
Z=zeta+(Param.H-zeta)*h/Param.H;
dZ=1-h/Param.H;
end

