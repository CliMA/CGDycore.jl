function Topo(x,y,zeta,Param)
if Param.TopoS == "AgnesiCartX"
    h=Param.hC/(((x-Param.x0C)/Param.aC)^2+1);
elseif Param.TopoS == "AgnesiCartY"
    h=Param.hC/(((y-Param.y0C)/Param.aC)^2+1);
else
    h=0;
end
Z=zeta+(Param.H-zeta)*h/Param.H;
dZ=1-h/Param.H;
return (Z,dZ)
end

