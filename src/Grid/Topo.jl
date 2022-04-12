function Topo(x,y,zeta,Topography)
if Topography.TopoS == "AgnesiCartX"
    h=Topography.hC/(((x-Topography.x0C)/Topography.aC)^2+1);
elseif Topography.TopoS == "AgnesiCartY"
    h=Topography.hC/(((y-Topography.y0C)/Topography.aC)^2+1);
else
    h=0;
end
Z=zeta+(Topography.H-zeta)*h/Topography.H;
dZ=1-h/Topography.H;
return (Z,dZ)
end

