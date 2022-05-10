function Topo(x,y,z,zeta,Topography)
if Topography.TopoS == "AgnesiCartX"
    h=Topography.hC/(((x-Topography.x0C)/Topography.aC)^2+1);
elseif Topography.TopoS == "AgnesiCartY"
    h=Topography.hC/(((y-Topography.y0C)/Topography.aC)^2+1);
elseif Topography.TopoS == "SchaerSphereCircle"
    (Lon,Lat,R)=cart2sphere(x,y,z);
    Z=max(R-Topography.RadEarth,0);
    GreatCircleR = Topography.RadEarth * acos(sin(Topography.PertLat) * sin(Lat) +
        cos(Topography.PertLat) * cos(Lat) * cos(Lon-Topography.PertLon))
    h = Topography.h0 * exp( -(GreatCircleR  / Topography.d0)^2) * cos(pi * GreatCircleR / Topography.ksi0)^2 
elseif Topography.TopoS == "SchaerSphereRidge"
    (Lon,Lat,R)=cart2sphere(x,y,z);
    Z=max(R-Topography.RadEarth,0);
    GreatCircleR = Topography.RadEarth * (Lon - Topography.PertLon) * cos(Lat)
    h = Topography.h0 * exp( -(GreatCircleR  / Topography.d0)^2) * cos(pi * GreatCircleR / Topography.ksi0)^2 * cos(Lat)
elseif Topography.TopoS == "SchaerCart"
   h = Topography.h0 * exp( -(x  / Topography.d0)^2) * cos(pi * x / Topography.ksi0)^2
elseif Topography.TopoS == "AdvectionSchaer"
   h = Topography.h0 * cos(pi * x /Topography.lambda)^2 
   if abs(x) <= Topography.a
     h = h * cos(0.5 * pi * x / Topography.a)^2
   else
     h = 0.0
   end  
else
    h=0;
end
Z=zeta+(Topography.H-zeta)*h/Topography.H;
dZ=1-h/Topography.H;
return (Z,dZ)
end

