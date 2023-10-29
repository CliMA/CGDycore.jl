function Topo(x,y,z,zeta,Topography,zs)
  if Topography.TopoS == "AgnesiCartX"
    x0C = Topography.P1
    aC= Topography.P2
    hC = Topography.P3
    h=hC/(((x-x0C)/aC)^2+1);
  elseif Topography.TopoS == "AgnesiCartY"
    y0C = Topography.P1
    aC= Topography.P2
    hC = Topography.P3
    h=hC/(((y-y0C)/aC)^2+1);
  elseif Topography.TopoS == "AgnesiCartXY"
    x0C = Topography.P1
    y0C = Topography.P2
    aC= Topography.P3
    hC = Topography.P4
    h=hC/(((x-x0C)/aC)^2+((y-y0C)/aC)^2+1);
  elseif Topography.TopoS == "GaussCartX"
    x0C = Topography.P1
    aC= Topography.P2
    hC = Topography.P3
    h=hC*exp(-((x-x0C)/aC)^2)
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
   d0 = Topography.P1
   ksi0 = Topography.P2
   h0 = Topography.P3
   h = h0 * exp( -(x  / d0)^2) * cos(pi * x / ksi0)^2
elseif Topography.TopoS == "AdvectionSchaer"
   h = Topography.h0 * cos(pi * x /Topography.lambda)^2 
   if abs(x) <= Topography.a
     h = h * cos(0.5 * pi * x / Topography.a)^2
   else
     h = 0.0
   end  
elseif Topography.TopoS == "EarthOrography"   
   h = zs
else
    h=0;
end
Z=zeta+(Topography.H-zeta)*h/Topography.H;
dZ=1-h/Topography.H;
return (Z,dZ)
end

