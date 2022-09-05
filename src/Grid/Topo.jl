function Topo(x,y,z,zeta,Topography,zs)
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
elseif Topography.TopoS == "EarthOrography"   
   h = zs
else
    h=0;
end
Z=zeta+(Topography.H-zeta)*h/Topography.H;
dZ=1-h/Topography.H;
return (Z,dZ)
end

function EarthTopography()
  # Load ETOPO1 ice-sheet surface data
  # Ocean values are considered 0
  data = NCDataset("/Users/knoth/Documents/GitHub/CliMA/CGDycore.jl/ETOPO1_Ice_g_gdal.grd")
  # Unpack information
  x_range = data["x_range"][:]
  y_range = data["y_range"][:]
  z_range = data["z_range"][:]
  spacing = data["spacing"][:]
  dimension = data["dimension"][:]
  elevation = data["z"][:]
  lon = collect(x_range[1]:spacing[1]:x_range[2])
  lat = collect(y_range[1]:spacing[2]:y_range[2])
  nlon = dimension[1]
  nlat = dimension[2]
  zlevels = reshape(elevation, (nlon, nlat))

  map_source = zlevels
  map_source[map_source .< 0.0] .= 0.0
  spline_2d=Spline2D(lon, lat, reverse(map_source, dims=2);kx=3, ky=3, s=0.0)
  h = evaluate(spline_2d,11.9,51.0)
  stop


  #(w,xw)=CGDycore.GaussLobattoQuad(OrdPoly)
  #NF = Grid.NumFaces
  #zs = zeros(OrdPoly+1,OrdPoly+1,NF)
  #for iF = 1 : NF
  #  F = Grid.Faces[iF]
  #  for j = 1 : OrdPoly+1
  #    for i = 1 : OrdPoly+1
  #      X1=0.25*((1-xw[i])*(1-xw[j])*F.P[1].x+
  #       (1+xw[i])*(1-xw[j])*F.P[2].x+
  #       (1+xw[i])*(1+xw[j])*F.P[3].x+
  #       (1-xw[i])*(1+xw[j])*F.P[4].x);
  #      X2=0.25*((1-xw[i])*(1-xw[j])*F.P[1].y+
  #       (1+xw[i])*(1-xw[j])*F.P[2].y+
  #       (1+xw[i])*(1+xw[j])*F.P[3].y+
  #       (1-xw[i])*(1+xw[j])*F.P[4].y);
  #      X3=0.25*((1-xw[i])*(1-xw[j])*F.P[1].z+
  #       (1+xw[i])*(1-xw[j])*F.P[2].z+
  #       (1+xw[i])*(1+xw[j])*F.P[3].z+
  #       (1-xw[i])*(1+xw[j])*F.P[4].z);
  #    (lam, phi)  = CGDycore.cart2sphere(X1,X2,X3)
  #    lamD = lam * 180.0 / pi
  #    phiD = phi * 180.0 / pi
  #    zs[i,j,iF] = evaluate(spline_2d,lamD,phiD)
  #  end
  #end
  return spline_2d
end  

function EarthTopography(OrdPoly,Grid)
  # Load ETOPO1 ice-sheet surface data
  # Ocean values are considered 0
  data = NCDataset("/Users/knoth/Documents/GitHub/CliMA/CGDycore.jl/ETOPO1_Ice_g_gdal.grd")
  # Unpack information
  x_range = data["x_range"][:]
  y_range = data["y_range"][:]
  z_range = data["z_range"][:]
  spacing = data["spacing"][:]
  dimension = data["dimension"][:]
  elevation = data["z"][:]
  lon = collect(x_range[1]:spacing[1]:x_range[2])
  lat = collect(y_range[1]:spacing[2]:y_range[2])
  nlon = dimension[1]
  nlat = dimension[2]
  zlevels = reshape(elevation, (nlon, nlat))

  map_source = zlevels
  map_source[map_source .< 0.0] .= 0.0
  spline_2d=Spline2D(lon, lat, reverse(map_source, dims=2);kx=3, ky=3, s=0.0)

  (w,xw)=GaussLobattoQuad(OrdPoly)
  NF = Grid.NumFaces
  zs = zeros(OrdPoly+1,OrdPoly+1,NF)
  for iF = 1 : NF
    F = Grid.Faces[iF]
    for j = 1 : OrdPoly+1
      for i = 1 : OrdPoly+1
        X1=0.25*((1-xw[i])*(1-xw[j])*F.P[1].x+
         (1+xw[i])*(1-xw[j])*F.P[2].x+
         (1+xw[i])*(1+xw[j])*F.P[3].x+
         (1-xw[i])*(1+xw[j])*F.P[4].x);
        X2=0.25*((1-xw[i])*(1-xw[j])*F.P[1].y+
         (1+xw[i])*(1-xw[j])*F.P[2].y+
         (1+xw[i])*(1+xw[j])*F.P[3].y+
         (1-xw[i])*(1+xw[j])*F.P[4].y);
        X3=0.25*((1-xw[i])*(1-xw[j])*F.P[1].z+
         (1+xw[i])*(1-xw[j])*F.P[2].z+
         (1+xw[i])*(1+xw[j])*F.P[3].z+
         (1-xw[i])*(1+xw[j])*F.P[4].z);
      (lam, phi)  = CGDycore.cart2sphere(X1,X2,X3)
      lamD = lam * 180.0 / pi - 180.0
      phiD = phi * 180.0 / pi 
      zs[i,j,iF] = evaluate(spline_2d,lamD,phiD)
    end
  end
end
return zs
end  

