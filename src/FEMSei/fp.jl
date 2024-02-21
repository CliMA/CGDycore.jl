function fp(x,y,z)
  lon,lat,r = Grids.cart2sphere(x,y,z)
  lat0 = 4.0*atan(1.0)
  lon0 = 2.0*atan(1.0)
  d = acos(sin(lat0)*sin(lat)+cos(lat0)*cos(lat)*cos(lon-lon0))
  if abs(d) <= 0.8
    h = cos(pi*d/0.8/2)^2
  else
    h = 0.0 
  end
end
function up(x,y,z)
  lon,lat,r = Grids.cart2sphere(x,y,z)
  VelSp = zeros(3)
  VelSp[1] = 40 * cos(lat)
  VelCa = VelSphere2Cart(VelSp,lon,lat)
  return VelCa
end
