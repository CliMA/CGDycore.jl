@inline function VelSphere2Cart(VelSp,lon,lat)
  sinlon = sin(lon)
  coslon = cos(lon)
  sinlat = sin(lat)
  coslat = cos(lat)
  rot11 = -sinlon
  rot21 = -sinlat * coslon
  rot31 = coslat * coslon
  rot12 = coslon
  rot22 = -sinlat * sinlon
  rot32 = coslat * sinlon
  rot13 = 0.0
  rot23 = coslat
  rot33 = sinlat
# VelCa = rot'*VelSp
  VelCa1 = rot11 * VelSp[1] + rot21 * VelSp[2] + rot31 * VelSp[3]
  VelCa2 = rot12 * VelSp[1] + rot22 * VelSp[2] + rot32 * VelSp[3]
  VelCa3 = rot13 * VelSp[1] + rot23 * VelSp[2] + rot33 * VelSp[3]
  return SVector{3}(VelCa1, VelCa2, VelCa3)
end

@inline function VelCart2Sphere(VelCa,lon,lat)
  sinlon = sin(lon)
  coslon = cos(lon)
  sinlat = sin(lat)
  coslat = cos(lat)
  rot11 = -sinlon
  rot21 = -sinlat * coslon
  rot31 = coslat * coslon
  rot12 = coslon
  rot22 = -sinlat * sinlon
  rot32 = coslat * sinlon
  rot13 = 0.0
  rot23 = coslat
  rot33 = sinlat
# VelSp = rot*VelCa
  VelSp1 = rot11 * VelCa[1] + rot12 * VelCa[2] + rot13 * VelCa[3]
  VelSp2 = rot21 * VelCa[1] + rot22 * VelCa[2] + rot23 * VelCa[3]
  VelSp3 = rot31 * VelCa[1] + rot32 * VelCa[2] + rot33 * VelCa[3]
  return SVector{3}(VelSp1, VelSp2, VelSp3)
end

@inline function MSphere2Cart(lon,lat)
  sinlon = sin(lon)
  coslon = cos(lon)
  sinlat = sin(lat)
  coslat = cos(lat)
  rot11 = -sinlon
  rot21 = -sinlat * coslon
  rot31 = coslat * coslon
  rot12 = coslon
  rot22 = -sinlat * sinlon
  rot32 = coslat * sinlon
  rot13 = 0.0
  rot23 = coslat
  rot33 = sinlat
  return  @SMatrix [rot11 rot21 rot31; rot12 rot22 rot32; rot13 rot23 rot33]
end

@inline function MCart2Sphere(lon,lat)
  sinlon = sin(lon)
  coslon = cos(lon)
  sinlat = sin(lat)
  coslat = cos(lat)
  rot11 = -sinlon
  rot21 = -sinlat * coslon
  rot31 = coslat * coslon
  rot12 = coslon
  rot22 = -sinlat * sinlon
  rot32 = coslat * sinlon
  rot13 = 0.0
  rot23 = coslat
  rot33 = sinlat
  return  @SMatrix [rot11 rot12 rot13; rot21 rot22 rot23; rot31 rot32 rot33]
end
