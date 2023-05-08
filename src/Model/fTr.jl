function fTr(x,time,Global,Param)
  Model=Global.Model
  Phys=Global.Phys
  str = lowercase(Model.ProfTr)
  if str == "cylinder"
    if abs(x[1] - Param.xC)<Param.xH && abs(x[3] - Param.zC) < Param.zH
      Tr = 1.0  
    else
      Tr = 0.0  
    end  
  elseif str == "advectionspheredcmipq1"
    (Lon,Lat,R) = cart2sphere(x[1],x[2],x[3])
    Z=max(R-Phys.RadEarth,0)
    zd = Z - Param.z_c
    # great circle distances
    rd1 = Phys.RadEarth * GreatCircle(Param.Lon_c1,Param.Lat_c,Lon,Lat)
    rd2 = Phys.RadEarth * GreatCircle(Param.Lon_c2,Param.Lat_c,Lon,Lat)
    d1 = min(1.0, (rd1 / Param.R_t)^2 + (zd / Param.Z_t)^2)
    d2 = min(1.0, (rd2 / Param.R_t)^2 + (zd / Param.Z_t)^2)
    Tr = 0.5 * (1 + cos(pi * d1)) + 0.5 * (1 + cos(pi * d2))
  elseif str == "advectionspheredcmipq2"
    (Lon,Lat,R) = cart2sphere(x[1],x[2],x[3])
    Z=max(R-Phys.RadEarth,0)
    zd = Z - Param.z_c
    # great circle distances
    rd1 = Phys.RadEarth * GreatCircle(Param.Lon_c1,Param.Lat_c,Lon,Lat)
    rd2 = Phys.RadEarth * GreatCircle(Param.Lon_c2,Param.Lat_c,Lon,Lat)
    d1 = min(1.0, (rd1 / Param.R_t)^2 + (zd / Param.Z_t)^2)
    d2 = min(1.0, (rd2 / Param.R_t)^2 + (zd / Param.Z_t)^2)
    q1 = 0.5 * (1 + cos(pi * d1)) + 0.5 * (1 + cos(pi * d2))  
    Tr = 0.9 - 0.8 * q1^2
  elseif str == "advectionschaer"
    r = sqrt( ((x[1] - Param.xC)/Param.Ax)^2 + ((x[3] - Param.zC)/Param.Az)^2)
    if r <= 1.0
      Tr = Param.q0*cos(0.5 * pi *r)^2  
    else    
      Tr = 0.0 
    end 
  elseif str == "advectiontestdeform"
    r1 = sqrt( (x[1] - Param.xB1)^2 + (x[2] - Param.yB1)^2 + (x[3] - Param.zB1)^2)
    r2 = sqrt( (x[1] - Param.xB2)^2 + (x[2] - Param.yB2)^2 + (x[3] - Param.zB2)^2)
    if r1 < Param.r0
      Tr = 0.1 + 0.9 * (1 / 2) * (1 + cospi(r1 / Param.r0))
    elseif r2 < Param.r0
      Tr = 0.1 + 0.9 * (1 / 2) * (1 + cospi(r2 / Param.r0))
    else
      Tr = 0.0
    end
#   Tr = 0.95 * (exp(-5.0 * (r1 / Param.r0)^2) + exp(-5.0 * (r2 / Param.r0)^2))
  elseif str == "advectionspheregaussian"
    lon1 = Param.lon1
    lat1 = Param.lat1
    lon2 = Param.lon2
    lat2 = Param.lat2
    hMax = Param.hMax
    b = Param.b
    (Lon,Lat,R)  =  cart2sphere(x[1],x[2],x[3])
    xLoc = x[1]/R
    yLoc = x[2]/R
    zLoc = x[3]/R
    x1 = cos(lat1)*cos(lon1)
    y1 = cos(lat1)*sin(lon1)
    z1 = sin(lat1)
    x2 = cos(lat2)*cos(lon2)
    y2 = cos(lat2)*sin(lon2)
    z2 = sin(lat2)
    Tr = hMax*exp(-b*((xLoc-x1)^2+(yLoc-y1)^2+(zLoc-z1)^2)) +
        hMax*exp(-b*((xLoc-x2)^2+(yLoc-y2)^2+(zLoc-z2)^2))
  elseif str == "advectionsphereslottedcylinder"
    (lon,lat,R)  =  cart2sphere(x[1],x[2],x[3])
    lon1 = Param.lon1
    lat1 = Param.lat1
    lon2 = Param.lon2
    lat2 = Param.lat2
    R = 1.0
    r = 0.5 * R
    r1 = R * GreatCircle(lon,lat,lon1,lat1)
    r2 = R * GreatCircle(lon,lat,lon2,lat2)
    if r1 <= r && abs(lon - lon1) >= r / (6.0 * R) 
      Tr = 1.0
    elseif r2 <= r && abs(lon - lon2) >= r / (6.0 * R)
      Tr = 1.0
    elseif r1 <= r && abs(lon - lon1) < r / (6.0 * R) && lat - lat1 < -5.0 / 12.0 * r / R
      Tr = 1.0
    elseif r2 <= r && abs(lon - lon2) < r / (6.0 * R) && lat - lat2 > 5.0 / 12.0 * r / R
      Tr = 1.0
    else
      Tr = .1
    end
  elseif str == "advectioncubecart"
    if x[1] >= Param.x1 && x[1] <= Param.x2 && x[2] >= Param.y1 && x[2] <= Param.y2
      Tr = 1
    else
      Tr = 0
    end
  else
    Tr  =  0.0  
  end
  return Tr
end

