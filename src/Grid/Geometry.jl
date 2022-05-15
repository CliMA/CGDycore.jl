function GreatCircle(Lon1,Lat1,Lon2,Lat2) 
  return acos(sin(Lat1) * sin(Lat2) +
        cos(Lat1) * cos(Lat2) * cos(Lon2-Lon1))
end  

