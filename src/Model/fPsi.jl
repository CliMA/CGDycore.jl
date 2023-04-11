function fPsi(x,time,Global,Param)
  Model=Global.Model
  Phys=Global.Phys
  ProfVel = Model.ProfVel
  if ProfVel == "Divergent"  
    (Lon,Lat,R) = cart2sphere(x[1],x[2],x[3])
    lonP = Lon - 2.0e0 * pi * time / Param.EndTime
    Psi = Param.FacVel * R^2 / Param.EndTime * sin(lonP)^2 * cos(Lat)^2 * 
      cos(pi * time / Param.EndTime) - 2 * pi * R^2 / Param.EndTime * sin(Lat)
  end   
  return Psi
end  
