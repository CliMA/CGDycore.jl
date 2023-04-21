function fPsi(x,time::Float64,Global,Param)
  Model=Global.Model
  Phys=Global.Phys
  ProfVel = Model.ProfVel
  Psi = 0.0
  if ProfVel == "Divergent"  
    (Lon,Lat,R) = cart2sphere(x[1],x[2],x[3])
    lonP = Lon - 2.0 * pi * time / Param.EndTime
    Psi = 10.0  / Param.EndTime * sin(lonP)^2 * cos(Lat)^2 * 
      cos(pi * time / Param.EndTime) - 2 * pi  / Param.EndTime * sin(Lat)
  end   
  return Psi
end  
