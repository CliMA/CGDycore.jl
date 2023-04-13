function fPsi(x,y,z,time,Global,Param)
  Model=Global.Model
  Phys=Global.Phys
  ProfVel = Model.ProfVel
  Psi = 0.0
  EndTime =  1.0368e6
  FacVel = 10.0
  if ProfVel == "Divergent"  
    (Lon,Lat,R) = cart2sphere(x,y,z)
    lonP = Lon - 2.0 * pi * time / EndTime
    Psi = FacVel * R^2 / EndTime * sin(lonP)^2 * cos(Lat)^2 * 
      cos(pi * time / EndTime) - 2 * pi * R^2 / EndTime * sin(Lat)
  end   
  return Psi
end  
