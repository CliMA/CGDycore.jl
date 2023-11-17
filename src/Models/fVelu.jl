function fVelu(x,time::Float64,Phys,Global,Param)
  Model = Global.Model
  ProfVel = Model.ProfVel
  if ProfVel == "Const"
    uS=Param.uMax
  elseif ProfVel == "Divergent"
    Lon,Lat,R = Grids.cart2sphere(x[1],x[2],x[3])
    lonP = Lon - 2.0e0 * pi * time / Param.EndTime
    uS = 10.0 / Param.EndTime * sin(lonP)^2 *
      sin(2.0 *Lat) * cos(pi * time / Param.EndTime) + 2.0e0 * pi / Param.EndTime * cos(Lat)
  else
    uS=0.0
  end
  return uS
end
