function fVelv(x,time::Float64,Phys,Global,Param)
  Model = Global.Model
  ProfVel = Model.ProfVel
  if ProfVel == "Const"
    vS=Param.vMax
    vS=0.0
  elseif ProfVel == "Divergent"
    Lon,Lat,R = Grids.cart2sphere(x[1],x[2],x[3])
    lonP = Lon - 2.0e0 * pi * time / Param.EndTime
    vS = 10.0 / Param.EndTime * sin(2.0 * lonP) * cos(Lat) * cos(pi * time / Param.EndTime)  
  else
    vS=0.0
  end
  return vS
end
