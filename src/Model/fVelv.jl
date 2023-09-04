function fVelv(x,time::Float64,Global,Param)
  Model = Global.Model
  Phys = Global.Phys
  ProfVel = Model.ProfVel
  if ProfVel == "Const"
    vS=Param.vMax
  else
    vS=0.0
  end
  return vS
end
