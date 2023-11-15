function fVelv(x,time::Float64,Phys,Global,Param)
  Model = Global.Model
  ProfVel = Model.ProfVel
  if ProfVel == "Const"
    vS=Param.vMax
    vS=0.0
  else
    vS=0.0
  end
  return vS
end
