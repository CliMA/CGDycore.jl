function fVelu(x,time::Float64,Global,Param)
  Model = Global.Model
  Phys = Global.Phys
  ProfVel = Model.ProfVel
  if ProfVel == "Const"
    uS=Param.uMax
  else
    uS=0.0
  end
  return uS
end
