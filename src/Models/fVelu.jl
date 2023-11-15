function fVelu(x,time::Float64,Phys,Global,Param)
  Model = Global.Model
  ProfVel = Model.ProfVel
  if ProfVel == "Const"
    uS=Param.uMax
  else
    uS=0.0
  end
  return uS
end
