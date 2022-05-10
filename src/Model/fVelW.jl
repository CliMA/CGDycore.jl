function fVelW(x,time,Global)
  Model=Global.Model
  Param=Global.Model.Param
  Phys=Global.Phys
  str = lowercase(Model.ProfVelW)
  if str == "linear"
    f=x[1];
  elseif str == "advectiontestdeform"
    f = Param.uMax * sinpi(x[3] / Param.H) * cospi(time / Param.EndTime)
  elseif str == "heldsuarezcart"
    f=sin(pi*x[3]/Param.H)
  else
    f=0;
  end
  return f
end

