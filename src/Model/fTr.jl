function fTr(x,time,Global)
  Model=Global.Model
  Param=Global.Model.Param
  Phys=Global.Phys
  str = lowercase(Model.ProfTr)
  if str == "cylinder"
    if abs(x[1] - Param.xC)<Param.xH && abs(x[3] - Param.zC) < Param.zH
      Tr = 1.0  
    else
      Tr = 0.0  
    end  
  elseif str == "advectionschaer"
    r = sqrt( ((x[1] - Param.xC)/Param.Ax)^2 + ((x[3] - Param.zC)/Param.Az)^2)
    if r <= 1.0
      Tr = Param.q0*cos(0.5 * pi *r)^2  
    else    
      Tr = 0.0 
    end 
  elseif str == "advectiontestdeform"
    r1 = sqrt( (x[1] - Param.xB1)^2 + (x[2] - Param.yB1)^2 + (x[3] - Param.zB1)^2)
    r2 = sqrt( (x[1] - Param.xB2)^2 + (x[2] - Param.yB2)^2 + (x[3] - Param.zB2)^2)
    if r1 < Param.r0
      Tr = 0.1 + 0.9 * (1 / 2) * (1 + cospi(r1 / Param.r0))
    elseif r2 < Param.r0
      Tr = 0.1 + 0.9 * (1 / 2) * (1 + cospi(r2 / Param.r0))
    else
      Tr = 0.0
    end
#   Tr = 0.95 * (exp(-5.0 * (r1 / Param.r0)^2) + exp(-5.0 * (r2 / Param.r0)^2))
  else
    Tr = 0.0  
  end
  return Tr
end

