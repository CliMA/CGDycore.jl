function fTr(x,Global)
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
  else
    Tr = 0.0  
  end
  return Tr
end

