function fTBGrd(x,time,Global,Param)
  Model=Global.Model
  Phys=Global.Phys

  if lowercase(Model.ProfTheta)== "baldaufcart"
    T=Param.T0;
  end
  return T
end




