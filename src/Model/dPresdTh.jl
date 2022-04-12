function dPresdTh!(p,RhoTh,Param)
  @. p=Param.Phys.Rd*(Param.Phys.Rd*RhoTh/Param.Phys.p0)^(Param.Phys.kappa/(1-Param.Phys.kappa));
end
