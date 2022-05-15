function dPresdTh!(p,RhoTh,Global)
  (; Rd,
     Cvd,
     Cpd,
     Rv,
     Cvv,
     Cpv,
     Cpl,
     p0,
     kappa) = Global.Phys

  Equation = Global.Model.Equation
  if Equation == "Compressible"
    @. p=Param.Phys.Rd*(Param.Phys.Rd*RhoTh/Param.Phys.p0)^(Param.Phys.kappa/(1-Param.Phys.kappa));
  elseif Equation == "CompressibleMoist"
    @views @. p = dPressureMoistdTh(RhoTh,Rho,Tr[:,Global.Model.RhoVPos],
      Tr[:,Global.Model.RhoCPos],Rd,Cpd,Rv,Cpv,Cpl,p0)
  end  
end


function dPressureMoistdTh(RhoTh,Rho,RhoV,RhoC,Rd,Cpd,Rv,Cpv,Cpl,p0)
  RhoD = Rho - RhoV - RhoC
  Cpml = Cpd * RhoD + Cpv * RhoV + Cpl * RhoC
  Rm  = Rd * RhoD + Rv * RhoV
  kappaM = Rm / Cpml
  dpdTh=Rd*(Rd*RhoTh/p0)^(kappaM/(1-kappaM));
end
