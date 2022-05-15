function Pressure!(p,RhoTh,Rho,Tr,Global)
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
    if Global.Model.Thermo == "Energy"
      p=(Rd/Cvd)*(RhoTh-Rho.*(KE+Grav*repmat(Grid.zP,1,size(Rho,1))'));
    else
      @. p = p0 * (Rd * RhoTh / p0)^(1.0e0 / (1.0e0-kappa));
    end
  elseif Equation == "CompressibleMoist"
    @views @. p = PressureMoist(RhoTh,Rho,Tr[:,Global.Model.RhoVPos],
      Tr[:,Global.Model.RhoCPos],Rd,Cpd,Rv,Cpv,Cpl,p0)
  elseif Equation == "Shallow"
      p .= 0.5 .* Grav .* RhoTh.^2;
  end
end

function Temperature!(T,RhoTh,Rho,Tr,Global)
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
    if Global.Model.Thermo == "Energy"
      T=(Rd/Cvd)*(RhoTh-Rho.*(KE+Grav*repmat(Grid.zP,1,size(Rho,1))'));
    else
      @. T = p0 * (Rd * RhoTh / p0)^(1.0e0 / (1.0e0-kappa));
    end
  elseif Equation == "CompressibleMoist"
    @views @. T = PressureMoist(RhoTh,Rho,Tr[:,Global.Model.RhoVPos],
      Tr[:,Global.Model.RhoCPos],Rd,Cpd,Rv,Cpv,Cpl,p0) / (Rho * Rd + Tr[:,Global.Model.RhoVPos] * Rv)
  end
end

function PressureMoist(RhoTh,Rho,RhoV,RhoC,Rd,Cpd,Rv,Cpv,Cpl,p0)
  RhoD = Rho - RhoV - RhoC
  Cpml = Cpd * RhoD + Cpv * RhoV + Cpl * RhoC
  Rm  = Rd * RhoD + Rv * RhoV
  kappaM = Rm / Cpml
  p = (Rd * RhoTh / p0^kappaM)^(1.0 / (1.0 - kappaM))
end


function dPresdTh!(dpdTh,RhoTh,Rho,Tr,Global)
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
    @. dpdTh=Rd*(Rd*RhoTh/p0)^(kappa/(1.0-kappa));
  elseif Equation == "CompressibleMoist"
    @views @. dpdTh = dPressureMoistdTh(RhoTh,Rho,Tr[:,Global.Model.RhoVPos],
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
