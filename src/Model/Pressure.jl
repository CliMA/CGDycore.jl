function Pressure(U,KE,zP,Global)
  Rd = Global.Phys.Rd
  Cvd = Global.Phys.Cvd
  Grav = Global.Phys.Grav
  if Global.Model.Equation == "Compressible"
    p=(Rd/Cvd)*(U[5]-U[1]*(KE+Grav*zP))
  elseif Global.Model.Equation == "CompressibleMoist"
  end
end 

function Pressure!(p,RhoTh,Rho,Tr,KE,zP,Global)
  (; Rd,
     Cvd,
     Cpd,
     Rv,
     Cvv,
     Cpv,
     Cpl,
     p0,
     Grav,
     kappa) = Global.Phys

  
  Equation = Global.Model.Equation
  if Equation == "Compressible"
     if Global.Model.Thermo == "TotalEnergy"
       @inbounds for i in eachindex(p)  
         p[i] = (Rd / Cvd) * (RhoTh[i] - Rho[i] * (KE[i] + Grav * zP[i]))
       end  
     elseif Global.Model.Thermo == "InternalEnergy"
       @inbounds for i in eachindex(p)  
         p[i] = (Rd / Cvd) * RhoTh[i] 
       end  
     else
       @inbounds for i in eachindex(p)  
         p[i] = p0 * (Rd * RhoTh[i] / p0)^(1.0 / (1.0 - kappa));
       end  
    end
  elseif Equation == "CompressibleMoist"
    @views TrRhoV = Tr[size(p)...,Global.Model.RhoVPos]
    @views TrRhoC = Tr[size(p)...,Global.Model.RhoCPos]
    @inbounds for i in eachindex(p)  
      RhoV = TrRhoV[i]
      RhoC = TrRhoC[i]
      RhoD = Rho[i] - RhoV - RhoC
      Cpml = Cpd * RhoD + Cpv * RhoV + Cpl * RhoC
      Rm  = Rd * RhoD + Rv * RhoV
      kappaM = Rm / Cpml
      p[i] = (Rd * RhoTh[i] / p0^kappaM)^(1.0 / (1.0 - kappaM))
    end  
  elseif Equation == "Shallow"
      p = 0.5 * Grav * RhoTh^2;
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
    if Global.Model.Thermo == "TotalEnergy" || Global.Model.Thermo == "InternalEnergy"
      @. dpdTh = Rd / Cvd  
    else  
      @. dpdTh=Rd*(Rd*RhoTh/p0)^(kappa/(1.0-kappa));
    end  
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
