function p=Pressure(RhoTh,Rho,KE,Param)
switch Param.Equation
  case 'Compressible'
    if strcmp(Param.Thermo,'Energy')
      p=(Param.Rd/Param.Cvd)*(RhoTh-Rho.*(KE+Param.Grav*repmat(Param.Grid.zP,1,size(Rho,1))'));
    else
      p=Param.p0*(Param.Rd*RhoTh/Param.p0).^(1.0e0/(1.0e0-Param.kappa));
    end
  case 'Shallow'
    p=0.5*Param.Grav*RhoTh.^2;
end