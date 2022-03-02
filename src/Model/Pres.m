function p=Pres(V,Param)
switch Param.Pressure
  case 'InternalEnergy'
    Rd=Param.Rd;
    Cvd=Param.Cvd;
    IEPos=Param.IEPos;
    p=Rd/Cvd*(V(IEPos,:));
  case 'ShallowWater'
    hPos=Param.hPos;
    p=0.5*Param.Grav*V(hPos,:).*V(hPos,:);
  case 'PotTemp'
    Rd=Param.Rd;
    p0=Param.p0;
    kappa=Param.kappa;
    p=(Rd*RhoTh/p0^kappa).^(1/(1-kappa));
end