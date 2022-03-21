function Pressure(RhoTh,Rho,KE,Param)
    (;Equation,
        Thermo,
        Rd,
        Cvd,
        Grav,
        Grid,
        p0,
        kappa) = Param
    str = Equation
    if str == "Compressible"
        if strcmp(Thermo,"Energy")
          p=(Rd/Cvd)*(RhoTh-Rho.*(KE+Grav*repmat(Grid.zP,1,size(Rho,1))'));
        else
          p=p0*(Rd*RhoTh/p0).^(1.0e0/(1.0e0-kappa));
        end
    elseif str == "Shallow"
        p=0.5*Grav*RhoTh.^2;
    end
    return p
end