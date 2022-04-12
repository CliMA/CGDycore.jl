function Pressure!(p,RhoTh,Rho,KE,Global)
    (;  Rd,
        Cvd,
        Grav,
        p0,
        kappa) = Global.Phys
    str = Global.Model.Equation
    if str == "Compressible"
        if Global.Model.Thermo == "Energy"
          p=(Rd/Cvd)*(RhoTh-Rho.*(KE+Grav*repmat(Grid.zP,1,size(Rho,1))'));
        else
          p .= p0 .* (Rd .* RhoTh ./ p0).^(1.0e0 ./ (1.0e0-kappa));
        end
    elseif str == "Shallow"
        p .= 0.5 .* Grav .* RhoTh.^2;
    end
end
