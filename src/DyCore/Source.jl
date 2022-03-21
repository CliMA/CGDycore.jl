function Source!(F,U,CG,Param)
nz=Param.Grid.nz;
if strcmp(Param.Thermo,"Energy")
else
  @views Pres=Pressure(U[:,:,Param.ThPos],U[:,:,Param.RhoPos],U[:,:,Param.RhoPos],Param);
end
sigma=Pres/Param.p0;
@views T=Pres./(Param.Rd*U[:,:,Param.RhoPos]);

matlab_max2(X) = max.(X, 0)
height_factor = matlab_max2((sigma .- Param.sigma_b) ./ (1-Param.sigma_b));

temp = (Param.T_equator .-
    Param.DeltaT_y*repmat(sin.(Param.latN).^2,1,nz) .-
    Param.DeltaTh_z*log.(sigma).*repmat(cos.(Param.latN).^2,1,nz)) .*
    sigma.^(Param.Rd/Param.Cpd)

matlab_max(X) = max.(X, Param.T_min)

@views DeltaRhoT=(Param.k_a .+ (Param.k_s - Param.k_a) .* height_factor .*
  repmat(cos.(Param.latN).^4,1,nz)) .* U[:,:,Param.RhoPos] .*
  (T .- matlab_max(temp));

@views F[:,:,Param.uPos:Param.vPos] .-= -Param.k_f*height_factor.*U[:,:,Param.uPos:Param.vPos];
if strcmp(Param.Thermo,"Energy")
else
@views  F[:,:,Param.ThPos] .-= DeltaRhoT.*(Param.p0./Pres).^Param.kappa;
end
end

