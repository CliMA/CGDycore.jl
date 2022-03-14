function F=Source(U,CG,Param)
F=zeros(size(U));
nz=Param.Grid.nz;
if strcmp(Param.Thermo,'Energy')
else
  Pres=Pressure(U(:,:,Param.ThPos),U(:,:,Param.RhoPos),U(:,:,Param.RhoPos),Param);
end
sigma=Pres/Param.p0;
T=Pres./(Param.Rd*U(:,:,Param.RhoPos));
height_factor=max(0,(sigma-Param.sigma_b)/(1-Param.sigma_b));

DeltaRhoT=(Param.k_a+(Param.k_s-Param.k_a)*height_factor...
   .*repmat(cos(Param.latN).^4,1,nz)).*U(:,:,Param.RhoPos)...
   .*(T-...
    max(Param.T_min...
    ,(Param.T_equator...
    -Param.DeltaT_y*repmat(sin(Param.latN).^2,1,nz)...
    -Param.DeltaTh_z*log(sigma).*repmat(cos(Param.latN).^2,1,nz))...
   .*sigma.^(Param.Rd/Param.Cpd)));

F(:,:,Param.uPos:Param.vPos)=-Param.k_f*height_factor.*U(:,:,Param.uPos:Param.vPos);
if strcmp(Param.Thermo,'Energy')
else
  F(:,:,Param.ThPos)=-DeltaRhoT.*(Param.p0./Pres).^Param.kappa;
end
end