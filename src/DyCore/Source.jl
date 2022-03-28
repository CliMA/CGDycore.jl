function Source!(F,U,CG,Param)
    str = lower(Param.ProfRho)
    Rho=view(U,:,:,Param.RhoPos)
    Th=view(U,:,:,Param.ThPos)
    uPos=Param.uPos
    vPos=Param.vPos
    ThPos=Param.ThPos
    Sigma=Param.Cache1
    height_factor=Param.Cache2
    if str == "heldsuarezsphere"
       nz=Param.Grid.nz;
       p=Pressure(Th,Rho,Th,Param);
       sigma = p / Param.p0
       height_factor = max.(0, (sigma .- Param.sigma_b) ./ (1 .- Param.sigma_b))
       ΔρT =
        (Param.k_a .+ (Param.k_s .- Param.k_a) .* height_factor .*
        repmat(cos.(Param.latN).^4,1,nz)) .*
        Rho .*
        ( # ᶜT - ᶜT_equil
            p ./ (Rho .* Param.Rd) .- max.(
                Param.T_min,
                (Param.T_equator .- Param.DeltaT_y .* repmat(sin.(Param.latN).^2,1,nz) .- 
                Param.DeltaTh_z .* log.(sigma) .* repmat(cos.(Param.latN).^2,1,nz)) .*
                sigma.^(Param.Rd / Param.Cpd),
            )
        )
        @views  F[:,:,uPos:vPos] .-= (Param.k_f * height_factor) .* U[:,:,uPos:vPos]
        @views  F[:,:,ThPos]  .-= ΔρT .*  (Param.p0 ./ p).^Param.kappa

 #      nz=Param.Grid.nz;
 #      if strcmp(Param.Thermo,"Energy")
 #      else
 #        @views Pres=Pressure(U[:,:,Param.ThPos],U[:,:,Param.RhoPos],U[:,:,Param.RhoPos],Param);
 #      end
 #      sigma=Pres/Param.p0;
 #      @views T=Pres./(Param.Rd*U[:,:,Param.RhoPos]);

 #      matlab_max2(X) = max.(X, 0)
 #      height_factor = matlab_max2((sigma .- Param.sigma_b) ./ (1-Param.sigma_b));

 #      temp = (Param.T_equator .-
 #          Param.DeltaT_y*repmat(sin.(Param.latN).^2,1,nz) .-
 #          Param.DeltaTh_z*log.(sigma).*repmat(cos.(Param.latN).^2,1,nz)) .*
 #          sigma.^(Param.Rd/Param.Cpd)

 #      matlab_max(X) = max.(X, Param.T_min)

 #      @views DeltaRhoT=(Param.k_a .+ (Param.k_s - Param.k_a) .* height_factor .*
 #        repmat(cos.(Param.latN).^4,1,nz)) .* U[:,:,Param.RhoPos] .*
 #        (T .- matlab_max(temp));
 #      @views F[:,:,Param.uPos:Param.vPos] .-= Param.k_f*height_factor.*U[:,:,Param.uPos:Param.vPos];
 #      if strcmp(Param.Thermo,"Energy")
 #      else
 #          @views  F[:,:,Param.ThPos] .-= DeltaRhoT.* (Param.p0./Pres).^Param.kappa;
 #      end
    elseif str == "heldsuarezcart"
       fac = Param.Rd / Param.p0
       sigma =  Pressure(Th,Rho,Th,Param) ./ Param.p0
       height_factor = max.(0, (sigma .- Param.sigma_b) ./ (1 .- Param.sigma_b))
       @views  F[:,:,ThPos]  .-=
            (Param.k_a .+ (Param.k_s .- Param.k_a) .* height_factor) .*
            Rho .*
            (sigma ./ (Rho .* fac) .- max.(
                Param.T_min,
                (Param.T_equator .- Param.DeltaT_y .- Param.DeltaTh_z .* log.(sigma)) .*
                sigma.^(Param.Rd / Param.Cpd),
                )
            ) ./ sigma.^Param.kappa

      @views  F[:,:,uPos:vPos] .-= (Param.k_f .* height_factor) .* U[:,:,uPos:vPos]
    end
end

