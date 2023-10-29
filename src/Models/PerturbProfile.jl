function PerturbMoistProfile(x, Rho, RhoTheta, RhoV, RhoC, Phys, Param)

  xc = Param.xc
  zc = Param.zc
  rc = Param.rc
  DeltaTh = Param.DeltaTh

  r = sqrt((x[1] - xc)^2 + (x[3] - zc)^2)
  Rho_d = Rho - RhoV - RhoC
  kappa_M = (Phys.Rd * Rho_d + Phys.Rv * RhoV) / (Phys.Cpd * Rho_d + Phys.Cpv * RhoV + Phys.Cpl * RhoC)
  p_loc = Phys.p0 *(Phys.Rd * RhoTheta / Phys.p0)^(1/(1-kappa_M))
  T_loc = p_loc / (Phys.Rd * Rho_d + Phys.Rv * RhoV)
  Rho_e = (Phys.Cvd * Rho_d + Phys.Cvv * RhoV + Phys.Cpl * RhoC) * T_loc + Phys.L00 * RhoV

  if r < rc && DeltaTh > 0 
    θ_dens = RhoTheta / Rho * (p_loc / Phys.p0)^(kappa_M - Phys.kappa)
    Theta_dens_new = θ_dens * (1 + DeltaTh * cospi(0.5*r/rc)^2 / 300)
    rt =(RhoV + RhoC) / Rho_d 
    rv = RhoV / Rho_d
    Theta_loc = Theta_dens_new * (1 + rt)/(1 + (Phys.Rv / Phys.Rd) * rv)
    if rt > 0 
      while true 
        T_loc = Theta_loc * (p_loc / Phys.p0)^Phys.kappa
        # SaturVapor
        pvs = fpvs(T_loc,Phys)
        Rho_d_new = (p_loc - pvs) / (Phys.Rd * T_loc)
        rvs = pvs / (Phys.Rv * Rho_d_new * T_loc)
        Theta_new = Theta_dens_new * (1 + rt) / (1 + (Phys.Rv / Phys.Rd) * rvs)
        if abs(Theta_new-Theta_loc) <= Theta_loc * 1.0e-12
          break
        else
          Theta_loc=Theta_new
        end
      end
    else
      rvs = 0
      T_loc = Theta_loc * (p_loc / Phys.p0)^Phys.kappa
      Rho_d_new = p_loc / (Phys.Rd * T_loc)
      Theta_new = Theta_dens_new * (1 + rt) / (1 + (Phys.Rv / Phys.Rd) * rvs)
    end
    RhoV = rvs * Rho_d_new
    RhoC = (rt - rvs) * Rho_d_new
    Rho = Rho_d_new * (1 + rt)
    Rho_d = Rho - RhoV - RhoC
    kappa_M = (Phys.Rd * Rho_d + Phys.Rv * RhoV) / (Phys.Cpd * Rho_d + Phys.Cpv * RhoV + Phys.Cpl * RhoC)
    RhoTheta = Theta_dens_new * (p_loc / Phys.p0)^(Phys.kappa - kappa_M)

  end
  return Rho, Theta, RhoV, RhoC
end
