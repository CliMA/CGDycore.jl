function PerturbMoistProfile(x, rho, rho_θ, rho_qv, rho_qc, Phys, Param)

  xc = Param.xc
  zc = Param.zc
  rc = Param.rc
  Δθ = Param.DeltaTh

  r = sqrt((x[1] - xc)^2 + (x[2] - zc)^2)
  rho_d = rho - rho_qv - rho_qc
  kappa_M = (Phys.Rd * rho_d + Phys.Rv * rho_qv) / (Phys.Cpd * rho_d + Phys.Cpv * rho_qv + Phys.Cpl * rho_qc)
  p_loc = Phys.p0 *(Phys.Rd * rho_θ / Phys.p0)^(1/(1-kappa_M))
  T_loc = p_loc / (Phys.Rd * rho_d + Phys.Rv * rho_qv)
  rho_e = (Phys.c_vd * rho_d + Phys.Cvv * rho_qv + Phys.Cpl * rho_qc) * T_loc + Phys.L00 * rho_qv

  if r < rc && Δθ > 0 
    θ_dens = rho_θ / rho * (p_loc / Phys.p0)^(kappa_M - Phys.kappa)
    θ_dens_new = θ_dens * (1 + Δθ * cospi(0.5*r/rc)^2 / 300)
    rt =(rho_qv + rho_qc) / rho_d 
    rv = rho_qv / rho_d
    θ_loc = θ_dens_new * (1 + rt)/(1 + (Phys.Rv / Phys.Rd) * rv)
    if rt > 0 
      while true 
        T_loc = θ_loc * (p_loc / Phys.p0)^Phys.kappa
        T_C = T_loc - 273.15
        # SaturVapor
        pvs = 611.2 * exp(17.62 * T_C / (243.12 + T_C))
        rho_d_new = (p_loc - pvs) / (Phys.Rd * T_loc)
        rvs = pvs / (Phys.Rv * rho_d_new * T_loc)
        θ_new = θ_dens_new * (1 + rt) / (1 + (Phys.Rv / Phys.Rd) * rvs)
        if abs(θ_new-θ_loc) <= θ_loc * 1.0e-12
          break
        else
          θ_loc=θ_new
        end
      end
    else
      rvs = 0
      T_loc = θ_loc * (p_loc / Phys.p0)^Phys.kappa
      rho_d_new = p_loc / (Phys.Rd * T_loc)
      θ_new = θ_dens_new * (1 + rt) / (1 + (Phys.Rv / Phys.Rd) * rvs)
    end
    rho_qv = rvs * rho_d_new
    rho_qc = (rt - rvs) * rho_d_new
    rho = rho_d_new * (1 + rt)
    rho_d = rho - rho_qv - rho_qc
    kappa_M = (Phys.Rd * rho_d + Phys.Rv * rho_qv) / (Phys.Cpd * rho_d + Phys.Cpv * rho_qv + Phys.Cpl * rho_qc)
    rho_θ = rho * θ_dens_new * (p_loc / Phys.p0)^(Phys.kappa - kappa_M)
    rho_e = (Phys.c_vd * rho_d + Phys.Cvv * rho_qv + Phys.Cpl * rho_qc) * T_loc + Phys.L00 * rho_qv

  end
  return rho, rho_e, rho_qv, rho_qc
end
