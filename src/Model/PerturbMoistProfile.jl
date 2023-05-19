#=
function PerturbMoistProfile(x, Rho, RhoTh, rho_qv, rho_qc, Phys, Param)

  xc = Param.xc
  zc = Param.zc
  rc = Param.rc
  r = sqrt((x[1] - Param.xc)^2 + (x[3] - Param.zc)^2)
  Rho_d = Rho - Rho_qv - Rho_qc
  kappa_M = (Phys.Rd * Rho_d + Phys.Rv * Rho_qv) / (Phys.Cpd * Rho_d + Phys.Cpv * Rho_qv + Phys.Cpl * Rho_qc)
  p_loc = Phys.p0 *(Phys.Rd * RhoTh / Phys.p0)^(1/(1-kappa_M))
  T_loc = p_loc / (Phys.Rd * Rho_d + Phys.Rv * Rho_qv)
  Rho_e = (Phys.Cvd * Rho_d + Phys.Cvv * Rho_qv + Phys.Cpl * Rho_qc) * T_loc + Phys.L00 * Rho_qv

  if r < rc && DeltaTh > 0 
    ThDens = RhoTh / Rho * (p_loc / Phys.p0)^(kappa_M - Phys.kappa)
    ThDens_new = ThDens * (1 + DeltaTh * cospi(0.5*r/rc)^2 / 300)
    rt =(Rho_qv + Rho_qc) / Rho_d 
    rv = Rho_qv / Rho_d
    ThLoc = ThDens_new * (1 + rt)/(1 + (Phys.Rv / Phys.Rd) * rv)
    if rt > 0 
      while true 
        T_loc = ThLoc * (p_loc / Phys.p0)^Phys.kappa
        T_C = T_loc - 273.15
        # SaturVapor
        pvs = 611.2 * exp(17.62 * T_C / (243.12 + T_C))
        Rho_d_new = (p_loc - pvs) / (Phys.Rd * T_loc)
        rvs = pvs / (Phys.Rv * Rho_d_new * T_loc)
        ThNew = ThDens_new * (1 + rt) / (1 + (Phys.Rv / Phys.Rd) * rvs)
        if abs(ThNew-ThLoc) <= ThLoc * 1.0e-12
          break
        else
          ThLoc=ThNew
        end
      end
    else
      rvs = 0
      T_loc = ThLoc * (p_loc / Phys.p0)^Phys.kappa
      Rho_d_new = p_loc / (Phys.Rd * T_loc)
      ThNew = ThDens_new * (1 + rt) / (1 + (Phys.Rv / Phys.Rd) * rvs)
    end
    Rho_qv = rvs * Rho_d_new
    Rho_qc = (rt - rvs) * Rho_d_new
    Rho = Rho_d_new * (1 + rt)
    Rho_d = Rho - Rho_qv - Rho_qc
    kappa_M = (Phys.Rd * Rho_d + Phys.Rv * Rho_qv) / (Phys.Cpd * Rho_d + Phys.Cpv * Rho_qv + Phys.Cpl * Rho_qc)
    RhoTh = Rho * ThDens_new * (p_loc / Phys.p0)^(Phys.kappa - kappa_M)
    Rho_e = (Phys.Cvd * Rho_d + Phys.Cvv * Rho_qv + Phys.Cpl * Rho_qc) * T_loc + Phys.L00 * Rho_qv
  end  
end    
=#
