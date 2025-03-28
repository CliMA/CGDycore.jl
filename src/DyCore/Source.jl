using NLsolve
function Source!(F,U,Pres,CG,Metric,Cache,Phys,Global,Param,iG)
  Model = Global.Model
  RhoPos = Model.RhoPos
  uPos = Model.uPos
  vPos = Model.vPos
  NumV = Model.NumV
  NumTr = Model.NumTr
  ThPos = Model.ThPos
  nz = Global.Grid.nz

  @views Rho = U[:,RhoPos]
  @views Th = U[:,ThPos]
  @views Tr = U[:,NumV+1:NumV+NumTr]
# @time if Global.Model.Problem  == "HeldSuarezSphere" || 
#   Global.Model.Problem == "HeldSuarezMoistSphere" || 
#   Global.Model.Problem  == "HeldSuarezSphereOro" || 
#   Global.Model.Problem == "HeldSuarezMoistSphereOro"
    @views SourceHeldSuarez!(F[:,ThPos],F[:,uPos:vPos],Rho,Th,U[:,uPos:vPos],Tr,Pres,
      Param.sigma_b,Param.k_s,Param.k_a,Param.k_f,Param.T_min,Param.T_equator,
      Param.DeltaT_y,Param.DeltaTh_z,
      Metric.lat[iG],Cache,Phys)
# end
end

function SourceHeldSuarez!(FTh,FV,Rho,Th,V,Tr,Pres,
  sigma_b,k_s,k_a,k_f,T_min,T_equator,DeltaT_y,DeltaTh_z,latN,Cache,Phys)
  @views Sigma = Cache.Cache1[:,1]
  @views height_factor = Cache.Cache2[:,1]
  @views Teq = Cache.Cache2[:,1]
  @views ΔρT = Cache.Cache3[:,1]
  @views kT = Cache.Cache3[:,1]
  @. Sigma = Pres
  @. Sigma = Sigma / Phys.p0
  @. height_factor = max(0.0, (Sigma - sigma_b) / (1.0 - sigma_b))
  @. FV -= (k_f * height_factor) * V 
  @. kT =
    (Sigma < 0.7) * (k_a + (k_s - k_a) * height_factor * cos(latN)^4) # Oswald 
  @. Teq = (T_equator - DeltaT_y * sin(latN)^2 - DeltaTh_z * log(Sigma) * cos(latN)^2) * Sigma^Phys.kappa
  @. Teq = max(T_min, Teq)
  @. ΔρT =  kT * (Phys.p0 * Sigma / (Rho * Phys.Rd) - Teq)
  @. FTh  -= Rho * ΔρT / Sigma^Phys.kappa

end  


function SourceHeldSuarez1!(FTh,FV,Rho,Th,V,Tr,latN,Global)
  Phys = Global.Phys
  @views Sigma = Global.Cache.Cache1[:,1]
  @views height_factor = Global.Cache.Cache2[:,1]
  @views ΔρT = Global.Cache.Cache3[:,1]
  day = 3600.0 * 24.0
  k_a=1.0/(40.0 * day)
  k_f=1.0/day
  k_s=1.0/(4.0 * day)
  DeltaT_y=0.0
  DeltaTh_z=-5.0
  T_equator=315.0
  T_min=200.0
  sigma_b=7.0/10.0
  Pressure!(Sigma,Th,Rho,Tr,Global)
  @. Sigma = Sigma / Phys.p0
  @. height_factor = max(0.0, (Sigma - sigma_b) / (1.0 - sigma_b))
  @. FV -= (k_f * height_factor) * V
  @. height_factor = (T_equator - DeltaT_y * sin(latN)^2 -   
      DeltaTh_z * log(Sigma) * cos(latN)^2) *
      Sigma^Phys.kappa
  @. ΔρT =
    (k_a + (k_s - k_a) * height_factor *
    cos(latN)^4) *
    Rho 
  @. ΔρT *=  Phys.p0 * Sigma / (Rho * Phys.Rd) - 
     max(T_min,height_factor)
  @. FTh  -= ΔρT / Sigma^Phys.kappa
end  


function SourceMicroPhysics(F,U,Pres,CG,Phys,Global,iG)
  (; Rd,
     Cpd,
     Rv,
     Cpv,
     Cpl,
     p0,
     kappa) = Phys
   ThPos=Global.Model.ThPos
   RhoPos=Global.Model.RhoPos
   RhoVPos=Global.Model.RhoVPos
   RhoCPos=Global.Model.RhoCPos
   RelCloud = Global.Model.RelCloud
   NumV=Global.Model.NumV
   nz = Global.Grid.nz
   @inbounds for i = 1:nz
     Rho = U[i,RhoPos]  
     RhoTh = U[i,ThPos]  
     RhoV = max(U[i,RhoVPos], 0.0)  
     RhoC = max(U[i,RhoCPos], 0.0)  
     RhoD = Rho - RhoV - RhoC
     Cpml = Cpd * RhoD + Cpv * RhoV + Cpl * RhoC
     Rm = Rd * RhoD + Rv * RhoV
     kappaM = Rm / Cpml
     p = Pres[i]
     T = p / Rm
     p_vs = Models.fpws(T,Phys.T0)
     a = p_vs / (Rv * T) - RhoV
     b = RhoC
     FPh = 0.5 * RelCloud * (a + b - sqrt(a * a + b * b))
     L = Models.LatHeat(T,Phys)
     FR = -FPh * Global.Model.Rain
     if Global.Model.Thermo == "InternalEnergy"
     else    
       FRhoTh = RhoTh*((-L/(Cpml*T) - log(p / p0) * (Rm / Cpml) * (Rv / Rm + (Cpl -Cpv) / Cpml)  + Rv / Rm) * FPh +
         (1.0 / Rho - log(p/p0) * (Rm / Cpml) * (Cpl / Cpml)) * FR)
       F[i,ThPos] += FRhoTh   
     end  
     F[i,RhoPos] += -FR
     F[i,RhoVPos+NumV] += FPh
     F[i,RhoCPos+NumV] += -FPh - FR
  end  
end

function Microphysics(RhoTh,Rho,RhoV,RhoC,Rd,
     Cpd,Rv,Cpv,Cpl,L00,p0,RelCloud)
  RhoD = Rho - RhoV - RhoC
  Cpml = Cpd * RhoD + Cpv * RhoV + Cpl * RhoC
  Rm = Rd * RhoD + Rv * RhoV
  kappaM = Rm / Cpml
  p = (Rd * RhoTh / p0^kappaM)^(1 / (1 - kappaM))
  T = p / Rm
  T_C = T - 273.15
  p_vs = 611.2 * exp(17.62 * T_C / (243.12 + T_C))
  a = p_vs / (Rv * T) - RhoV
  b = 0.0
  FRhoV = 0.5 * RelCloud * (a + b - sqrt(a * a + b * b))
  L = L00 - (Cpl - Cpv) * T
  FRhoTh = RhoTh*(-L/(Cpml*T) - log(p / p0) * (Rm / Cpml) *(Rv / Rm  - Cpv / Cpml) + Rv / Rm)*FRhoV
# PotM*( (-Lv/(Cpml*T+Eps) &
#                       -LOG((p+Eps)/P0)*(Rm/Cpml)*(Rv/Rm-Cpv/Cpml) &
#                       +Rv/Rm                                   &
#                       )*Qcond                                   &
#                       +(LOG((p+Eps)/P0)*(Rm/Cpml)*(Cpl/Cpml)        &
#                       )*(-Qcond)                                   &
#                     )
  FRho = FRhoV
  FRhoC = 0.0
  return (FRhoTh,FRho,FRhoV,FRhoC)
end  


function MicrophysicsDCMIP(RhoTh,Rho,RhoV,RhoC,Rd,
     Cpd,Rv,Cpv,Cpl,L00,p0,RelCloud)
  RhoD = Rho - RhoV - RhoC
  Cpml = Cpd * RhoD + Cpv * RhoV + Cpl * RhoC
  Rm = Rd * RhoD + Rv * RhoV
  kappaM = Rm / Cpml
  p = (Rd * RhoTh / p0^kappaM)^(1 / (1 - kappaM))
  T = p / Rm
  T_C = T - 273.15
  p_vs = 611.2 * exp(17.62 * T_C / (243.12 + T_C))
  RhoVS = p_vs / (Rv * T) 
  if RhoV > RhoVS
    L = L00 - (Cpl - Cpv) * T
    Cond = (RhoV - RhoVS) / (1.0+(L/Cpd)*(L*(RhoVS/Rho)/(Rv*T^2)))
    FRhoV = -RelCloud * Cond
    FRhoTh = RhoTh*(-L/(Cpml*T) - log(p / p0) * (Rm / Cpml) *(Rv / Rm - Cpv / Cpml) + Rv / Rm)*FRhoV
    FRho = FRhoV
  else
    FRhoV = 0.0
    FRhoTh = 0.0
    FRho = 0.0
  end  
  FRhoC = 0.0
  return (FRhoTh,FRho,FRhoV,FRhoC)
end

function MicrophysicsSat(RhoTh,Rho,RhoV,RhoC,Rd,
     Cpd,Rv,Cpv,Cpl,L00,p0,RelCloud)
  RhoD = Rho - RhoV - RhoC
  Cpml = Cpd * RhoD + Cpv * RhoV + Cpl * RhoC
  Rm = Rd * RhoD + Rv * RhoV
  kappaM = Rm / Cpml
  p = (Rd * RhoTh / p0^kappaM)^(1 / (1 - kappaM))
  T = p / Rm
  T_C = T - 273.15
  p_vs = 611.2 * exp(17.62 * T_C / (243.12 + T_C))
  RhoVS = p_vs / (Rv * T) 
  if RhoV > RhoVS
    F=SetResT(T,Rho,RhoV,Rd,Cpd,Rv,Cpv,Cpl,L00,p0)
    res=nlsolve(F,[0.0,0.0])
    Cond = -res.zero[2]
    L = L00 - (Cpl - Cpv) * T
    FRhoV = -RelCloud * Cond
    FRhoTh = RhoTh*(-L/(Cpml*T) - log(p / p0) * (Rm / Cpml) *(Rv / Rm - Cpv / Cpml) + Rv / Rm)*FRhoV
    FRho = FRhoV
  else
    FRhoV = 0.0
    FRhoTh = 0.0
    FRho = 0.0
  end  
  FRhoC = 0.0
  return (FRhoTh,FRho,FRhoV,FRhoC)
end

function SetRes(RhoTh0, Rho0, RhoV0, Rd, Cpd, Rv, Cpv, Cpl, L00, p0) 
  function ResMoisture(y)
    dRhoTh = y[1]
    dRhoV = y[2]
    RhoTh = RhoTh0 + dRhoTh
    RhoV = RhoV0 + dRhoV
    Rho = Rho0 + dRhoV
    RhoC = -dRhoV
    RhoD = Rho - RhoV - RhoC
    Cpml = Cpd * RhoD + Cpv * RhoV + Cpl * RhoC
    Rm = Rd * RhoD + Rv * RhoV
    kappaM = Rm / Cpml
    p = (Rd * RhoTh / p0^kappaM)^(1 / (1 - kappaM))
    T = p / Rm
    T_C = T - 273.15
    p_vs = 611.2 * exp(17.62 * T_C / (243.12 + T_C))
    RhoVS = p_vs / (Rv * T)
    L = L00 - (Cpl - Cpv) * T
    FdRhoV = RhoVS - RhoV
    FdRhoTh = dRhoTh - RhoTh*(-L/(Cpml*T) - log(p / p0) * (Rm / Cpml) *(Rv / Rm - Cpv / Cpml) + Rv / Rm) * dRhoV
    return [FdRhoTh,FdRhoV]
  end
end

function SetResT(T0, Rho0, RhoV0, Rd, Cpd, Rv, Cpv, Cpl, L00, p0) 
  function ResMoisture(y)
    dT = y[1]
    dRhoV = y[2]
    T = T0 + dT
    RhoV = RhoV0 + dRhoV
    T_C = T - 273.15
    p_vs = 611.2 * exp(17.62 * T_C / (243.12 + T_C))
    RhoVS = p_vs / (Rv * T)
    L = L00 - (Cpl - Cpv) * T
    FdRhoV = RhoVS - RhoV
    FdT = dT + L/Cpd * dRhoV
    return [FdT,FdRhoV]
  end
end

