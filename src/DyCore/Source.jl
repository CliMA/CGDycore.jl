function Source!(F,U,CG,Global,iG)
  Model = Global.Model
  Param = Global.Model.Param
  Phys = Global.Phys
  RhoPos = Model.RhoPos
  uPos = Model.uPos
  vPos = Model.vPos
  ThPos = Model.ThPos
  nz = Global.Grid.nz

  str = lowercase(Model.ProfRho)
  @views Rho = U[:,Model.RhoPos]
  @views Th = U[:,Model.ThPos]
  @views Tr = U[:,Model.NumV+1:end]
  @views Sigma = Global.Cache.Cache1[:,1]
  @views height_factor = Global.Cache.Cache2[:,1]
  @views ΔρT = Global.Cache.Cache3[:,1]
  if str == "heldsuarezsphere"
    Pressure!(Sigma,Th,Rho,Tr,Global)
    Sigma = Sigma / Phys.p0
    @. height_factor = max(0.0, (Sigma - Param.sigma_b) / (1.0 - Param.sigma_b))
    @. ΔρT =
      (Param.k_a + (Param.k_s - Param.k_a) * height_factor *
      cos(Global.latN[iG])^4) *
      Rho *
      (                  # ᶜT - ᶜT_equil
         Phys.p0 * Sigma / (Rho * Phys.Rd) - 
         max(Param.T_min,
         (Param.T_equator - Param.DeltaT_y * sin(Global.latN[iG])^2 - 
         Param.DeltaTh_z * log(Sigma) * cos(Global.latN[iG])^2) *
         Sigma^Phys.kappa)
      )
     @views @. F[:,uPos:vPos] -= (Param.k_f * height_factor) * U[:,uPos:vPos]
     @views @. F[:,ThPos]  -= ΔρT / Sigma^Phys.kappa
  end
end

function SourceMicroPhysics(F,U,CG,Global,iG)
  (; Rd,
     Cpd,
     Rv,
     Cpv,
     Cpl,
     p0,
     L00,
     kappa) = Global.Phys
  ThPos=Global.Model.ThPos
  RhoPos=Global.Model.RhoPos
  RhoVPos=Global.Model.RhoVPos
  RhoCPos=Global.Model.RhoCPos
  RelCloud = Global.Model.Param.RelCloud
  NumV=Global.Model.NumV
  nz = Global.Grid.nz
  for i = 1:nz
    (FTh,FRho,FRhoV,FRhoC) = 
      Microphysics( U[i,ThPos],U[i,RhoPos], U[i,NumV+RhoVPos], U[i,NumV+RhoCPos], 
        Rd,Cpd,Rv,Cpv,Cpl,L00,p0,RelCloud)
    F[i,ThPos] += FTh   
    F[i,RhoPos] += FRho
    F[i,RhoVPos+NumV] += FRhoV
    F[i,RhoCPos+NumV] += FRhoC
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
  FRhoV = RelCloud * (a + b - sqrt(a * a + b * b))
  L = L00 - (Cpl - Cpv) * T
  FRhoTh = RhoTh * ((-L / (Cpml * T) -
                        log(p / p0) * kappaM * (Rv / Rm - Cpv / Cpml) +
                        Rv / Rm) * FRhoV +
                        (log(p / p0) * kappaM * (Cpl / Cpml)) * (-FRhoV))
  FRho = FRhoV
  FRhoC = 0.0
  return (FRhoTh,FRho,FRhoV,FRhoC)
end  


