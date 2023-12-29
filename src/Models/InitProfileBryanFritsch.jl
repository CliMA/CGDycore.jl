function  LLv(T,Phys)

  Phys.L00 -(Phys.Cpl - Phys.Cpv) *(T - Phys.T0)
end

function  SatVap(T,Phys)
  Phys.p0 * (T / Phys.T0)^((Phys.Cpv - Phys.Cpl) / Phys.Rv) *
    exp((Phys.L00 / Phys.Rv) *(1.0 / Phys.T0 - 1.0 / T))
# T_C = T - 273.15
# 611.2 * exp(17.62 * T_C / (243.12 + T_C))  
end

function ResMoisture(z::FT,y::Array{FT,1},yPrime::Array{FT,1},Phys) where {FT<:Real}

  p=y[1]
  rho=y[2]
  T=y[3]
  r_t=y[4]
  r_v=y[5]
  rho_v=y[6]
  theta_e=y[7]
  pPrime=yPrime[1]
  F=zeros(7,1)

# rho_d=rho-rho_v
  rho_d=rho / (1 + r_t)
  p_d=Phys.Rd*rho_d*T
  Δcp=Phys.Cpv-Phys.Cpl
  p_vs=Thermodynamics.fpvs(T,Phys.T0)
  L=LLv(T,Phys)
  F[1]=pPrime+Phys.Grav*rho
  F[2]=p-(Phys.Rd*rho_d+Phys.Rv*rho_v)*T
  F[3]=theta_e-T*(p_d/Phys.p0)^(-Phys.Rd/(Phys.Cpd+Phys.Cpl*r_t))*exp(L*r_v/((Phys.Cpd+Phys.Cpl*r_t)*T))
  F[4]=r_t-r_tP(z)
  F[5]=rho_v-rho_d*r_v
  F[6]=theta_e-theta_eP(z)
  a = p_vs / (Phys.Rv * T) - rho_v
  b = rho - rho_v - rho_d
  F[7]=a+b-sqrt(a*a+b*b)
# @show F[7],T,a,p_vs,rho_v
# F[7]=rho_v-p_vs/(Phys.Rv*T)
  return F
end

function r_tP(z::FT)  where {FT<:Real}
 return 2.e-2
end

function theta_eP(z::FT)  where {FT<:Real}
320.e0
end

function SetImplEuler(z::FT,dz::FT,y0::Array{FT,1},Phys) where {FT<:Real}
  function ImplEuler(y::Array{FT,1}) where {FT<:Real}
    ResMoisture(z,y,(y-y0)/dz,Phys)
  end
end

function TestRes(Phys)
  y=zeros(7)
  p=1.e5
  rho=1.4
  r_t=r_tP(0)
  r_v=r_t
  rho_v=rho*r_v
  theta_e=theta_eP(0)
  T=theta_e

  yPrime=zeros(7)
  y0=zeros(7)
  y0[1]=p
  y0[2]=rho
  y0[3]=T
  y0[4]=r_t
  y0[5]=r_v
  y0[6]=rho_v
  y0[7]=theta_e

  z=0.0
  dz=0.01
  y=deepcopy(y0);

  F=SetImplEuler(z,dz,y0,Phys)
  res=nlsolve(F,y0)
  res.zero

  Prof = zeros(1001,5)
  p = y0[1]
  Rho = y0[2]
  T = y0[3]
  r_t = y0[4]
  r_v = y0[5]
  rho_v = y0[6]
  theta_e =y0[7]
  r_c = r_t - r_v
# Rho = RhoD + RhoD * r_t
  RhoD = Rho / (1.0 + r_t)
  RhoV = RhoD * r_v
  RhoC = RhoD * r_c
  ThetaV = Thermodynamics.fThetaV(Rho,RhoV,RhoC,T,Phys)
  Prof[1,1] = z
  Prof[1,2] = Rho
  Prof[1,3] = ThetaV
  Prof[1,4] = RhoV
  Prof[1,5] = RhoC
  dz=10.0
  for i=1:1000
    y0=deepcopy(res.zero)
    F=SetImplEuler(z,dz,y0,Phys)
    res=nlsolve(F,y0)
    z=z+dz
    p = y0[1]
    Rho = y0[2]
    T = y0[3]
    r_t = y0[4]
    r_v = y0[5]
    rho_v = y0[6]
    theta_e =y0[7]
    r_c = r_t - r_v
#   Rho = RhoD + RhoD * r_t
    RhoD = Rho / (1.0 + r_t)
    RhoV = RhoD * r_v
    RhoC = RhoD * r_c
    ThetaV = Thermodynamics.fThetaV(Rho,RhoV,RhoC,T,Phys)
    Prof[i+1,1] = z
    Prof[i+1,2] = Rho
    Prof[i+1,3] = ThetaV
    Prof[i+1,4] = RhoV
    Prof[i+1,5] = RhoC
#   FB=FischerBurmeisterRhoTh(Rho,RhoV,RhoC,Rho*ThetaV,Phys)
#   @show i,FB
#   stop
  end
  return Prof
end

function PerturbMoistProfile(x, Rho, RhoTheta, RhoV, RhoC, Phys, Param)

  xc = Param.xC0
  zc = Param.zC0
  rc = Param.rC0
  DeltaTh = Param.DeltaTh

  r = sqrt((x[1] - xc)^2 + (x[3] - zc)^2)
  Rho_d = Rho - RhoV - RhoC
  kappa_M = (Phys.Rd * Rho_d + Phys.Rv * RhoV) / (Phys.Cpd * Rho_d + Phys.Cpv * RhoV + Phys.Cpl * RhoC)
  p_loc = Phys.p0 *(Phys.Rd * RhoTheta / Phys.p0)^(1/(1-kappa_M))
  T_loc = p_loc / (Phys.Rd * Rho_d + Phys.Rv * RhoV)
  Theta = RhoTheta / Rho
  qv = RhoV / Rho
  qc = RhoC / Rho

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
    Theta = Theta_dens_new * (p_loc / Phys.p0)^(Phys.kappa - kappa_M)
    qv = RhoV / Rho
    qc = RhoC / Rho
    

  end
  return Rho, Theta, qv, qc
end

