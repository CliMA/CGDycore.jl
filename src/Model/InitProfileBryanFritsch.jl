function  LLv(T,Phys)

  Phys.L00 -(Phys.Cpl - Phys.Cpv) *(T - Phys.T0)
end

function  SatVap(T,Phys)
  Phys.p0 * (T / Phys.T0)^((Phys.Cpv - Phys.Cpl) / Phys.Rv) *
    exp((Phys.L00 / Phys.Rv) *(1.0 / Phys.T0 - 1.0 / T))
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

  rho_d=rho-rho_v
  p_d=Phys.Rd*rho_d*T
  Δcp=Phys.Cpv-Phys.Cpl
  p_vs=SatVap(T,Phys)
  L=LLv(T,Phys)
  F[1]=pPrime+Phys.Grav*rho
  F[2]=p-(Phys.Rd*rho_d+Phys.Rv*rho_v)*T
  F[3]=theta_e-T*(p_d/Phys.p0)^(-Phys.Rd/(Phys.Cpd+Phys.Cpl*r_t))*exp(L*r_v/((Phys.Cpd+Phys.Cpl*r_t)*T))
  F[4]=r_t-r_tP(z)
  F[5]=rho_v-rho_d*r_v
  F[6]=theta_e-theta_eP(z)
  F[7]=rho_v-p_vs/(Phys.Rv*T)
  return F
end

function r_tP(z::FT)  where {FT<:Real}
 return 2.e-2
end

function theta_eP(z::FT)  where {FT<:Real}
320.e0
end

function SetImplEuler(z::FT,dz::FT,y0::Array{Float64,1},Phys) where {FT<:Real}
  function ImplEuler(y::Array{Float64,1}) where {FT<:Real}
    ResMoisture(z,y,(y-y0)/dz,Phys)
  end
end

function TestRes(Phys)
y=zeros(7)
p=1.e5
rho=1.4
r_t=r_tP(0)
println("r_t  ",r_t)
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

Prof = zeros(1001,8)
@. Prof[1,2:end] = y0
Prof[1,1] = z
dz=10.0
for i=1:1000
  y0=deepcopy(res.zero)
  F=SetImplEuler(z,dz,y0,Phys)
  res=nlsolve(F,y0)
  z=z+dz
  @. Prof[i+1,2:end] = y0
  Prof[i+1,1] = z
end
  
  return Prof

end

