using NLsolve
using CLIMA.MoistThermodynamics
using CLIMA.PlanetParameters


function ResMoisture(z::FT,y::Array{Float64,1},yPrime::Array{Float64,1}) where {FT<:Real}

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
p_d=R_d*rho_d*T
Δcp=cp_v-cp_l
p_vs=saturation_vapor_pressure(T,FT(LH_v0),Δcp)
L=latent_heat_fusion(T)
F[1]=pPrime+grav*rho
F[2]=p-(R_d*rho_d+R_v*rho_v)*T
F[3]=theta_e-T*(p_d/MSLP)^(-R_d/(cp_d+cp_l*r_t))*exp(L*r_v/((cp_d+cp_l*r_t)*T))
F[4]=r_t-r_tP(z)
F[5]=rho_v-rho_d*r_v
F[6]=theta_e-theta_eP(z)
F[7]=rho_v-p_vs/(R_v*T)
F
end

function r_tP(z::FT)  where {FT<:Real}
 return 2.e-2
end

function theta_eP(z::FT)  where {FT<:Real}
320.e0
end

function SetImplEuler(z::FT,dz::FT,y0::Array{Float64,1}) where {FT<:Real}
  function ImplEuler(y::Array{Float64,1}) where {FT<:Real}
    ResMoisture(z,y,(y-y0)/dz)
  end
end

function TestRes()
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

F=SetImplEuler(z,dz,y0)
res=nlsolve(F,y0)
res.zero

dz=10.0
for i=1:1000
  y0=deepcopy(res.zero)
  F=SetImplEuler(z,dz,y0)
  res=nlsolve(F,y0)
  z=z+dz
end
y0


end

