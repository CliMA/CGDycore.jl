using NLsolve 
#using CLIMA.MoistThermodynamics
#using CLIMA.PlanetParameters

cp_d = 1004.0
cv_d = 717.0
cp_v = 1885.0 
cv_v = 1424.0
cp_l = 4186.0
R_d = cp_d - cv_d
R_v = cp_v - cv_v
MSLP = 1e5
grav = 9.81
z_i = 795.0

function saturation_vapor_pressure(T)
T_C=T-273.15
6.112e2*exp(17.62*T_C/(243.12+T_C))
end

function ResMoisture(z::FT,y::Array{Float64,1},yPrime::Array{Float64,1}) where {FT<:Real}

p=y[1]
ρ=y[2]
T=y[3]
q_t=y[4]
θ=y[5]
θ_l=y[6]
pPrime=yPrime[1]
F=zeros(6,1)


ρ_v=q_t*ρ
ρ_d=ρ-ρ_v
p_d=R_d*ρ_d*T
Δcp=cp_v-cp_l
p_vs=saturation_vapor_pressure(T)
q_v=q_t
q_l=0
F[1]=pPrime+grav*ρ
F[2]=p-(R_d*ρ_d+R_v*ρ_v)*T
F[3]=θ-T*(p/1.e5)^(-R_d/cp_d)
F[4]=θ_l-θ*exp(-FT(2.47e6)/(cp_d*T)*q_l/(1-q_t))
F[5]=q_t-q_tP(z)
F[6]=θ_l-θ_lP(z)
return F
end

function q_tP(z::FT)  where {FT<:Real}
  if z<=z_i
    q_t=9.45e-3
  else
    q_t=(5.0-3.0*(1-exp(z/z_i))/500.0)*1.e-3  
  end
end

function θ_lP(z::FT)  where {FT<:Real}
  if z<=z_i
    θ_l=288.3
  else
    θ_l=297.5+(z-z_i)^(1.0/3.0)
  end
end

function uP(z::FT)  where {FT<:Real}
  uP=3.0+4.3*z/1000.0
end

function vP(z::FT)  where {FT<:Real}
  vP=-9.0+5.6*z/1000.0
end

function SetImplEuler(z::FT,dz::FT,y0::Array{Float64,1}) where {FT<:Real}
  function ImplEuler(y::Array{Float64,1}) where {FT<:Real}
    ResMoisture(z+dz,y,(y-y0)/dz)
  end
end

function TestRes()
y=zeros(7)
p=MSLP
ρ=1.4
q_t=q_tP(0)
q_v=q_t
θ_l=θ_lP(0)
T=θ_l
θ=θ_l

yPrime=zeros(7)
y0=zeros(7)
y0[1]=p
y0[2]=ρ
y0[3]=T
y0[4]=q_t
y0[5]=θ
y0[6]=θ_l

z=0.0
dz=0.01
F=SetImplEuler(z,dz,y0)
res=nlsolve(F,y0)
#res.zero

println("    p         z         T         T         rH     qV     WS    WD    qt    u     v")
p=y0[1]
ρ=y0[2]
T=y0[3]
q_t=y0[4]
rH=R_v*ρ*q_v*T/saturation_vapor_pressure(T)
println(p," ",z," ",T," ",T," ",rH," ",sqrt(uP(z)^2+vP(z)^2)," ",asind(uP(z)/sqrt(uP(z)^2+vP(z)^2))," ",q_t," ",
  uP(z)," ",vP(z))
dz=z_i/84
for i=1:84
  y0=deepcopy(res.zero)
  F=SetImplEuler(z,dz,y0)
  res=nlsolve(F,y0)
  z=z+dz
# println("z ",z," p ",y0[1], " q_l ",y0[4]-y0[5]," q_v ",y0[5]," T ",y0[3])
# println(z," ",y0[1]," ",max(y0[4]-y0[5],0.0)," ",y0[5]," ",y0[3]," ",sqrt(uP(z)^2+vP(z)^2)," ",asind(uP(z)/sqrt(uP(z)^2+vP(z)^2)))
  p=y0[1]
  ρ=y0[2]
  T=y0[3]
  q_t=y0[4]
  rH=R_v*ρ*q_v*T/saturation_vapor_pressure(T)
  println(p," ",z," ",T," ",T," ",rH," ",sqrt(uP(z)^2+vP(z)^2)," ",asind(uP(z)/sqrt(uP(z)^2+vP(z)^2))," ",q_t," ",
  uP(z)," ",vP(z))
end
z=z_i
dz=0.01
y0=deepcopy(res.zero)
y0[3]=T
y0[4]=q_tP(z+dz)
y0[5]=θ_lP(z+dz)
y0[6]=θ_lP(z+dz)
F=SetImplEuler(z,dz,y0)
res=nlsolve(F,y0)
#println("z ",z," p ",y0[1], " q_l ",y0[4]-y0[5]," q_v ",y0[5])
dz=10.0
for i=1:71
  y0=deepcopy(res.zero)
  F=SetImplEuler(z,dz,y0)
  res=nlsolve(F,y0)
  z=z+dz
# println("z ",z," p ",y0[1], " q_l ",y0[4]-y0[5]," q_v ",y0[5])
  p=y0[1]
  ρ=y0[2]
  T=y0[3]
  q_t=y0[4]
  rH=R_v*ρ*q_v*T/saturation_vapor_pressure(T)
   println(p," ",z," ",T," ",T," ",rH," ",sqrt(uP(z)^2+vP(z)^2)," ",asind(uP(z)/sqrt(uP(z)^2+vP(z)^2))," ",q_t," ",
   uP(z)," ",vP(z))
end


end

TestRes()
