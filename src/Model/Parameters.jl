Base.@kwdef struct ParamGalewskiSphere
  H0G = 10000.0
  hH = 120.0
  alphaG = 1.0/3.0
  betaG = 1.0/15.0
  lat0G = pi/7.0
  lat1G = pi/2.0-lat0G
  eN = exp(-4.0/(lat1G-lat0G)^2.0)
  uM = 80.0
  Omega = 2*pi/24.0/3600.0 
end

Base.@kwdef struct ParamBaroWaveDrySphere
  T0E=310.0
  T0P=240.0
  B=2.0
  K=3.0
  LapseRate=0.005
  U0=-0.5
  PertR=1.0/6.0
  Up=1.0
  PertExpR=0.1
  PertLon=pi/9.0
  PertLat=2.0*pi/9.0
  PertZ=15000.0
  NBr=1.e-2
  DeltaT=1
  ExpDist=5
  T0=300
  TEq=300
  T_init = 315
  lapse_rate = -0.008
  Deep=false
  pert = 0.1
  uMax = 1.0
  vMax = 0.0
  DeltaT_y=0
  DeltaTh_z=-5
  T_equator=315
  T_min=200
  sigma_b=7/10
  z_D=20.0e3
  #      Moist
  q_0 = 0.018                # Maximum specific humidity (default: 0.018)
  q_t = 1.0e-12
end

Base.@kwdef struct ParamHeldSuarezDrySphere
  day = 3600.0 * 24.0
  k_a=1.0/(40.0 * day)
  k_f=1.0/day
  k_s=1.0/(4.0 * day)
  DeltaT_y=60.0
  DeltaTh_z=10.0
  T_equator=315.0
  T_min=200.0
  sigma_b=7.0/10.0
  CE = 0.0044
  CH = 0.0044
  CTr = 0.004
  p_pbl = 85000.0
  p_strato = 10000.0
  T_virt_surf = 290.0
  T_min_ref = 220.0
  H_t = 8.e3
  q_0 = 0.018                # Maximum specific humidity (default: 0.018)
  q_t = 1e-12
  T0E = 310.0
  T0P = 240.0
  B=2.0
  K=3.0
  LapseRate=0.005
  DeltaTS = 29.0
  TSMin = 271.0
  DeltaLat = 26.0 * pi / 180.0
  uMax = 0.0
  vMax = 0.0
end

Base.@kwdef struct ParamHillSchaerCart
  Example = "HillSchaerCart"
  Deep=false
  NBr=1.e-2
  Th0=300.0
  uMax=10
  vMax=0
  TEq=300.0
  Stretch = false
end

Base.@kwdef struct ParamHillAgnesiXCart
  Example = "HillAgnesiXCart"
  Deep=false
  NBr=1.e-2
  Th0=300.0
  uMax=10
  vMax=0
  wMax=0
  TEq=300.0
  a = 1000.0
  h = 400.0
  xc = 0.0
  Stretch = false
end

Base.@kwdef struct ParamHillAgnesiYCart
  Example = "HillAgnesiYCart"
  Deep=false
  NBr=1.e-2
  Th0=300.0
  uMax=0
  vMax=10
  wMax=0
  TEq=300.0
  a = 1000.0
  h = 400.0
  yc = 0.0
  Stretch = false
end

Base.@kwdef struct ParamWarmBubble2DXCart
  Example = "WarmBubble2DXCart"
  Th0=300.0
  uMax=0.0
  vMax=0
  wMax=0
  DeltaTh = 2.0
  xC0 = 10000.0
  zC0 = 2000.0
  rC0 = 2000.0
  Stretch = false
end

Base.@kwdef struct ParamHillGaussCart
  Deep=false
  NBr=1.e-2
  Th0=300.0
  uMax=0
  vMax=0
  TEq=300.0
  Stretch = false
end

H = 1.2e4
Cpd=1004.0e0
Cvd=717.0e0
Rd=Cpd-Cvd
T_0 = 300.0
Grav = 9.81e0
ScaleHeight = Float64(Rd * T_0 / Grav)
tau = 1036800.0
omega_0 = 23000 * pi / tau
RadEarth = 6.37122e+6
p0 = 1.e5
p_top = p0 * exp(-H / ScaleHeight)
Base.@kwdef struct ParamAdvectionSphereGaussian
  TimeDependent = true
  hMax = 0.95
  b = 5.0
  lon1 = 5.0e0 / 6.0e0 * pi
  lat1 = 0.0e0
  lon2 = 7.0e0 / 6.0e0 * pi
  lat2 = 0.0e0
  EndTime = 12.0 * 24.0 * 3600.0
  FacVel = 10.0
  StreamFun = true
end  
Base.@kwdef struct ParamAdvectionSphereSlottedCylinder
  TimeDependent = true
  hMax = 0.95
  b = 5.0
  lon1 = 5.0e0 / 6.0e0 * pi
  lat1 = 0.0e0
  lon2 = 7.0e0 / 6.0e0 * pi
  lat2 = 0.0e0
  EndTime = 1.0368e6 
  FacVel = 10.0
  StreamFun = true
end  
Base.@kwdef struct ParamAdvectionSphereDCMIP
  xC = 0.0
  H = H
  R_t = RadEarth / 2.0
  Z_t = 1000.0
  z_c = 5.0e3
  p_top = p_top
  T_0 = T_0
  ScaleHeight = ScaleHeight
  b = 0.2
  Lon_c1 = 150.0/360.0 * 2 * pi
  Lon_c2 = 210.0/360.0 * 2 * pi
  Lat_c = 0.0
  tau = tau
  omega_0 = omega_0
  TimeDependent = true
end

Base.@kwdef struct ParamAdvectionCubeCart
  Example = "AdvectionCubeCart"
  StreamFun = false
  uMax = 0.25
  vMax = 0.25
  x1 = 400.0
  x2 = 600.0
  y1 = 400.0
  y2 = 600.0
end  

Base.@kwdef struct ParamAdvectionCart
  Example = "AdvectionCart"
  xC = 0.0
  H = H
end  



function Parameters(Problem::String)
  if Problem == "BaroWaveDrySphere" || Problem == "BaroWaveDrySphereOro"
    @show Problem
    Param = ParamBaroWaveDrySphere()
  elseif Problem == "GalewskiSphere"
    @show Problem
    Param = ParamGalewskiSphere()
  elseif Problem == "HeldSuarezDrySphere" || Problem == "HeldSuarezDrySphereOro"
    @show Problem
    Param = ParamHeldSuarezDrySphere()
  elseif Problem == "HillSchaerCart"
    @show Problem
    Param = ParamHillSchaerCart()
  elseif Problem == "HillAgnesiXCart"
    @show Problem
    Param = ParamHillAgnesiXCart()
  elseif Problem == "HillAgnesiYCart"
    @show Problem
    Param = ParamHillAgnesiYCart()
  elseif Problem == "HillGaussCart"
    @show Problem
    Param = ParamHillGaussCart()
  elseif Problem == "AdvectionSphereDCMIP"
    @show Problem
    Param = ParamAdvectionSphereDCMIP()
  elseif Problem == "AdvectionSphereGaussian"
    @show Problem
    Param = ParamAdvectionSphereGaussian()
  elseif Problem == "AdvectionSphereSlottedCylinder"
    @show Problem
    Param = ParamAdvectionSphereSlottedCylinder()
  elseif Problem == "AdvectionCart"
    @show Problem
    Param = ParamAdvectionCart()
  elseif Problem == "AdvectionCubeCart"
    @show Problem
    Param = ParamAdvectionCubeCart()
  elseif Problem == "WarmBubble2DXCart"
    @show Problem
    Param = ParamWarmBubble2DXCart()
  else
    @show "False Problem",Problem  
  end
end

