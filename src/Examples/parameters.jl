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

Base.@kwdef struct ParamHaurwitzSphere
  ω = 7.848e-6 # Hz
  K = 7.848e-6 # Hz
  h0 = 8000 # m
  R = 4
  Ω = 2*pi/24.0/3600.0
end

Base.@kwdef struct ParamLinearBlob
  lat0 = 4.0*atan(1.0)
  lon0 = 2.0*atan(1.0)
  Width = 0.8
end

Base.@kwdef struct ParamBaroWaveDrySphere{FT}
  T0E::FT = 310.0
  T0P::FT = 240.0
  B::FT = 2.0
  K::FT = 3.0
  LapseRate::FT = 0.005
  U0::FT = -0.5
  PertR::FT = 1.0/6.0
  Up::FT = 1.0
  PertExpR::FT = 0.1
  PertLon::FT = pi/9.0
  PertLat::FT = 2.0 * pi / 9.0
  PertZ::FT = 15000.0
  NBr::FT = 1.e-2
  DeltaT::FT = 1
  ExpDist::FT = 5
  T0::FT = 300
  TEq::FT = 300
  T_init::FT  = 315
  lapse_rate::FT  = -0.008
  Deep::Bool = false
  pert::FT = 0.1
  uMax::FT = 1.0
  vMax::FT = 0.0
  DeltaT_y::FT = 0
  DeltaTh_z::FT = -5
  T_equator::FT = 315
  T_min::FT = 200
  sigma_b::FT = 7/10
  z_D::FT = 20.0e3
  #      Moist
  q_0::FT = 0.018                # Maximum specific humidity (default: 0.018)
  q_t::FT = 1.0e-12
  # Surface flux
  CMom::FT = 1.e-3
end

Base.@kwdef struct ParamBaroWaveHillDrySphere{FT}
  T0E::FT = 310.0
  T0P::FT = 240.0
  B::FT = 2.0
  K::FT = 3.0
  LapseRate::FT = 0.005
  U0::FT = -0.5
  PertR::FT = 1.0/6.0
  Up::FT = 0.0
  PertExpR::FT = 0.1
  PertLon::FT = pi/9.0
  PertLat::FT = 2.0 * pi / 9.0
  PertZ::FT = 15000.0
  NBr::FT = 1.e-2
  DeltaT::FT = 1
  ExpDist::FT = 5
  T0::FT = 300
  TEq::FT = 300
  T_init::FT  = 315
  lapse_rate::FT  = -0.008
  Deep::Bool = false
  pert::FT = 0.1
  uMax::FT = 1.0
  vMax::FT = 0.0
  DeltaT_y::FT = 0
  DeltaTh_z::FT = -5
  T_equator::FT = 315
  T_min::FT = 200
  sigma_b::FT = 7/10
  z_D::FT = 20.0e3
  #      Moist
  q_0::FT = 0.018                # Maximum specific humidity (default: 0.018)
  q_t::FT = 1.0e-12
  # Surface flux
  CMom::FT = 1.e-3
end

Base.@kwdef struct ParamHeldSuarezDrySphere{FT}
  day::FT = 3600.0 * 24.0
  k_a::FT= 1.0 / (40.0 * day)
  k_f::FT = 1.0 / day
  k_s::FT = 1.0 / (4.0 * day)
  DeltaT_y::FT = 60.0
  DeltaTh_z::FT = 10.0
  T_equator::FT = 315.0
  T_min::FT = 200.0
  sigma_b::FT = 7.0/10.0
  CM::FT  = 0.01 #0.0044
  CE::FT  = 0.0044
  CH::FT = 0.0044
  CTr::FT = 0.004
  p_pbl::FT = 85000.0
  p_strato::FT = 10000.0
  T_virt_surf::FT = 290.0
  T_min_ref::FT = 220.0
  H_t::FT = 8.e3
  q_0::FT = 0.018                # Maximum specific humidity (default: 0.018)
  q_t::FT = 1e-12
  T0E::FT = 310.0
  T0P::FT = 240.0
  B::FT = 2.0
  K::FT = 3.0
  LapseRate::FT = 0.005
  DeltaTS::FT = 29.0
  TSMin::FT = 271.0
  DeltaLat::FT = 26.0 * pi / 180.0
  uMax::FT = 0.0
  vMax::FT = 0.0
  CMom::FT = 1.e-3
  Deep::Bool = false
  T_Init::FT = 300.0
end

Base.@kwdef struct ParamHillSchaerCart
  Deep::Float64 = false
  NBr::Float64 = 1.e-2
  Th0::Float64 = 300.0
  uMax::Float64 = 10
  vMax::Float64 = 0
  TEq::Float64 = 300.0
  Stretch::Bool = false
end

Base.@kwdef struct ParamHillAgnesiXCart
  Deep::Bool = false
  NBr::Float64 = 1.e-2
  Th0::Float64 =300.0
  uMax::Float64 =10
  vMax::Float64 =0
  wMax::Float64 =0
  TEq::Float64 =300.0
  a::Float64  = 1000.0
  h::Float64  = 400.0
  xc::Float64  = 0.0
  Stretch::Bool = false
  CMom::Float64 = 1.e-3
end

Base.@kwdef struct ParamHillAgnesiYCart
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
  Th0::Float64 = 300.0
  uMax::Float64 = 0.0
  vMax::Float64 = 0
  wMax::Float64 = 0
  DeltaTh::Float64 = 2.0
  xC0::Float64 = 10000.0
  zC0::Float64 = 2000.0
  rC0::Float64 = 2000.0
end

Base.@kwdef struct ParamBryanFritschCart
  Th0::Float64 = 300.0
  uMax::Float64 = 0.0
  vMax::Float64 = 0
  wMax::Float64 = 0
  DeltaTh::Float64 = 2.0
  xC0::Float64 = 10000.0
  zC0::Float64 = 2000.0
  rC0::Float64 = 2000.0
end

Base.@kwdef struct ParamDensityCurrent2DXCart
  T0::Float64 = 300.0
  uMax::Float64 = 0.0
  vMax::Float64 = 0.0
  wMax::Float64 = 0.0
  DeltaT::Float64  = -15.0
  xC0::Float64  = 0.0
  zC0::Float64  = 2000.0
  xrC0::Float64  = 4000.0
  zrC0::Float64  = 2000.0
end  

Base.@kwdef struct ParamTestGradient
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
  TimeDependent::Bool = true
  hMax::Float64 = 0.95
  b::Float64 = 5.0
  lon1::Float64 = 5.0e0 / 6.0e0 * pi
  lat1::Float64 = 0.0e0
  lon2::Float64 = 7.0e0 / 6.0e0 * pi
  lat2::Float64 = 0.0e0
  EndTime::Float64 = 5.0
  FacVel::Float64 = 10.0
  StreamFun::Bool = false
end  
Base.@kwdef struct ParamAdvectionSphereDCMIP{FT}
  xC::FT = 0.0
  H::FT = H
  R_t::FT = RadEarth / 2.0
  Z_t::FT = 1000.0
  z_c::FT = 5.0e3
  p_top::FT = p_top
  T_0::FT = T_0
  ScaleHeight::FT = ScaleHeight
  b::FT = 0.2
  Lon_c1::FT = 150.0/360.0 * 2 * pi
  Lon_c2::FT = 210.0/360.0 * 2 * pi
  Lat_c::FT = 0.0
  tau::FT = tau
  omega_0::FT = omega_0
  TimeDependent::Bool = true
end

Base.@kwdef struct ParamAdvectionCubeCart
  StreamFun::Bool = false
  uMax::Float64 = 1.0
  vMax::Float64 = 1.0
  x1::Float64 = 399.0
  x2::Float64 = 601.0
  y1::Float64 = 399.0
  y2::Float64 = 601.0
end  

Base.@kwdef struct LimAdvectionCart{FT}
  xmin::FT = -2π              # domain x lower bound
  xmax::FT = 2π               # domain x upper bound
  ymin::FT = -2π              # domain y lower bound
  ymax::FT = 2π               # domain y upper bound
  zmin::FT = 0                # domain z lower bound
  zmax::FT = 4π               # domain z upper bound
  ρ₀::FT = 1.0                # air density
  D₄::FT = 0.0                # hyperdiffusion coefficient
  u0::FT = π / 2              # angular velocity
  r0::FT = (xmax - xmin) / 6  # bells radius
  end_time::FT = 2π           # simulation period in seconds
  centers1xC::FT = xmin + (xmax - xmin) / 4
  centers1yC::FT = ymin + (ymax - ymin) / 2
  centers1zC::FT = zmin + (zmax - zmin) / 2
  centers2xC::FT = xmin + 3 * (xmax - xmin) / 4
  centers2yC::FT = ymin + (ymax - ymin) / 2
  centers2zC::FT = zmin + (zmax - zmin) / 2
end

Base.@kwdef struct ParamAdvectionCubeRotCart
  StreamFun::Bool = false
  uMax::Float64 = 1.0
  vMax::Float64 = 0.0
  xC::Float64 = 500.0
  zC::Float64 = 500.0
  x1::Float64 = 299.0
  x2::Float64 = 501.0
  z1::Float64 = 299.0
  z2::Float64 = 501.0
  EndTime::Float64 = 1000.0
  H::Float64 = 1000.0
end

Base.@kwdef struct ParamAdvectionCart
  xC = 0.0
  H = H
end  

Base.@kwdef struct ParamSchaerSphericalSphere{FT}
  TEq::FT  = 300.0
  X::FT = 166.7
  H::FT = 20000
  uEq::FT = 20
end  


function Parameters(FT,Problem::String)
  if Problem == "BaroWaveDrySphere" || Problem == "BaroWaveDrySphereOro" || Problem == "BaroWaveMoistSphere"
    @show Problem
    Param = ParamBaroWaveDrySphere{FT}()
  elseif Problem == "BaroWaveHillDrySphere" || Problem == "BaroWaveHillMoistSphere"
    @show Problem
    Param = ParamBaroWaveHillDrySphere{FT}()
  elseif Problem == "SchaerSphericalSphere"
    @show Problem
    Param = ParamSchaerSphericalSphere{FT}()
  elseif Problem == "GalewskiSphere"
    @show Problem
    Param = ParamGalewskiSphere()
  elseif Problem == "HaurwitzSphere"
    @show Problem
    Param = ParamHaurwitzSphere()
  elseif Problem == "LinearBlob"
    @show Problem
    Param = ParamLinearBlob()
  elseif Problem == "HeldSuarezDrySphere" || Problem == "HeldSuarezDrySphereOro" || 
    Problem == "HeldSuarezMoistSphere" || Problem == "HeldSuarezMoistSphereOro"
    @show Problem
    Param = ParamHeldSuarezDrySphere{FT}()
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
    Param = ParamAdvectionSphereDCMIP{FT}()
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
  elseif Problem == "AdvectionCubeRotCart"
    @show Problem
    Param = ParamAdvectionCubeRotCart()
  elseif Problem == "LimAdvectionCart"
    @show Problem
    Param = LimAdvectionCart{FT}()  
  elseif Problem == "WarmBubble2DXCart"
    @show Problem
    Param = ParamWarmBubble2DXCart()
  elseif Problem == "HeldSuarezDrySphere" || Problem == "HeldSuarezDrySphereOro" || 
    Problem == "HeldSuarezMoistSphere" || Problem == "HeldSuarezMoistSphereOro"
    @show Problem
    Param = ParamHeldSuarezDrySphere{FT}()
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
  elseif Problem == "AdvectionDCMIP"
    @show Problem
    Param = ParamAdvectionSphereDCMIP{FT}()
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
  elseif Problem == "AdvectionCubeRotCart"
    @show Problem
    Param = ParamAdvectionCubeRotCart()
  elseif Problem == "LimAdvectionCart"
    @show Problem
    Param = LimAdvectionCart{FT}()  
  elseif Problem == "WarmBubble2DXCart"
    @show Problem
    Param = ParamWarmBubble2DXCart()
  elseif Problem == "BryanFritschCart"
    @show Problem
    Param = ParamBryanFritschCart()
  elseif Problem == "DensityCurrent2DXCart"
    @show Problem
    Param = ParamDensityCurrent2DXCart()
  elseif Problem == "TestGradient"
    @show Problem
    Param = ParamTestGradient()
  else
    @show "False Problem",Problem  
  end
end

