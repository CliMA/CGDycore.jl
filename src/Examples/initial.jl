abstract type Example end

Base.@kwdef struct LinearGravityExample <: Example end

function (profile::LinearGravityExample)(Param,Phys)
  @inline function local_profile(x,time)
    binv = (1 + (x[1] - Param.xc)^2 / Param.A^2)
    b = Param.b0 * sin(pi * x[3] / Param.H) / binv
    p = eltype(x)(0)
    u = eltype(x)(0)
    v = eltype(x)(0)
    w = eltype(x)(0)

    return (p,u,v,w,b)
  end
  return local_profile
end

Base.@kwdef struct InertiaGravityExample <: Example end

function (profile::InertiaGravityExample)(Param,Phys)
  @inline function local_profile(x,time)
    FT = eltype(x)
    ThPert = Param.DeltaTh * sin(pi * x[3] / Param.H) / (FT(1) + ((x[1] - Param.xC) / Param.a)^2)
    S = Param.NBr * Param.NBr / Phys.Grav
    ThBG = Param.Th0 * exp(x[3] *  S)
    Th =  ThBG + ThPert
    pLoc = Phys.p0 * (FT(1) - Phys.Grav / (Phys.Cpd * Param.Th0 * S) * 
      (FT(1) - exp(-S * x[3])))^(Phys.Cpd / Phys.Rd)
    Rho = pLoc / ((pLoc / Phys.p0)^Phys.kappa * Phys.Rd * Th)
    u = Param.uMax
    v = eltype(x)(0)
    w = eltype(x)(0)
    T = pLoc / (Phys.Rd * Rho)
    E = Phys.Cvd * T + FT(0.5) * (u * u + v * v) + Phys.Grav * x[3]
    IE = Phys.Cvd * T
    qv = eltype(x)(0)
    qc = eltype(x)(0)

    return (Rho,u,v,w,Th,E,IE,qv,qc,ThBG)
  end

  return local_profile
end

"""
  BickleyJetExample <: Example

Implements the classic "Bickley jet" initial condition, widely used to study barotropic instability and nonlinear jet evolution in atmospheric and oceanic flows.

# Usage

  profile = BickleyJetExample()
  initial_conditions = profile(Param, Phys)

# Arguments

- `Param`: Parameter object containing jet and perturbation parameters (e.g., `L`, `Ly`, `l`, `k`, `ϵ`).
- `Phys`: Physical parameters (not directly used in this function).

# Returns

A function `local_profile(x, time)` that computes the initial conditions at position `x` and time `time`:
- `Rho`: Density (set to 1.0).
- `u`: Zonal velocity, consisting of the jet profile and perturbations.
- `v`: Meridional velocity, consisting of perturbations.
- `w`: Vertical velocity (set to 0.0).
- `Th`: Potential temperature (set to 1.0).

# Physical Background

The Bickley jet is a model for a narrow, east-west (zonal) jet centered at y=0, defined by the velocity profile `U = sech(y)^2`.  
A streamfunction `Ψ = -tanh(y)` describes the main flow structure.  
Sinusoidal tracer fields and localized vortical perturbations (Gaussian-modulated cosines) are added to trigger barotropic instability and jet meandering.  
Velocity perturbations are computed from the streamfunction derivatives, resulting in small-scale turbulence superimposed on the jet.

This setup is commonly used to investigate the nonlinear evolution of jets, the formation of coherent structures, and the transition to turbulence (see Oceananigans.jl example and Beron-Vera et al. 2003).

# References

- Beron-Vera, F. J., Olascoaga, M. J., & Brown, M. G. (2003).  
  "The Nonlinear Evolution of Bickley Jets". Journal of Physical Oceanography, 33(10), 2173–2188.  
  [Link](https://journals.ametsoc.org/view/journals/phoc/33/10/1520-0485_2003_033_2173_tneobu_2.0.co_2.xml)

"""
Base.@kwdef struct BickleyJetExample <: Example end

function (profile::BickleyJetExample)(Param,Phys)
  @inline function local_profile(x,time)
    # Definition of the "Bickley jet": a sech(y)^2 jet with sinusoidal tracer
    Ψ = - tanh(x[2])
    U = sech(x[2])^2

    # A sinusoidal tracer
    C = sin(2π * x[2] / Param.L)

    # Slightly off-center vortical perturbations
    ψ̃ = exp(-(x[2] + Param.l/10)^2 / 2 * Param.l^2) * cos(Param.k * x[1]) * cos(Param.k * x[2])

    # Vortical velocity fields (ũ, ṽ) = (-∂_y, +∂_x) ψ̃
    ũ = + ψ̃ * (Param.k * tan(Param.k * x[2]) + x[2] / Param.l^2)
    ṽ = - ψ̃ * Param.k * tan(Param.k * x[1])

    # total initial conditions
    u = U + Param.ϵ * ũ
    v = Param.ϵ * ṽ
    w = 0.0
    Rho = sin(2π * x[2] / Param.Ly)
    Rho = 1.0
    Th = 1.0
    return (Rho,u,v,w,Th)
  end
  return local_profile
end

"""
  ModonCollisionExample <: Example

Implements the initial condition for a modon collision on the sphere, following the setup in 
McRae & Dritschel (2018), "The modon test: A barotropic test case for numerical methods on the sphere".

# Usage

    profile = ModonCollisionExample()
    initial_conditions = profile(Param, Phys)

# Arguments

- `Param`: Parameter object containing modon parameters (e.g., modon centers, radii, amplitudes).
- `Phys`: Physical parameters (e.g., gravity, planetary radius).

# Returns

A function `local_profile(x, time)` that computes the initial conditions at position `x` and time `time`:
- `Rho`: Layer depth (height field).
- `u`: Zonal velocity.
- `v`: Meridional velocity.
- `w`: Vertical velocity (set to 0.0).
- `Th`: Passive tracer or potential temperature (set to 1.0).

# Reference

- McRae, A. T. T., & Dritschel, D. G. (2018). 
  "The modon test: A barotropic test case for numerical methods on the sphere." 
  Geoscientific Model Development, 11(2), 645–655. 
  [Link](https://doi.org/10.5194/gmd-11-645-2018)
- Lin, S.-J., Chen, J.-H., Harris, L. M., & Zhou, L. (2017). 
  "Colliding modons on the sphere: A barotropic test case for numerical methods." 
  Geoscientific Model Development, 10(10), 3801–3817.
  [Link](https://doi.org/10.1002/2017MS000965)
"""
Base.@kwdef struct ModonCollisionExample <: Example end

function (profile::ModonCollisionExample)(Param, Phys)
    @inline function local_profile(x, time)
        FT = eltype(x)
        (lon, lat, r) = Grids.cart2sphere(x[1], x[2], x[3])
        R = Phys.RadEarth

        # Modon centers (in radians)
        lonC1, latC1 = Param.lonC1, Param.latC1
        lonC2, latC2 = Param.lonC2, Param.latC2

        r1 = Grids.SizeGreatCircle(lon, lat, lonC1, latC1) * R
        r2 = Grids.SizeGreatCircle(lon, lat, lonC2, latC2) * R

        # Modon amplitude profiles (Lin et al. 2017, Eq. 2.2)
        M1 = Param.u0 * exp(-(r1^2 / Param.r0^2))
        M2 = Param.u0 * exp(-(r2^2 / Param.r0^2))

        Rho = Param.h0
        u = M1 - M2
        v = FT(0)
        w = FT(0)
        Th = FT(1)

        return (Rho, u, v, w, Th)
    end
    return local_profile
end

Base.@kwdef struct DivergentSphereExample <: Example end

function (profile::DivergentSphereExample)(Param,Phys)
  @inline function local_profile(x,time)
    FT = eltype(x)
    Rho = FT(1)
    Lon,Lat,R = Grids.cart2sphere(x[1],x[2],x[3])
    lonP = Lon - FT(2) * pi * time / Param.EndTime
    uS = FT(10) / Param.EndTime * sin(lonP) * sin(lonP) *
      sin(FT(2) * Lat) * cos(pi * time / Param.EndTime) + FT(2) * pi / Param.EndTime * cos(Lat)
    vS = FT(10) / Param.EndTime * sin(FT(2) * lonP) * cos(Lat) * cos(pi * time / Param.EndTime)
    w = FT(0)
    lon1 = Param.lon1
    lat1 = Param.lat1
    lon2 = Param.lon2
    lat2 = Param.lat2
    R = FT(1)
    r = FT(0.5) * R
    r1 = R * Grids.SizeGreatCircle(Lon,Lat,lon1,lat1)
    r2 = R * Grids.SizeGreatCircle(Lon,Lat,lon2,lat2)
    if r1 <= r && abs(Lon - lon1) >= r / (FT(6.0) * R)
      Tr = FT(1.0)
    elseif r2 <= r && abs(Lon - lon2) >= r / (FT(6.0) * R)
      Tr = FT(1.0)
    elseif r1 <= r && abs(Lon - lon1) < r / (FT(6.0) * R) && Lat - lat1 < FT(-5.0 / 12.0) * r / R
      Tr = FT(1.0)
    elseif r2 <= r && abs(Lon - lon2) < r / (FT(6.0) * R) && Lat - lat2 > FT(5.0 / 12.0) * r / R
      Tr = FT(1.0)
    else
      Tr = FT(.1)
    end
    return (Rho,uS,vS,w,Tr)
  end
  return local_profile
end

Base.@kwdef struct RotationalCartExample <: Example end

function (profile::RotationalCartExample)(Param,Phys)
    @inline function local_profile(x,time)
      FT = eltype(x)
      Rho = FT(1)
      u = Param.uMax
      v = Param.vMax
      w = sinpi(x[3] / Param.H) * cospi(time / Param.EndTime)
      w = FT(0)
      if x[1] >= Param.x1 && x[1] <= Param.x2 && 
         x[2] >= Param.y1 && x[2] <= Param.y2 && 
         x[3] >= Param.z1 && x[3] <= Param.z2
        Tr = FT(1)
      else
        Tr = FT(0)
      end
      return (Rho,u,v,w,Tr)
    end
    return local_profile
end

Base.@kwdef struct AdvectionSphereSpherical <: Example end

function (profile::AdvectionSphereSpherical)(Param,Phys)
  function local_profile(x,time)
    FT = eltype(x)
    (Lon,Lat,R) = Grids.cart2sphere(x[1],x[2],x[3])
    Z = max(R-Phys.RadEarth,FT(0))
    uS = Param.uMax * cos(Lat)
    vS = 0.0
    w = 0.0
    RhoZ = 1.0
    d = acos(sin(Param.lat0)*sin(Lat)+cos(Param.lat0)*cos(Lat)*cos(Lon-Param.lon0))
    if abs(d) <= Param.Width
      Tr1 = cos(pi*d/Param.Width/2)^2 + 1.0
    else
      Tr1 = 0.0 + 1.0
    end
#   Tr1 = 1.0
    return (RhoZ,uS,vS,w,Tr1)
  end
  return local_profile
end

#
Base.@kwdef struct AdvectionVelocity<: Example end

function (profile::AdvectionVelocity)(Param,Phys)
  function local_profile(x,time)
    FT = eltype(x)
    (Lon,Lat,R) = Grids.cart2sphere(x[1],x[2],x[3])
    Z = max(R-Phys.RadEarth,FT(0))
    uS = Param.uMax * cos(Lat)
    vS = 0.0
    w = 0.0
    RhoZ = 1.0
    d = acos(sin(Param.lat0)*sin(Lat)+cos(Param.lat0)*cos(Lat)*cos(Lon-Param.lon0))
    if abs(d) <= Param.Width
      uS = cos(pi*d/Param.Width/2)^2 
    else
      uS = 0.0 
    end
    vS = 0.0
    return (0.0,uS,vS,0.0)
  end
  return local_profile
end

Base.@kwdef struct AdvectionSphereDCMIP <: Example end

function (profile::AdvectionSphereDCMIP)(Param,Phys)
  @inline function local_profile(x,time)
    FT = eltype(x)
    (Lon,Lat,R) = Grids.cart2sphere(x[1],x[2],x[3])
    Z=max(R-Phys.RadEarth,FT(0))
    pZ = Phys.p0 * exp(-Z / Param.ScaleHeight)
    LonP = Lon - FT(2)* FT(pi) * time / Param.tau
    k = FT(10) * Phys.RadEarth / Param.tau
    ua = k * sin(LonP)^2 * sin(2 * Lat) * cos(FT(pi) * time / Param.tau) +
      FT(2) * FT(pi) * Phys.RadEarth / Param.tau * cos(Lat)
    ud = Param.omega_0 * Phys.RadEarth / Param.b / Param.p_top *
      cos(LonP) *
      cos(Lat)^2 *
      cos(2 * pi * time / Param.tau) *
      (-exp((pZ - Phys.p0) / Param.b / Param.p_top) + exp((Param.p_top - pZ) / Param.b / Param.p_top))
    uS = ua + ud
    vS = k * sin(FT(2) * LonP) * cos(Lat) * cos(FT(pi) * time / Param.tau)
    RhoZ = pZ / Phys.Rd / Param.T_0
    sp =
       FT(1) + exp((Param.p_top - Phys.p0) / Param.b / Param.p_top) - exp((pZ - Phys.p0) / Param.b / Param.p_top) -
         exp((Param.p_top - pZ) / Param.b / Param.p_top)
    omega = Param.omega_0 * sin(LonP) * cos(Lat) * cos(FT(2) * FT(pi) * time / Param.tau) * sp
    w = -omega / RhoZ / Phys.Grav
    zd = Z - Param.z_c
    # great circle distances
    rd1 = Phys.RadEarth * Grids.SizeGreatCircle(Param.Lon_c1,Param.Lat_c,Lon,Lat)
    rd2 = Phys.RadEarth * Grids.SizeGreatCircle(Param.Lon_c2,Param.Lat_c,Lon,Lat)
    d1 = min(FT(1), (rd1 / Param.R_t)^2 + (zd / Param.Z_t)^2)
    d2 = min(FT(1), (rd2 / Param.R_t)^2 + (zd / Param.Z_t)^2)
    Tr1 = FT(0.5) * (FT(1) + cos(FT(pi) * d1)) + FT(0.5) * (FT(1) + cos(FT(pi) * d2))
    Tr2 = FT(0.9) - FT(0.8) * Tr1^2
#   return (RhoZ,uS,vS,w,[Tr1,Tr2])
    return (RhoZ,uS,vS,w,Tr1)
  end
  return local_profile
end

Base.@kwdef struct LimAdvectionCartExample <: Example end

function (profile::LimAdvectionCartExample)(Param,Phys)
  @inline function local_profile(x,time)
    FT = eltype(x)
    Rho = FT(1)
    rd1 = sqrt((x[1] - Param.centers1xC)^2 +
      (x[2] - Param.centers1yC)^2 +
      (x[3] - Param.centers1zC)^2)
    rd2 = sqrt((x[1] - Param.centers2xC)^2 +
      (x[2] - Param.centers2yC)^2 +
      (x[3] - Param.centers2zC)^2)
    if rd1 <= Param.r0 && abs(x[1] - Param.centers1xC) >= Param.r0 / FT(6)
      Tr = FT(1.0)
    elseif rd2 <= Param.r0 && abs(x[1] - Param.centers2xC) >= Param.r0 / FT(6)
      Tr = FT(1.0)
    elseif rd1 <= Param.r0 &&
      abs(x[1] - Param.centers1xC) < Param.r0 / FT(6) &&
      (x[2] - Param.centers1yC) < -FT(5) * Param.r0 / FT(12)
      Tr = FT(1.0)
    elseif rd2 <= Param.r0 &&
      abs(x[1] - Param.centers2xC) < Param.r0 / FT(6) &&
      (x[2] - Param.centers2yC) > FT(5) * Param.r0 / FT(12)
      Tr = FT(1.0)
    else
      Tr = FT(0.1)
    end

    u = -Param.u0  #* x[2]  * cospi(time / Param.end_time)
    v = Param.u0  #* x[1] * cospi(time / Param.end_time)
    w = Param.u0 * sinpi(x[3] / Param.zmax) * cospi(time / Param.end_time)

    return (Rho,u,v,w,Tr)
  end
  return local_profile
end

Base.@kwdef struct BryanFritsch <: Example end

function (profile::BryanFritsch)(Param,Phys,ProfileBF)
  @inline function local_profile(x,time)
    FT = eltype(x)
    z = x[3]
    @views zP = ProfileBF[:,1]
    iz = 1000
    for i = 2:size(zP,1)
      if z <= zP[i]
        iz = i - 1
        break
      end
    end
    z_l = zP[iz]
    Rho_l = ProfileBF[iz,2]
    Theta_l = ProfileBF[iz,3]
    RhoV_l = ProfileBF[iz,4]
    RhoC_l = ProfileBF[iz,5]
    p_l = ProfileBF[iz,6]
    z_r = zP[iz+1]
    Rho_r = ProfileBF[iz+1,2]
    Theta_r = ProfileBF[iz+1,3]
    RhoV_r = ProfileBF[iz+1,4]
    RhoC_r = ProfileBF[iz+1,5]
    p_r = ProfileBF[iz+1,6]
    Rho = (Rho_r * (z - z_l) + Rho_l * (z_r - z)) / (z_r - z_l)
    Theta = (Theta_r * (z - z_l) + Theta_l * (z_r - z)) / (z_r - z_l)
    ThBGrd = Theta
    RhoV = (RhoV_r * (z - z_l) + RhoV_l * (z_r - z)) / (z_r - z_l)
    RhoC = (RhoC_r * (z - z_l) + RhoC_l * (z_r - z)) / (z_r - z_l)
    pLoc = (p_r * (z - z_l) + p_l * (z_r - z)) / (z_r - z_l)

    Rho, Th, qV::FT, qC::FT = PerturbMoistProfile(x, Rho, Rho*Theta, RhoV, RhoC, Phys, Param)
    u = FT(0)
    v = FT(0)
    w = FT(0)
    T = pLoc / (Rho * (Phys.Rd * (FT(1.0) - qV - qC) + Phys.Rv * qV))
    IE = Thermodynamics.InternalEnergyW(FT(1.0),qV,qC,T,Phys) 
    E = IE + FT(0.5) * (u * u + v * v) + Phys.Grav * z
    return (Rho,u,v,w,Th,E,IE,qV,qC,ThBGrd)
  end
  return local_profile
end

Base.@kwdef struct WarmBubbleCartExample <: Example end

function (profile::WarmBubbleCartExample)(Param,Phys)
  @inline function local_profile(x,time)
    FT = eltype(x)
    u = Param.uMax
    v = Param.vMax
    w = FT(0)
    z = x[3]
    Grav = Phys.Grav
    p0 = Phys.p0
    Rd = Phys.Rd
    Rv = Phys.Rv
    Cvd = Phys.Cvd
    kappa = Phys.kappa
    Th0 = Param.Th0
    DeltaTh = Param.DeltaTh
    xC0 = Param.xC0
    zC0 = Param.zC0
    rC0 = Param.rC0
    x3 = x[3]
    x1 = x[1]
    pLoc = p0 * (FT(1) - Grav * x3 * kappa / (Rd * Th0))^(FT(1) / kappa)
    rr = sqrt((x1 - xC0)^2 + (x3 - zC0)^2)
    Th = Th0
    if rr < rC0
      Th = Th + DeltaTh * cos(FT(0.5) * pi * rr /rC0)^2
    end
    Rho = pLoc / ((pLoc / p0)^kappa * Rd * Th)
    T = pLoc / (Rd * Rho)
    E = Cvd * T + FT(0.5) * (u * u + v * v) + Grav * z
    IE = Cvd * T 
    qV = FT(0)
    qC = FT(0)
    ThBGrd = Th0

    return (Rho,u,v,w,Th,E,IE,qV,qC,ThBGrd)
  end
  return local_profile
end

Base.@kwdef struct NeutralCartExample <: Example end

function (profile::NeutralCartExample)(Param,Phys)
  @inline function local_profile(x,time)
    FT = eltype(x)
    u = Param.uMax
    v = Param.vMax
    w = FT(0)
    z = x[3]
    Grav = Phys.Grav
    p0 = Phys.p0
    Rd = Phys.Rd
    Rv = Phys.Rv
    Cvd = Phys.Cvd
    kappa = Phys.kappa
    Th0 = Param.Th0
    pLoc = p0 * (FT(1) - Grav * z * kappa / (Rd * Th0))^(FT(1) / kappa)
    Th = Th0
    Rho = pLoc / ((pLoc / p0)^kappa * Rd * Th)
    T = pLoc / (Rd * Rho)
    E = Cvd * T + FT(0.5) * (u * u + v * v) + Grav * z
    IE = Cvd * T 
    qV = FT(0)
    qC = FT(0)
    ThBGrd = Th0

    return (Rho,u,v,w,Th,E,IE,qV,qC,ThBGrd)
  end
  return local_profile
end

Base.@kwdef struct StratifiedExample <: Example end

function (profile::StratifiedExample)(Param,Phys)
  @inline function local_profile(x,time)
    FT = eltype(x)
    z = x[3]
    u = Param.uMax
    v = Param.vMax
    w = FT(0)
    NBr = Param.NBr
    Grav = Phys.Grav
    p0 = Phys.p0
    Cpd = Phys.Cpd
    Cvd = Phys.Cvd
    Rd = Phys.Rd
    kappa = Phys.kappa
    Th0 = Param.Th0
    S = NBr * NBr / Grav
    Th = Th0 * exp(z *  S)
    pLoc = p0 * (FT(1) - Grav / (Cpd * Th0 * S) * (FT(1) - exp(-S * z))).^(Cpd / Rd)
    Rho = pLoc / ((pLoc / p0)^kappa * Rd * Th)
    T = pLoc / (Rd * Rho)
    E = Cvd * T + FT(0.5) * (u * u + v * v) + Grav * z
    IE = Cvd * T 
    qv = FT(0)
    qc = FT(0)

    return (Rho,u,v,w,Th,E,IE,qv,qc,Th)
  end
  return local_profile
end

Base.@kwdef struct StratifiedSphereExample <: Example end

function (profile::StratifiedSphereExample)(Param,Phys)
  @inline function local_profile(x,time)
    FT = eltype(x)
    (lon,lat,r)= Grids.cart2sphere(x[1],x[2],x[3])
    z = max(r - Phys.RadEarth, FT(0))
    u = Param.uMax
    v = Param.vMax
    w = FT(0)
    NBr = Param.NBr
    Grav = Phys.Grav
    p0 = Phys.p0
    Cpd = Phys.Cpd
    Cvd = Phys.Cvd
    Rd = Phys.Rd
    kappa = Phys.kappa
    Th0 = Param.Th0
    S = NBr * NBr / Grav
    Th = Th0 * exp(z *  S)
    pLoc = p0 * (FT(1) - Grav / (Cpd * Th0 * S) * (FT(1) - exp(-S * z))).^(Cpd / Rd)
    Rho = pLoc / ((pLoc / p0)^kappa * Rd * Th)
    T = pLoc / (Rd * Rho)
    E = Cvd * T + FT(0.5) * (u * u + v * v) + Grav * z
    IE = Cvd * T 
    qv = FT(0)
    qc = FT(0)

    return (Rho,u,v,w,Th,E,IE,qv,qc,Th)
  end
  return local_profile
end

"""
  GalewskyExample <: Example

Implements the Galewsky barotropic instability initial condition, a standard test case for shallow water models on the sphere.

# Usage

  profile = GalewskyExample()
  initial_conditions = profile(Param, Phys)

# Arguments

- `Param`: Parameter object containing Galewsky test parameters (e.g., `H0G`, `hH`, `alphaG`, `betaG`, `lat0G`, `lat1G`, `uM`, `eN`).
- `Phys`: Physical parameters (e.g., `Grav`, `Omega`, `RadEarth`).

# Returns

A function `local_profile(x, time)` that computes the initial conditions at position `x` and time `time`:
- `Rho`: Layer thickness (height field).
- `u`: Zonal velocity.
- `v`: Meridional velocity (set to 0).
- `w`: Vertical velocity (set to 0).
- `Th`: Potential temperature (set to 1).

# Physical Background

The Galewsky test case is designed to study barotropic instability and the nonlinear evolution of a midlatitude jet on the sphere. The initial zonal velocity is confined between two latitudes and is in geostrophic balance with the height field. A localized height perturbation is added to trigger instability. This setup is widely used to benchmark shallow water models and study the development of jets and eddies (see Galewsky et al. 2004).

# References

- Galewsky, J., Scott, R. K., & Polvani, L. M. (2004).  
  "An initial-value problem for testing numerical models of the global shallow-water equations." Tellus A: Dynamic Meteorology and Oceanography, 56(5), 429–440.  
  [Link](https://doi.org/10.1111/j.1600-0870.2004.00049.x)

"""
Base.@kwdef struct GalewskyExample <: Example end

function (profile::GalewskyExample)(Param,Phys)
  @inline function local_profile(x,time)
    FT = eltype(x)
    Grav=Phys.Grav 
    Omega=Phys.Omega
    (lon,lat,r)= Grids.cart2sphere(x[1],x[2],x[3])
    r=Phys.RadEarth
    Rho=(Grav*Param.H0G-(simpson(-FT(0.5)*pi,lat,r,pi/FT(100.0),integrandG,Param,Phys)))/Grav +
      Param.hH*cos(lat)*exp(-((lon-pi)/Param.alphaG)^2)*exp(-((pi/FT(4.0)-lat)/Param.betaG)^2)
    Th = FT(1)   
    if (lat<=Param.lat0G) || (lat>=Param.lat1G)
      u=FT(0)
    else
      u=Param.uM/Param.eN*exp(FT(1.0)/((lat-Param.lat0G)*(lat-Param.lat1G)))
    end
#   Rho=Grav*Param.H0G
#   u = FT(0)
    v = FT(0)
    w = FT(0)
    Th = FT(1)

    return (Rho,u,v,w,Th)
  end
  return local_profile
end

"""
  HaurwitzExample <: Example

Implements the Haurwitz wave initial condition, a standard test case for the shallow water equations on the sphere, as described in Williamson et al. (1992).

# Usage

  profile = HaurwitzExample()
  initial_conditions = profile(Param, Phys)

# Arguments

- `Param`: Parameter object containing Haurwitz wave parameters:
  - `ω`: Angular velocity of the wave (rad/s)
  - `K`: Amplitude parameter (m/s)
  - `R`: Zonal wavenumber (integer, e.g., 4)
  - `h0`: Mean fluid depth (m)
- `Phys`: Physical parameters:
  - `Grav`: Gravitational acceleration (m/s²)
  - `Omega`: Planetary rotation rate (rad/s)
  - `RadEarth`: Planetary radius (m)

# Returns

A function `local_profile(x, time)` that computes the initial conditions at position `x` and time `time`:
- `Rho`: Layer thickness (height field).
- `u`: Zonal velocity.
- `v`: Meridional velocity.
- `w`: Vertical velocity (set to 0).
- `Th`: Potential temperature (set to 1).

# Physical Background

The Haurwitz wave is a classic analytic solution for the barotropic vorticity equation on the sphere, featuring a wavenumber-R pattern that propagates westward. It is widely used to test the accuracy and stability of global shallow water models.

# References

- Williamson, D. L., Drake, J. B., Hack, J. J., Jakob, R., & Swarztrauber, P. N. (1992).
  "A standard test set for numerical approximations to the shallow water equations in spherical geometry."
  Journal of Computational Physics, 102(1), 211–224.
  [Link](https://doi.org/10.1016/S0021-9991(05)80016-6)

"""
Base.@kwdef struct HaurwitzExample <: Example end

function (profile::HaurwitzExample)(Param,Phys)
  @inline function local_profile(x,time)
    FT = eltype(x)
    g = Phys.Grav
    Ω = Phys.Omega
    ω = Param.ω
    R = Param.R
    K = Param.K
    h0 = Param.h0
    (lon,lat,r) = Grids.cart2sphere(x[1],x[2],x[3])
    a = Phys.RadEarth
    Th = FT(1)
    aKcos = a * K * cos(lat)^(R-1)
    # eq. 143 from Williamson et al. 1992
    u = a * ω * cos(lat) +
      aKcos * (R * sin(lat)^2 - cos(lat)^2) * cos(R*lon)
    # eq. 144 from Williamson et al. 1992
    v = -aKcos * R * sin(lat) * sin(R * lon)
    w = FT(0)
    Th = FT(1)

    c2lat = cos(lat)^2
    cRlat = cos(lat)^R
    c2Rlat = cos(lat)^(2*R)
    # eq. 147 from Williamson et al. 1992
    A = ω / 2 * (2 * Ω + ω) * c2lat +
      K^2 / 4 * c2Rlat * ((R + 1) * c2lat + (2 * R^2 - R - 2) - 2 * R^2 / c2lat) 
    # eq. 148 from Williamson et al. 1992
    B =  2 * (Ω + ω) * K /((R + 1) * (R+2)) * cRlat *((R^2 + 2 * R + 2) - (R + 1)^2 * c2lat)
    # eq. 149 from Williamson et al. 1992
    C = K^2 / 4 * c2Rlat * ((R + 1) * c2lat - (R + 2))
    # eq. 146 from Williamson et al. 1992
    Rho = h0 + a^2 / g *(A + B * cos(R * lon) + C * cos(2 * R * lon))

    return (Rho,u,v,w,Th)
  end
  return local_profile
end

Base.@kwdef struct LinearBlob <: Example end

function (profile::LinearBlob)(Param,Phys)
  @inline function local_profile(x,time)
    FT = eltype(x)
    (lon,lat,r)= Grids.cart2sphere(x[1],x[2],x[3])
    d = acos(sin(Param.lat0)*sin(lat)+cos(Param.lat0)*cos(lat)*cos(lon-Param.lon0))
    if abs(d) <= Param.Width
      h = Param.H*cos(pi*d/Param.Width/2)^2 + 1.0
    else
      h = 0.0 + 1.0
    end
    u = FT(0)
    v = FT(0)
    w = FT(0)
    Th = FT(1)
    return (h,u,v,w,Th)
  end
  return local_profile
end

Base.@kwdef struct BaroWaveDryCart <: Example end

function (profile::BaroWaveDryCart)(Param,Phys)
  @inline function local_profile(x,time)
    FT = eltype(x)
    u0 = Param.u0
    Ly = Param.Ly
    b = Param.b
    up = Param.up
    xC = Param.xC
    yC = Param.yC
    Lp = Param.Lp
    p0 = Phys.p0
    Rd = Phys.Rd
    Cpd = Phys.Cpd
    Cvd = Phys.Cvd
    Grav = Phys.Grav
    z = x[3]
    eta = z2Eta(x[1],x[2],x[3],Param,Phys)
    Temperature = 300.0
    uC = -u0 * sin(pi * x[2] / Ly)^2 * log(eta) * exp(-(log(eta)/b)^2)
    uC += up * exp(-((x[1] - xC)^2 + (x[2] - yC)^2)/ Lp^2)
    vC = 0.0
    w = 0.0
    T = TemperatureBaroWaveCart(x[1],x[2],eta,Param,Phys)
    p = eta * p0
    Rho = p / (Rd * T)
    Th = T * (p0 / p)^(Rd / Cpd)
    E = Cvd * T + FT(0.5) * (uC * uC + vC * vC) + Grav * z
    return (Rho,uC,vC,w,Th,E)
  end
  return local_profile
end

Base.@kwdef struct BaroWaveDryExample <: Example end

function (profile::BaroWaveDryExample)(Param,Phys)
  @inline function local_profile(x,time)
    FT = eltype(x)
    (Lon,Lat,R)= Grids.cart2sphere(x[1],x[2],x[3])
    Z = max(R - Phys.RadEarth, FT(0))
    T0 = FT(0.5) * (Param.T0E + Param.T0P)
    ConstA = FT(1.0) / Param.LapseRate
    ConstB = (T0 - Param.T0P) / (T0 * Param.T0P)
    ConstC = FT(0.5) * (Param.K + FT(2.0)) * (Param.T0E - Param.T0P) / (Param.T0E * Param.T0P)
    ConstH = Phys.Rd * T0 / Phys.Grav
    ScaledZ = Z / (Param.B * ConstH)
    Tau1 = ConstA * Param.LapseRate / T0 * exp(Param.LapseRate / T0 * Z) +
    ConstB * (FT(1.0) - FT(2.0) * ScaledZ * ScaledZ) * exp(-ScaledZ * ScaledZ)
    Tau2 = ConstC * (FT(1.0) - FT(2.0) * ScaledZ * ScaledZ) * exp(-ScaledZ * ScaledZ)
    IntTau1 = ConstA * (exp(Param.LapseRate / T0 * Z) - FT(1.0)) +
    ConstB * Z * exp(-ScaledZ * ScaledZ)
    IntTau2 = ConstC * Z * exp(-ScaledZ * ScaledZ)
    if Param.Deep
      RRatio = R / Phys.RadEarth
    else
      RRatio = FT(1.0)
    end
    InteriorTerm = (RRatio * cos(Lat))^Param.K -
    Param.K / (Param.K + FT(2.0)) * (RRatio * cos(Lat))^(Param.K + FT(2.0))
    Temperature = FT(1.0) / (RRatio * RRatio) / (Tau1 - Tau2 * InteriorTerm)
    Pressure = Phys.p0 * exp(-Phys.Grav/Phys.Rd *
      (IntTau1 - IntTau2 * InteriorTerm))
    Rho = Pressure / (Phys.Rd * Temperature)
    Th = Temperature * (Phys.p0 / Pressure)^(Phys.Rd / Phys.Cpd)

    InteriorTermU = (RRatio * cos(Lat))^(Param.K - FT(1.0)) -
         (RRatio * cos(Lat))^(Param.K + FT(1.0))
    BigU = Phys.Grav / Phys.RadEarth * Param.K *
      IntTau2 * InteriorTermU * Temperature
    if Param.Deep
      RCosLat = R * cos(Lat)
    else
      RCosLat =Phys.RadEarth * cos(Lat)
    end
    OmegaRCosLat = Phys.Omega*RCosLat

    #                 if (dOmegaRCosLat * dOmegaRCosLat + dRCosLat * dBigU < 0.0) {
    #                         _EXCEPTIONT("Negative discriminant detected.")
    #                 }

    uS = -OmegaRCosLat + sqrt(OmegaRCosLat * OmegaRCosLat + RCosLat * BigU)
    vS = 0
    # Exponential perturbation
    GreatCircleR = acos(sin(Param.PertLat) * sin(Lat) +
      cos(Param.PertLat) * cos(Lat) * cos(Lon - Param.PertLon))

    GreatCircleR = GreatCircleR / Param.PertExpR

    # Tapered perturbation with height
    if Z < Param.PertZ
      PertTaper = FT(1.0) - FT(3.0) * Z * Z / (Param.PertZ * Param.PertZ) +
         FT(2.0) * Z * Z * Z / (Param.PertZ * Param.PertZ * Param.PertZ)
    else
      PertTaper = FT(0.0)
    end

    # Apply perturbation in zonal velocity
    if GreatCircleR < FT(1.0)
      uSPert = Param.Up * PertTaper * exp(-GreatCircleR * GreatCircleR)
    else
      uSPert = FT(0.0)
    end
    uS = uS + uSPert
    w = FT(0)

    qV = FT(0)
    qC = FT(0)
    IE = Thermodynamics.InternalEnergyW(FT(1),qV,qC,Temperature,Phys)
    E = IE + FT(0.5) * (uS * uS + vS * vS) + Phys.Grav * Z
    return (Rho,uS,vS,w,Th,E,IE,qV,qC)
  end
  return local_profile
end
Base.@kwdef struct BaroWaveMoistExample <: Example end

function (profile::BaroWaveMoistExample)(Param,Phys)
  @inline function local_profile(x,time)
    FT = eltype(x)
    (Lon,Lat,R)= Grids.cart2sphere(x[1],x[2],x[3])
    Z = max(R - Phys.RadEarth, FT(0))
    T0 = FT(0.5) * (Param.T0E + Param.T0P)
    ConstA = FT(1.0) / Param.LapseRate
    ConstB = (T0 - Param.T0P) / (T0 * Param.T0P)
    ConstC = FT(0.5) * (Param.K + FT(2.0)) * (Param.T0E - Param.T0P) / (Param.T0E * Param.T0P)
    ConstH = Phys.Rd * T0 / Phys.Grav
    ScaledZ = Z / (Param.B * ConstH)
    Tau1 = ConstA * Param.LapseRate / T0 * exp(Param.LapseRate / T0 * Z) +
    ConstB * (FT(1.0) - FT(2.0) * ScaledZ * ScaledZ) * exp(-ScaledZ * ScaledZ)
    Tau2 = ConstC * (FT(1.0) - FT(2.0) * ScaledZ * ScaledZ) * exp(-ScaledZ * ScaledZ)
    IntTau1 = ConstA * (exp(Param.LapseRate / T0 * Z) - FT(1.0)) +
    ConstB * Z * exp(-ScaledZ * ScaledZ)
    IntTau2 = ConstC * Z * exp(-ScaledZ * ScaledZ)
    if Param.Deep
      RRatio = R / Phys.RadEarth
    else
      RRatio = FT(1.0)
    end
    InteriorTerm = (RRatio * cos(Lat))^Param.K -
    Param.K / (Param.K + FT(2.0)) * (RRatio * cos(Lat))^(Param.K + FT(2.0))
    Temperature = FT(1.0) / (RRatio * RRatio) / (Tau1 - Tau2 * InteriorTerm)
    Pressure = Phys.p0 * exp(-Phys.Grav/Phys.Rd *
      (IntTau1 - IntTau2 * InteriorTerm))
    Rho = Pressure / (Phys.Rd * Temperature)
    Th = Temperature * (Phys.p0 / Pressure)^(Phys.Rd / Phys.Cpd)

    InteriorTermU = (RRatio * cos(Lat))^(Param.K - FT(1.0)) -
         (RRatio * cos(Lat))^(Param.K + FT(1.0))
    BigU = Phys.Grav / Phys.RadEarth * Param.K *
      IntTau2 * InteriorTermU * Temperature
    if Param.Deep
      RCosLat = R * cos(Lat)
    else
      RCosLat =Phys.RadEarth * cos(Lat)
    end
    OmegaRCosLat = Phys.Omega*RCosLat

    #                 if (dOmegaRCosLat * dOmegaRCosLat + dRCosLat * dBigU < 0.0) {
    #                         _EXCEPTIONT("Negative discriminant detected.")
    #                 }

    uS = -OmegaRCosLat + sqrt(OmegaRCosLat * OmegaRCosLat + RCosLat * BigU)
    vS = 0
    # Exponential perturbation
    GreatCircleR = acos(sin(Param.PertLat) * sin(Lat) +
      cos(Param.PertLat) * cos(Lat) * cos(Lon - Param.PertLon))

    GreatCircleR = GreatCircleR / Param.PertExpR

    # Tapered perturbation with height
    if Z < Param.PertZ
      PertTaper = FT(1.0) - FT(3.0) * Z * Z / (Param.PertZ * Param.PertZ) +
         FT(2.0) * Z * Z * Z / (Param.PertZ * Param.PertZ * Param.PertZ)
    else
      PertTaper = FT(0.0)
    end

    # Apply perturbation in zonal velocity
    if GreatCircleR < FT(1.0)
      uSPert = Param.Up * PertTaper * exp(-GreatCircleR * GreatCircleR)
    else
      uSPert = FT(0.0)
    end
    uS = uS + uSPert
    w = FT(0)

    eta = Pressure / Phys.p0
    eta_crit = Param.p_w / Phys.p0
    if eta > eta_crit
      qV = Param.q_0 * exp(-(Lat / Param.lat_w)^4) * exp(-((eta-FT(1)) * Param.p0 / Param.p_w)^2)
    else
      qV = Param.q_t  
    end 
    qC = FT(0)

    IE = Thermodynamics.InternalEnergyW(FT(1),qV,qC,Temperature,Phys)
    E = IE + FT(0.5) * (uS * uS + vS * vS) + Phys.Grav * Z
    return (Rho,uS,vS,w,Th,E,IE,qV,qC)
  end
  return local_profile
end

Base.@kwdef struct HeldSuarezDryExample <: Example end

function (::HeldSuarezDryExample)(Param,Phys)
  @inline function profile(x,time)
    FT = eltype(x)
    (Lon,Lat,R)= Grids.cart2sphere(x[1],x[2],x[3])
    Z = max(R - Phys.RadEarth, FT(0))
    T0 = FT(0.5) * (Param.T0E + Param.T0P)
    ConstA = FT(1.0) / Param.LapseRate
    ConstB = (T0 - Param.T0P) / (T0 * Param.T0P)
    ConstC = FT(0.5) * (Param.K + FT(2.0)) * (Param.T0E - Param.T0P) / (Param.T0E * Param.T0P)
    ConstH = Phys.Rd * T0 / Phys.Grav
    ScaledZ = Z / (Param.B * ConstH)
    Tau1 = ConstA * Param.LapseRate / T0 * exp(Param.LapseRate / T0 * Z) +
    ConstB * (FT(1.0) - FT(2.0) * ScaledZ * ScaledZ) * exp(-ScaledZ * ScaledZ)
    Tau2 = ConstC * (FT(1.0) - FT(2.0) * ScaledZ * ScaledZ) * exp(-ScaledZ * ScaledZ)
    IntTau1 = ConstA * (exp(Param.LapseRate / T0 * Z) - FT(1.0)) +
    ConstB * Z * exp(-ScaledZ * ScaledZ)
    IntTau2 = ConstC * Z * exp(-ScaledZ * ScaledZ)
    if Param.Deep
      RRatio = R / Phys.RadEarth
    else
      RRatio = FT(1.0)
    end
    InteriorTerm = (RRatio * cos(Lat))^Param.K -
    Param.K / (Param.K + FT(2.0)) * (RRatio * cos(Lat))^(Param.K + FT(2.0))
    Temperature = FT(1.0) / (RRatio * RRatio) / (Tau1 - Tau2 * InteriorTerm)
    Pressure = Phys.p0 * exp(-Phys.Grav/Phys.Rd *
      (IntTau1 - IntTau2 * InteriorTerm))
    Rho = Pressure / (Phys.Rd * Temperature)
    Th = Temperature * (Phys.p0 / Pressure)^(Phys.Rd / Phys.Cpd)

    InteriorTermU = (RRatio * cos(Lat))^(Param.K - FT(1.0)) -
         (RRatio * cos(Lat))^(Param.K + FT(1.0))
    BigU = Phys.Grav / Phys.RadEarth * Param.K *
      IntTau2 * InteriorTermU * Temperature
    if Param.Deep
      RCosLat = R * cos(Lat)
    else
      RCosLat =Phys.RadEarth * cos(Lat)
    end
    OmegaRCosLat = Phys.Omega*RCosLat

    #                 if (dOmegaRCosLat * dOmegaRCosLat + dRCosLat * dBigU < 0.0) {
    #                         _EXCEPTIONT("Negative discriminant detected.")
    #                 }

    uS = -OmegaRCosLat + sqrt(OmegaRCosLat * OmegaRCosLat + RCosLat * BigU)
    vS = 0
    w = FT(0)

    E = Phys.Cvd * Temperature + FT(0.5) * (uS * uS + vS * vS) + Phys.Grav * Z
    IE = Phys.Cvd * Temperature
    return (Rho,uS,vS,w,Th,E,IE)
  end
  @inline function Force(F,U,p,lat)
    FT = eltype(U)
    Sigma = p / Phys.p0
    SigmaPowKappa = fast_powGPU(Sigma,Phys.kappa)
    height_factor = max(FT(0), (Sigma - Param.sigma_b) / (FT(1) - Param.sigma_b))
    coslat = cos(lat)
    sinlat = sin(lat)
    F[2] += -(Param.k_f * height_factor) * U[2]
    F[3] += -(Param.k_f * height_factor) * U[3]
    kT = Param.k_a + (Param.k_s - Param.k_a) * height_factor * coslat * coslat * coslat * coslat
    Teq = (Param.T_equator - Param.DeltaT_y * sinlat * sinlat -
      Param.DeltaTh_z * log(Sigma) * coslat * coslat) * SigmaPowKappa
    Teq = max(Param.T_min, Teq)
    T = p / (U[1] * Phys.Rd)
    DeltaT =  kT * (T - Teq)
    F[5]  += - U[5] / T * DeltaT 
  end
  return profile,Force
end

Base.@kwdef struct HeldSuarezMoistExample <: Example end

function (::HeldSuarezMoistExample)(Param,Phys)
  @inline function profile(x,time)
    FT = eltype(x)
    (Lon,Lat,R)= Grids.cart2sphere(x[1],x[2],x[3])
    Z = max(R - Phys.RadEarth, FT(0))
    T0 = FT(0.5) * (Param.T0E + Param.T0P)
    ConstA = FT(1.0) / Param.LapseRate
    ConstB = (T0 - Param.T0P) / (T0 * Param.T0P)
    ConstC = FT(0.5) * (Param.K + FT(2.0)) * (Param.T0E - Param.T0P) / (Param.T0E * Param.T0P)
    ConstH = Phys.Rd * T0 / Phys.Grav
    ScaledZ = Z / (Param.B * ConstH)
    Tau1 = ConstA * Param.LapseRate / T0 * exp(Param.LapseRate / T0 * Z) +
    ConstB * (FT(1.0) - FT(2.0) * ScaledZ * ScaledZ) * exp(-ScaledZ * ScaledZ)
    Tau2 = ConstC * (FT(1.0) - FT(2.0) * ScaledZ * ScaledZ) * exp(-ScaledZ * ScaledZ)
    IntTau1 = ConstA * (exp(Param.LapseRate / T0 * Z) - FT(1.0)) +
    ConstB * Z * exp(-ScaledZ * ScaledZ)
    IntTau2 = ConstC * Z * exp(-ScaledZ * ScaledZ)
    if Param.Deep
      RRatio = R / Phys.RadEarth
    else
      RRatio = FT(1.0)
    end
    InteriorTerm = (RRatio * cos(Lat))^Param.K -
    Param.K / (Param.K + FT(2.0)) * (RRatio * cos(Lat))^(Param.K + FT(2.0))
    Temperature = FT(1.0) / (RRatio * RRatio) / (Tau1 - Tau2 * InteriorTerm)
    Pressure = Phys.p0 * exp(-Phys.Grav/Phys.Rd *
      (IntTau1 - IntTau2 * InteriorTerm))
    Rho = Pressure / (Phys.Rd * Temperature)
    Th = Temperature * (Phys.p0 / Pressure)^(Phys.Rd / Phys.Cpd)

    InteriorTermU = (RRatio * cos(Lat))^(Param.K - FT(1.0)) -
         (RRatio * cos(Lat))^(Param.K + FT(1.0))
    BigU = Phys.Grav / Phys.RadEarth * Param.K *
      IntTau2 * InteriorTermU * Temperature
    if Param.Deep
      RCosLat = R * cos(Lat)
    else
      RCosLat =Phys.RadEarth * cos(Lat)
    end
    OmegaRCosLat = Phys.Omega*RCosLat

    #                 if (dOmegaRCosLat * dOmegaRCosLat + dRCosLat * dBigU < 0.0) {
    #                         _EXCEPTIONT("Negative discriminant detected.")
    #                 }

    uS = -OmegaRCosLat + sqrt(OmegaRCosLat * OmegaRCosLat + RCosLat * BigU)
    uS = FT(0)
    vS = FT(0)
    w = FT(0)
    qV = FT(0)
    qC = FT(0)
    IE = Thermodynamics.InternalEnergyW(FT(1.0),qV,qC,Temperature,Phys) 
    E = IE + FT(0.5) * (uS * uS + vS * vS) + Phys.Grav * Z

    return (Rho,uS,vS,w,Th,E,IE,qV,qC)
  end
  @inline function Force(F,U,p,lat)
    FT = eltype(U)
    Sigma = p / Phys.p0
    SigmaPowKappa = fast_powGPU(Sigma,Phys.kappa)
    height_factor = max(FT(0), (Sigma - Param.sigma_b) / (FT(1) - Param.sigma_b))
    coslat = cos(lat)
    sinlat = sin(lat)
    F[2] += -(Param.k_f * height_factor) * U[2]
    F[3] += -(Param.k_f * height_factor) * U[3]
    kT = Param.k_a + (Param.k_s - Param.k_a) * height_factor * coslat * coslat * coslat * coslat
    Teq = (Param.T_equator - Param.DeltaT_y * sinlat * sinlat -
      Param.DeltaTh_z * log(Sigma) * coslat * coslat) * SigmaPowKappa
    Teq = max(Param.T_min, Teq)
    T = p / (U[1] * Phys.Rd)
    DeltaT =  kT * (T - Teq)
    F[5]  += - U[5] / T * DeltaT 
  end
  return profile, Force
end

Base.@kwdef struct SchaerSphereExample <: Example end

function (profile::SchaerSphereExample)(Param,Phys)
  @inline function local_profile(x,time)
    FT = eltype(x)
    (Lon,Lat,R)= Grids.cart2sphere(x[1],x[2],x[3])
    z = max(FT(0), R - Phys.RadEarth / Param.X)
    uS = Param.uEq * cos(Lat)
    vS = FT(0.0)
    w = FT(0.0)
    pLoc = Phys.p0 * exp(-Param.uEq * Param.uEq / (2.0 * Phys.Rd * Param.TEq) * sin(Lat)^2 -
      Phys.Grav * z / (Phys.Rd * Param.TEq))
    Th = Param.TEq * (Phys.p0 / pLoc)^(Phys.Rd / Phys.Cpd)
    Rho = pLoc / (Phys.Rd * Param.TEq)
    return (Rho,uS,vS,w,Th)
  end
  return local_profile
end

Base.@kwdef struct GapSphereExample <: Example end

function (profile::GapSphereExample)(Param,Phys)
  @inline function local_profile(x,time)
    FT = eltype(x)
    (Lon,Lat,R)= Grids.cart2sphere(x[1],x[2],x[3])
    z = max(FT(0), R - Phys.RadEarth / Param.X)
    uS = Param.uEq * cos(Lat)
    vS = FT(0.0)
    w = FT(0.0)
    pLoc = Phys.p0 * exp(-Param.uEq * Param.uEq / (2.0 * Phys.Rd * Param.TEq) * sin(Lat)^2 -
      Phys.Grav * z / (Phys.Rd * Param.TEq))
    Th = Param.TEq * (Phys.p0 / pLoc)^(Phys.Rd / Phys.Cpd)
    Rho = pLoc / (Phys.Rd * Param.TEq)
    return (Rho,uS,vS,w,Th)
  end
  return local_profile
end

@inline function integrandG(tau,RadEarth,Param,Phys)
  FT = eltype(tau)
  f=FT(2.0)*Phys.Omega*sin(tau)
  if (tau<=Param.lat0G) || (tau>=Param.lat1G)
    uStart=FT(0.0)
  else
    uStart=Param.uM/Param.eN*exp(FT(1.0)/((tau-Param.lat0G)*(tau-Param.lat1G)))
  end
  if abs(tau)<FT(0.5)*pi
    intG=(RadEarth*f+uStart*tan(tau))*uStart
  else
    intG=FT(0.0)
  end
  return intG
end

@inline function simpson(x0,xN,r,dx,f,Param,Phys)
  FT = eltype(x0)
# n=floor(Int32,(xN-x0)/dx) + 1
  n=50
  h=(xN-x0)/n
  res=FT(0.5)*(f(x0,r,Param,Phys)+f(xN,r,Param,Phys))
  xi=x0
  for i=1:(n-1)
    xi=xi+h
    res=res+f(xi,r,Param,Phys)
  end
  xi=x0-FT(0.5)*h
  for i=1:n
    xi=xi+h
    res=res+FT(2.0) *f(xi,r,Param,Phys)
  end
  res=res*h/FT(3.0)
  return res
end

@inline function GeoPotentialDeviation(y,Param)

  u0 = Param.u0
  f0 = Param.f0
  beta0 = Param.beta0
  y0 = Param.y0
  Ly = Param.Ly

  0.5 * u0 *((f0 + beta0 * y0) * (y - 0.5 * Ly - 0.5 * Ly / pi * sin(2 * pi * y / Ly)) +
    0.5* beta0 * (y^2 - Ly *y / pi *sin(2 * pi * y / Ly) -
    0.5 * Ly^2 / pi^2 * cos(2 * pi *y / Ly) - Ly^2 / 3 - 0.5 *Ly^2 / pi^2)) 

end

@inline function GeoPotentialMean(eta,Param,Phys)
  T0 = Param.T0
  LapseRate = Param.LapseRate
  Grav = Phys.Grav
  Rd = Phys.Rd

  T0 * Grav / LapseRate * (1 - eta^(Rd * LapseRate / Grav))
end  

@inline function GeoPotentialBaroWaveCart(x,y,eta,Param,Phys)
  b = Param.b

  GeoM = GeoPotentialMean(eta,Param,Phys)
  GeoD = GeoPotentialDeviation(y,Param)

  G = GeoM + GeoD * log(eta) * exp(-(log(eta) / b)^2)
end

@inline function TemperatureMean(eta,Param,Phys)
  T0 = Param.T0
  LapseRate = Param.LapseRate
  Grav = Phys.Grav
  Rd = Phys.Rd

  T0 * eta^(Rd * LapseRate / Grav)
end

@inline function TemperatureBaroWaveCart(x,y,eta,Param,Phys)
  Rd = Phys.Rd
  b = Param.b

  TMean = TemperatureMean(eta,Param,Phys)
  GeoD = GeoPotentialDeviation(y,Param)

  T = TMean + GeoD / Rd *(2 / b^2 * log(eta)^2 - 1) * exp(-(log(eta) / b)^2)
end  

@inline function z2Eta(x,y,z,Param,Phys)

  Eta = 1.e-2
  for i = 1 : 10
    F = -Grav * z + GeoPotentialBaroWaveCart(x,y,Eta,Param,Phys)
    dFdEta = -Rd / Eta * TemperatureBaroWaveCart(x,y,Eta,Param,Phys)
    Eta = Eta - F / dFdEta
  end

  return Eta

end

# we may be hitting a slow path:
# https://stackoverflow.com/questions/14687665/very-slow-stdpow-for-bases-very-close-to-1
@inline fast_powGPU(x::FT, y::FT) where {FT <: AbstractFloat} = exp(y * log(x))
