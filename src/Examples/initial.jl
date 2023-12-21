abstract type Example end

Base.@kwdef struct DivergentSphereExample <: Example end

function (profile::DivergentSphereExample)(Param,Phys)
  function local_profile(x,time)
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
    function local_profile(x,time)
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

Base.@kwdef struct AdvectionSphereDCMIP <: Example end

function (profile::AdvectionSphereDCMIP)(Param,Phys)
  function local_profile(x,time)
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
  function local_profile(x,time)
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

Base.@kwdef struct WarmBubbleCartExample <: Example end

function (profile::WarmBubbleCartExample)(Param,Phys)
  function local_profile(x,time)
    FT = eltype(x)
    u = Param.uMax
    v = Param.vMax
    w = FT(0)
    Grav = Phys.Grav
    p0 = Phys.p0
    Rd = Phys.Rd
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
    return (Rho,u,v,w,Th)
  end
  return local_profile
end

Base.@kwdef struct StratifiedExample <: Example end

function (profile::StratifiedExample)(Param,Phys)
  function local_profile(x,time)
    FT = eltype(x)
    z = x[3]
    u = Param.uMax
    v = Param.vMax
    w = FT(0)
    NBr = Param.NBr
    Grav = Phys.Grav
    p0 = Phys.p0
    Cpd = Phys.Cpd
    Rd = Phys.Rd
    kappa = Phys.kappa
    Th0 = Param.Th0
    S = NBr * NBr / Grav
    Th = Th0 * exp(z *  S)
    pLoc = p0 * (FT(1) - Grav / (Cpd * Th0 * S) * (FT(1) - exp(-S * z))).^(Cpd / Rd)
    Rho = pLoc / ((pLoc / p0)^kappa * Rd * Th)
    return (Rho,u,v,w,Th)
  end
  return local_profile
end

Base.@kwdef struct GalewskiExample <: Example end

function (profile::GalewskiExample)(Param,Phys)
  function local_profile(x,time)
    FT = eltype(x)
    Grav=Phys.Grav 
    Omega=Phys.Omega
    (lon,lat,r)= Grids.cart2sphere(x[1],x[2],x[3])
    r=Phys.RadEarth
    Rho=(Grav*Param.H0G-(simpson(-0.5*pi,lat,r,pi/100.0,integrandG,Param)))/Grav +
      Param.hH*cos(lat)*exp(-((lon-pi)/Param.alphaG)^2.0)*exp(-((pi/4.0-lat)/Param.betaG)^2.0)
    Th = FT(1)   
    if (lat<=Param.lat0G) || (lat>=Param.lat1G)
      u=FT(0)
    else
      u=Param.uM/Param.eN*exp(1.0/((lat-Param.lat0G)*(lat-Param.lat1G)))
    end
    v = FT(0)
    w = FT(0)

    return (Rho,u,v,w,Th)
  end
  return local_profile
end

Base.@kwdef struct BaroWaveExample <: Example end

function (profile::BaroWaveExample)(Param,Phys)
  function local_profile(x,time)
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
    return (Rho,uS,vS,w,Th)
  end
  return local_profile
end

Base.@kwdef struct HeldSuarezDryExample <: Example end

function (::HeldSuarezDryExample)(Param,Phys)
  function profile(x,time)
    FT = eltype(x)
    (Lon,Lat,R)= Grids.cart2sphere(x[1],x[2],x[3])
    z=max(R-Phys.RadEarth,FT(0));
    temp = Param.T_Init + Param.LapseRate * z #+ rand(FT) * FT(0.1) * (z < FT(5000))
    pres = Phys.p0 * (FT(1) + Param.LapseRate / Param.T_Init * z)^(-Phys.Grav / Phys.Rd / Param.LapseRate)
    Rho = pres / Phys.Rd / temp
    Th = temp * (Phys.p0 / pres)^Phys.kappa
    uS = FT(0)
    vS = FT(0)
    w = FT(0)
    return (Rho,uS,vS,w,Th)
  end
  function Force(U,p,lat)
    FT = eltype(U)
    Sigma = p / Phys.p0
    height_factor = max(FT(0), (Sigma - Param.sigma_b) / (FT(1) - Param.sigma_b))
    Fu = -(Param.k_f * height_factor) * U[2]
    Fv = -(Param.k_f * height_factor) * U[3]
    if Sigma < FT(0.7)
      kT = Param.k_a + (Param.k_s - Param.k_a) * height_factor * cos(lat) * cos(lat) * cos(lat) * cos(lat)
    else
      kT = FT(0)
    end
    Teq = (Param.T_equator - Param.DeltaT_y * sin(lat) * sin(lat) -
      Param.DeltaTh_z * log(Sigma) * cos(lat) * cos(lat)) * Sigma^Phys.kappa
    Teq = max(Param.T_min, Teq)
    DeltaT =  kT * (Phys.p0 * Sigma / (U[1] * Phys.Rd) - Teq)
    FRhoTh  = -U[1] * DeltaT / Sigma^Phys.kappa
    return FT(0),Fu,Fv,FT(0),FRhoTh
  end
  return profile,Force
end

Base.@kwdef struct HeldSuarezMoistExample <: Example end

function (profile::HeldSuarezMoistExample)(Param,Phys)
  function local_profile(x,time)
    FT = eltype(x)
    (Lon,Lat,R)= Grids.cart2sphere(x[1],x[2],x[3])
    z=max(R-Phys.RadEarth,FT(0));
    temp = Param.T_Init + Param.LapseRate * z #+ rand(FT) * FT(0.1) * (z < FT(5000))
    pres = Phys.p0 * (FT(1) + Param.LapseRate / Param.T_Init * z)^(-Phys.Grav / Phys.Rd / Param.LapseRate)
    Rho = pres / Phys.Rd / temp
    Th = temp * (Phys.p0 / pres)^Phys.kappa
    uS = FT(0)
    vS = FT(0)
    w = FT(0)
    if z <= FT(1000)
      qv = FT(1.e-2)
    else
      qv = FT(1.e-8)
    end  
    qv = FT(0)
    qc = FT(0)
    return (Rho,uS,vS,w,Th,qv,qc)
  end
  function Force(U,p,lat)
    FT = eltype(U)
    Sigma = p / Phys.p0
    height_factor = max(FT(0), (Sigma - Param.sigma_b) / (FT(1) - Param.sigma_b))
    Fu = -(Param.k_f * height_factor) * U[2]
    Fv = -(Param.k_f * height_factor) * U[3]
    if Sigma < FT(0.7)
      kT = Param.k_a + (Param.k_s - Param.k_a) * height_factor * cos(lat) * cos(lat) * cos(lat) * cos(lat)
    else
      kT = FT(0)
    end
    Teq = (Param.T_equator - Param.DeltaT_y * sin(lat) * sin(lat) -
      Param.DeltaTh_z * log(Sigma) * cos(lat) * cos(lat)) * Sigma^Phys.kappa
    Teq = max(Param.T_min, Teq)
    DeltaT =  kT * (Phys.p0 * Sigma / (U[1] * Phys.Rd) - Teq)
    FRhoTh  = -U[1] * DeltaT / Sigma^Phys.kappa
    return FT(0),Fu,Fv,FT(0),FRhoTh
  end
  function Eddy(uStar,p,dz)
    K = Param.CE * uStar * dz / 2
    if p < Param.p_pbl
      dpR = (Param.p_pbl - p) / Param.p_strato
      K = K * exp(-dpR * dpR)
    end
    return K
  end
  function TSurf(x)
    FT = eltype(x)
    (Lon,Lat,R)= Grids.cart2sphere(x[1],x[2],x[3])
    TS = Param.DeltaTS * exp(-FT(0.5) * Lat^2 / Param.DeltaLat^2) + Param.TSMin
  end
  return local_profile,Force,Eddy,TSurf
end
