abstract type Example end

Base.@kwdef struct RotationalCartExample <: Example end


function (profile::RotationalCartExample)(Param,Phys)
    function local_profile(x,time)
      FT = eltype(x)
      Rho = FT(1)
      u = Param.uMax
      v = Param.vMax
      w = sinpi(x[3] / Param.H) * cospi(time / Param.EndTime)
      if x[1] >= Param.x1 && x[1] <= Param.x2 && x[3] >= Param.z1 && x[3] <= Param.z2
        Tr = 1
      else
        Tr = 0
      end
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
    pLoc = p0 * (1 - Grav * x3 * kappa / (Rd * Th0))^(1 / kappa)
    rr = sqrt((x1 - xC0)^2 + (x3 - zC0)^2)
    Th = Th0
    if rr < rC0
      Th = Th + DeltaTh * cos(0.5 * pi * rr /rC0)^2
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


Base.@kwdef struct DCMIPAdvectionExample <: Example end


function (profile::DCMIPAdvectionExample)(Param,Phys)
  function local_profile(x,time)
    FT = eltype(x)
    (Lon,Lat,R) = cart2sphere(x[1],x[2],x[3])
    Z=max(R-Phys.RadEarth,0.0)
    pZ = Phys.p0 * exp(-Z / Param.ScaleHeight)
    LonP = Lon - 2* pi * time / Param.tau
    k = 10 * Phys.RadEarth / Param.tau
    ua = k * sin(LonP)^2 * sin(2 * Lat) * cos(pi * time / Param.tau) +
    2 * pi * Phys.RadEarth / Param.tau * cos(Lat)
    ud = Param.omega_0 * Phys.RadEarth / Param.b / Param.p_top *
      cos(LonP) *
      cos(Lat)^2 *
      cos(2 * pi * time / Param.tau) *
      (-exp((pZ - Phys.p0) / Param.b / Param.p_top) + exp((Param.p_top - pZ) / Param.b / Param.p_top))
    u = ua + ud
    v = k * sin(2 * LonP) * cos(Lat) * cos(pi * time / Param.tau)

    RhoZ = pZ / Phys.Rd / Param.T_0
    sp =
       1 + exp((Param.p_top - Phys.p0) / Param.b / Param.p_top) - exp((pZ - Phys.p0) / Param.b / Param.p_top) -
         exp((Param.p_top - pZ) / Param.b / Param.p_top)
    omega = Param.omega_0 * sin(LonP) * cos(Lat) * cos(2 * pi * time / Param.tau) * sp
    w = -omega / RhoZ / Phys.Grav

    zd = Z - Param.z_c
    # great circle distances
    rd1 = Phys.RadEarth * GreatCircle(Param.Lon_c1,Param.Lat_c,Lon,Lat)
    rd2 = Phys.RadEarth * GreatCircle(Param.Lon_c2,Param.Lat_c,Lon,Lat)
    d1 = min(1.0, (rd1 / Param.R_t)^2 + (zd / Param.Z_t)^2)
    d2 = min(1.0, (rd2 / Param.R_t)^2 + (zd / Param.Z_t)^2)
    Tr = 0.5 * (1 + cos(pi * d1)) + 0.5 * (1 + cos(pi * d2))
    return (RhoZ,u,v,w,Tr)
  end
  return local_profile
end

Base.@kwdef struct GalewskiExample <: Example end

function (profile::GalewskiExample)(Param,Phys)
  function local_profile(x,time)
    FT = eltype(x)
    Grav=Phys.Grav 
    Omega=Phys.Omega
    (lon,lat,r)=cart2sphere(x[1],x[2],x[3])
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
    (Lon,Lat,R)=cart2sphere(x[1],x[2],x[3])
    Z = max(R - Phys.RadEarth, 0)
    T0 = 0.5 * (Param.T0E + Param.T0P)
    ConstA = 1.0 / Param.LapseRate
    ConstB = (T0 - Param.T0P) / (T0 * Param.T0P)
    ConstC = 0.5 * (Param.K + 2.0) * (Param.T0E - Param.T0P) / (Param.T0E * Param.T0P)
    ConstH = Phys.Rd * T0 / Phys.Grav
    ScaledZ = Z / (Param.B * ConstH)
    Tau1 = ConstA * Param.LapseRate / T0 * exp(Param.LapseRate / T0 * Z) +
    ConstB * (1.0 - 2.0 * ScaledZ * ScaledZ) * exp(-ScaledZ * ScaledZ)
    Tau2 = ConstC * (1.0 - 2.0 * ScaledZ * ScaledZ) * exp(-ScaledZ * ScaledZ)
    IntTau1 = ConstA * (exp(Param.LapseRate / T0 * Z) - 1.0) +
    ConstB * Z * exp(-ScaledZ * ScaledZ)
    IntTau2 = ConstC * Z * exp(-ScaledZ * ScaledZ)
    if Param.Deep
      RRatio = R / Phys.RadEarth
    else
      RRatio = 1.0
    end
    InteriorTerm = (RRatio * cos(Lat))^Param.K -
    Param.K / (Param.K + 2.0) * (RRatio * cos(Lat))^(Param.K + 2.0)
    Temperature = 1.0 / (RRatio * RRatio) / (Tau1 - Tau2 * InteriorTerm)
    Pressure = Phys.p0 * exp(-Phys.Grav/Phys.Rd *
      (IntTau1 - IntTau2 * InteriorTerm))
    Rho = Pressure / (Phys.Rd * Temperature)
    Th = Temperature * (Phys.p0 / Pressure)^(Phys.Rd / Phys.Cpd)

    InteriorTermU = (RRatio * cos(Lat))^(Param.K - 1.0) -
         (RRatio * cos(Lat))^(Param.K + 1.0)
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
        PertTaper = 1.0 - 3.0 * Z * Z / (Param.PertZ * Param.PertZ) +
           2.0 * Z * Z * Z / (Param.PertZ * Param.PertZ * Param.PertZ)
      else
        PertTaper = 0.0
      end

      # Apply perturbation in zonal velocity
      if GreatCircleR < 1.0
        uSPert = Param.Up * PertTaper * exp(-GreatCircleR * GreatCircleR)
      else
        uSPert = 0.0
      end
      uS = uS + uSPert
      w = FT(0)

    return (Rho,uS,vS,w,Th)
  end
  return local_profile
end

