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
    w = 0
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
    Rho = pLoc / ((pLoc / p0)^kappa * Rd * ThLoc)
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
