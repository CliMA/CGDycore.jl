abstract type Example end

Base.@kwdef struct RotationalCartExample <: Example end


function RotationalCart(profile::RotationalCartExample)
    function local_profile(x,time,Param,Phys)
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


function WarmBubbleCart(profile::WarmBubbleCartExample)
  function local_profile(x,time,Param,Phys)
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

