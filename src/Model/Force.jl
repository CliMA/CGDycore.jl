abstract type AbstractForcing end

Base.@kwdef struct HeldSuarezForcing <: AbstractForcing end

function (Force::HeldSuarezForcing)(Param,Phys)
  function local_force(U,p,lat)
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
  return local_force
end

struct NoForcing <: AbstractForcing end
function (::NoForcing)(Param,Phys)
  function local_force(U,p,lat)
    return FT(0),FT(0),FT(0),FT(0),FT(0)
  end
  return local_force
end



