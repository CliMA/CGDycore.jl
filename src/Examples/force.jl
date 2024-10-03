abstract type AbstractForcing end

Base.@kwdef struct HeldSuarezForcing <: AbstractForcing end

function (Force::HeldSuarezForcing)(Param,Phys)
  function local_force(F,U,p,lat)
    FT = eltype(U)
    Sigma = p / Phys.p0
    SigmaPowKappa = fast_powGPU(Sigma,Phys.kappa)
    height_factor = max(FT(0), (Sigma - Param.sigma_b) / (FT(1) - Param.sigma_b))
    coslat = cos(lat)
    sinlat = sin(lat)
    F[2] += -(Param.k_f * height_factor) * U[2]
    F[3] += -(Param.k_f * height_factor) * U[3]
    if Sigma < FT(0.7)
      kT = Param.k_a + (Param.k_s - Param.k_a) * height_factor * coslat * coslat * coslat * coslat
    else
      kT = FT(0)
    end  
    Teq = (Param.T_equator - Param.DeltaT_y * sinlat * sinlat - 
      Param.DeltaTh_z * log(Sigma) * coslat * coslat) * SigmaPowKappa
    Teq = max(Param.T_min, Teq)
    DeltaT =  kT * (Phys.p0 * Sigma / (U[1] * Phys.Rd) - Teq)
    F[5]  += -U[1] * DeltaT / SigmaPowKappa
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



