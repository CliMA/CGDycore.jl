abstract type AbstractForcing end

struct HeldSuarezForcing <: AbstractForcing end

function (::HeldSuarezForcing)(Param,Phys)
  function local_force(U,p,lat)
    FT = eltype(U)
    Sigma = p / Phys.p0
    height_factor = max(FT(0), (Sigma - Param.sigma_b) / (FT(1) - Param.sigma_b))
    Fu = -(Param.k_f * height_factor) * U[2]
    Fv = -(Param.k_f * height_factor) * U[3]
    kT = (Sigma < 0.7) * (Param.k_a + (Param.k_s - Param.k_a) * height_factor * cos(lat)^4) 
    Teq = (Param.T_equator - Param.DeltaT_y * sin(lat)^2 - 
      Param.DeltaTh_z * log(Sigma) * cos(lat)^2) * Sigma^Phys.kappa
    Teq = max(Param.T_min, Teq)
    DeltaT =  kT * (Phys.p0 * Sigma / (U[1] * Phys.Rd) - Teq)
    FRhoTh  = -U[1] * DeltaT / Sigma^Phys.kappa
    return SVector{5}(FT(0),Fu,Fv,FT(0),FRhoTh)  
  end  
  return local_force
end

struct NoForcing <: AbstractForcing end
function (::NoForcing)(Param,Phys)
  function local_force(U,p,lat)
  end
  return local_force
end




