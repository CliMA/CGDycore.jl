abstract type EquationType end
    
struct ShallowWater  <: EquationType end
struct Advection  <: EquationType end
struct CompressibleShallow  <: EquationType end
struct CompressibleDeep  <: EquationType  end

abstract type State end
struct ShallowWaterState  <: State  end
struct Dry  <: State  end
struct Moist  <: State end
  
function (::ShallowWaterState)(Phys)
  function Pressure(U)
    FT = eltype(U)
    p = FT(0.5) * Phys.Grav * U[5]^2
    return p
  end
  return Pressure
end 

function (::Dry)(Phys)
  function Pressure(U)
    FT = eltype(U)
    p = Phys.p0 * fast_powGPU(Phys.Rd * U[5] / Phys.p0, FT(1) / (FT(1) - Phys.kappa))
    return p
  end
  return Pressure
end 

function (::Moist)(Phys,RhoPos,ThPos,RhoVPos,RhoCPos)
  function Pressure(U)
    FT = eltype(U)
    RhoV = U[RhoVPos]
    RhoC = U[RhoCPos]
    RhoD = U[RhoPos] - RhoV - RhoC
    Cpml = Phys.Cpd * RhoD + Phys.Cpv * RhoV + Phys.Cpl * RhoC
    Rm  = Phys.Rd * RhoD + Phys.Rv * RhoV
    kappaM = Rm / Cpml
    p = (Phys.Rd * U[ThPos] / Phys.p0^kappaM)^(FT(1) / (FT(1) - kappaM))
    return p
  end
  return Pressure
end

# we may be hitting a slow path:
# https://stackoverflow.com/questions/14687665/very-slow-stdpow-for-bases-very-close-to-1
fast_powGPU(x::FT, y::FT) where {FT <: AbstractFloat} = exp(y * log(x))
