abstract type EquationType end
    
struct ShallowWater  <: EquationType end
struct Advection  <: EquationType end
struct CompressibleShallow  <: EquationType end
struct CompressibleDeep  <: EquationType  end

abstract type State end
struct ShallowWaterState  <: State  end
struct Dry  <: State  end
struct DryDG  <: State  end
struct Moist  <: State end
struct DryTotalEnergy  <: State  end
struct DryInternalEnergy  <: State  end
struct MoistInternalEnergy  <: State  end
struct IceInternalEnergy  <: State  end
  
function (::ShallowWaterState)(Phys)
  @inline function Pressure(Thermo,U,wL,wR,z)
    FT = eltype(U)
    p = FT(0.5) * Phys.Grav * U[5]^2
    T = FT(0.0)
    PotT = U[5] / U[1]
    Thermo[1] = p
    Thermo[2] = T
    Thermo[3] = PotT
  end
  return Pressure
end 

function (::DryDG)(Phys)
  @inline function Pressure(RhoTh)
    FT = eltype(RhoTh)
    p = Phys.p0 * fast_powGPU(Phys.Rd * RhoTh / Phys.p0, FT(1) / (FT(1) - Phys.kappa))
  end
  @inline function dPresdRhoTh(RhoTh)
    dpdRhoTh = Phys.Rd * (Phys.Rd * RhoTh / Phys.p0)^(Phys.kappa / (eltype(RhoTh)(1) - Phys.kappa))
    return dpdRhoTh
  end
  @inline function dPresdRho()
    dpdRho = eltype(Phys.Rd)(0)
    return dpdRho
  end
  return Pressure,dPresdRhoTh,dPresdRho
end

function (::Dry)(Phys)
  @inline function Pressure(Thermo,U,wL,wR,z;T=300.0)
    FT = eltype(U)
    p = Phys.p0 * fast_powGPU(Phys.Rd * U[5] / Phys.p0, FT(1) / (FT(1) - Phys.kappa))
    TLoc = p / (Phys.Rd * U[1])
    PotT = (Phys.p0/p)^(Phys.Rd/Phys.Cpd)*TLoc
    Thermo[1] = p
    Thermo[2] = TLoc
    Thermo[3] = PotT
  end
  @inline function dPresdRhoTh(RhoTh)
    dpdRhoTh = Phys.Rd * (Phys.Rd * RhoTh / Phys.p0)^(Phys.kappa / (eltype(RhoTh)(1) - Phys.kappa))
    return dpdRhoTh
  end  
  @inline function dPresdRho()
    dpdRho = eltype(Phys.Rd)(0)
    return dpdRho
  end
  return Pressure,dPresdRhoTh,dPresdRho
end 

function (::DryTotalEnergy)(Phys)
  @inline function Pressure(Thermo,U,wL,wR,z)
    FT = eltype(U)
    KE = FT(0.5) * (U[2]^2 + U[3]^2 + FT(0.5) * (wL^2 + wR^2))
    p = (Phys.Rd / Phys.Cvd) * (U[5] - U[1] * (KE + Phys.Grav * z))
    T = p / (Phys.Rd * U[1])
    PotT = (Phys.p0/p)^(Phys.Rd/Phys.Cpd)*T
    Thermo[1] = p
    Thermo[2] = T
    Thermo[3] = PotT
  end
  @inline function dPresdRhoE(RhoE)
    dpdRhoE = Phys.Rd / Phys.Cvd
    return dpdRhoE
  end  
  return Pressure,dPresdRhoE
end 

function (::DryInternalEnergy)(Phys)
  @inline function Pressure(Thermo,U,wL,wR,z;T=300.0)
    FT = eltype(U)
    TLoc = Thermodynamics.fTempIE(U[1],U[5],Phys)
    p = Phys.Rd * U[1] * TLoc
    if p < 0.0
      @show p,T
      @show Thermo
      @show U
      @show z
    end  
    PotT = (Phys.p0 / p)^(Phys.Rd / Phys.Cpd) * TLoc
    Thermo[1] = p
    Thermo[2] = TLoc
    Thermo[3] = PotT
  end
  @inline function dPresdRhoIE(RhoE)
    dpdRhoE = Phys.Rd / Phys.Cvd
    return dpdRhoE
  end
  @inline function dPresdRho()
    dpdRho = Phys.Rd * Phys.Cpd * Phys.T0 / Phys.Cvd 
    return dpdRho
  end
  return Pressure,dPresdRhoIE,dPresdRho
end

function (::Moist)(Phys,RhoPos,ThPos,RhoVPos,RhoCPos)
  @inline function Pressure(Thermo,U,wL,wR,z;T=300.0)
    FT = eltype(U)
    RhoV = U[RhoVPos]
    RhoC = U[RhoCPos]
    RhoD = U[RhoPos] - RhoV - RhoC
    Cpml = Phys.Cpd * RhoD + Phys.Cpv * RhoV + Phys.Cpl * RhoC
    Rm  = Phys.Rd * RhoD + Phys.Rv * RhoV
    kappaM = Rm / Cpml
    p = (Phys.Rd * U[ThPos] / Phys.p0^kappaM)^(FT(1) / (FT(1) - kappaM))
    TLoc = p / Rm
    PotT = (Phys.p0/p)^kappaM * TLoc
    Thermo[1] = p
    Thermo[2] = TLoc
    Thermo[3] = PotT
  end
  @inline function dPresdRhoTh(RhoTh)
    dpdRhoTh = Phys.Rd * (Phys.Rd * RhoTh / Phys.p0)^(Phys.kappa / (eltype(RhoTh)(1) - Phys.kappa))
    return dpdRhoTh
  end  
  @inline function dPresdRho()
    dpdRho = eltype(Phys.Rd)(0)
    return dpdRho
  end
  return Pressure,dPresdRhoTh,dPresdRho
end

function (::MoistInternalEnergy)(Phys,RhoPos,RhoIEPos,RhoTPos)
  @inline function Pressure(Thermo,U,wL,wR,z;T=300.0)
    FT = eltype(U)
    TLoc = T
    RhoV, RhoC, TLoc = SaturationAdjustmentIEW(U[RhoPos],U[RhoIEPos],U[RhoTPos],TLoc,Phys)
    p = ((U[RhoPos] - RhoV - RhoC) * Phys.Rd + RhoV * Phys.Rv) * TLoc
    PotT = (Phys.p0 / p)^(Phys.Rd / Phys.Cpd) * TLoc
    Thermo[1] = p
    Thermo[2] = TLoc
    Thermo[3] = PotT
    Thermo[5] = RhoV
    Thermo[6] = RhoC
  end
  @inline function dPresdRhoIE(RhoE)
    dpdRhoE = Phys.Rd / Phys.Cvd
    return dpdRhoE
  end
  @inline function dPresdRho()
    dpdRho = Phys.Rd * Phys.Cpd * Phys.T0 / Phys.Cvd 
    return dpdRho
  end
  return Pressure,dPresdRhoIE,dPresdRho
end


function (::IceInternalEnergy)(Phys,RhoPos,RhoIEPos,RhoTPos)
  @inline function Pressure(Thermo,U,wL,wR,z;T=300.0)
    FT = eltype(U)
    @inbounds RhoV, RhoC, RhoI, T = SaturationAdjustmentIEI(U[RhoPos],U[RhoIEPos],U[RhoTPos],T,Phys)
    p = ((U[RhoPos] - RhoV - RhoC - RhoI) * Phys.Rd + RhoV * Phys.Rv) * T
    PotT = (Phys.p0/p)^(Phys.Rd/Phys.Cpd)*T
    @inbounds Thermo[1] = p
    @inbounds Thermo[2] = T
    @inbounds Thermo[3] = PotT
    @inbounds Thermo[5] = RhoV
    @inbounds Thermo[6] = RhoC
    @inbounds Thermo[7] = RhoI
  end
  @inline function dPresdRhoIE(RhoE)
    dpdRhoE = Phys.Rd / Phys.Cvd
    return dpdRhoE
  end
  return Pressure,dPresdRhoIE
end

# we may be hitting a slow path:
# https://stackoverflow.com/questions/14687665/very-slow-stdpow-for-bases-very-close-to-1
@inline fast_powGPU(x::FT, y::FT) where {FT <: AbstractFloat} = exp(y * log(x))
