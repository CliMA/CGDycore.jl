################################################################################
# Two-point volume fluxes for the split-form DGSEM.
#
# A flux is a function
#   (VL, VR, AuxL, AuxR, mL, mR) -> SVector{5}
# with the primitive states VL, VR = (rho, u, v, w, theta), the auxiliary
# states AuxL, AuxR (pressure at pPos, geopotential at GPPos) and the metric
# rows mL, mR = J a^dir at both nodes. It returns the contravariant flux in
# the order (rho, rho u, rho v, rho w, rho theta). The volume kernels are
# independent of the flux formula.
#
# Gravity is discretized as the nonconservative product rhoM * [[phi]] and
# added to the pressure average, {{p}} + 1/2 rhoM [[phi]].
#
# Potential temperature fluxes after Artiano, Knoth, Spichtinger, Ranocha,
# arXiv:2509.10311: TEC (total energy conservative), EC (entropy
# conservative) and ETEC (entropy and total energy conservative). The density
# mean rhoM (ArithMean, LogMean, GammaMean) is used in the density flux and in
# the gravity term, as required for total energy conservation; LogMean is
# well-balanced for T = const, GammaMean for theta = const.
################################################################################

abstract type AverageFlux end
abstract type RiemannSolver end

# Formulation of the prognostic variables, decides which components are
# divided by rho when the volume kernels load the primitive state
abstract type Formulation end
struct ConservativeForm <: Formulation end      # (rho, rho v, rho theta)
struct VectorInvariantForm <: Formulation end   # (rho, v, rho theta)

Formulation(ModelType::String) =
  ModelType == "VectorInvariant" ? VectorInvariantForm() : ConservativeForm()

@inline PrimScale(::ConservativeForm, iv, invRho) = invRho
@inline PrimScale(::VectorInvariantForm, iv, invRho) = iv == 5 ? invRho : one(invRho)

@inline function contravariant_velocity(VL, VR, mL, mR)
  FT = eltype(VL)
  m1 = FT(0.5) * (mL[1] + mR[1])
  m2 = FT(0.5) * (mL[2] + mR[2])
  m3 = FT(0.5) * (mL[3] + mR[3])
  uAv = FT(0.5) * (VL[2] + VR[2])
  vAv = FT(0.5) * (VL[3] + VR[3])
  wAv = FT(0.5) * (VL[4] + VR[4])
  qHat = m1 * uAv + m2 * vAv + m3 * wAv
  return qHat, uAv, vAv, wAv, m1, m2, m3
end

# Momentum flux f_rho {{v}} + m pM (kinetic energy preserving)
@inline function flux_vector(fRho, fRhoTh, pM, uAv, vAv, wAv, m1, m2, m3)
  return SVector(fRho,
    fRho * uAv + m1 * pM,
    fRho * vAv + m2 * pM,
    fRho * wAv + m3 * pM,
    fRhoTh)
end

Base.@kwdef struct KennedyGruberGrav <: AverageFlux end

function (::KennedyGruberGrav)(::Val{pPos}, ::Val{GPPos}) where {pPos, GPPos}
  @inline function FluxNonLinAver(VL, VR, AuxL, AuxR, mL, mR)
    FT = eltype(VL)
    qHat, uAv, vAv, wAv, m1, m2, m3 = contravariant_velocity(VL, VR, mL, mR)
    RhoAv = FT(0.5) * (VL[1] + VR[1])
    pM = FT(0.5) * ((AuxL[pPos] + AuxR[pPos]) + RhoAv * (AuxR[GPPos] - AuxL[GPPos]))
    fRho = RhoAv * qHat
    fRhoTh = fRho * FT(0.5) * (VL[5] + VR[5])
    return flux_vector(fRho, fRhoTh, pM, uAv, vAv, wAv, m1, m2, m3)
  end
  return FluxNonLinAver
end

# TEC: f_rho = rhoM {{v}}, f_rhotheta = {{rho theta}}_gamma {{v}}
Base.@kwdef struct PotTempTEC <: AverageFlux end

function (::PotTempTEC)(::Val{pPos}, ::Val{GPPos}, Mean::DensityMean, Gamma) where {pPos, GPPos}
  @inline function FluxNonLinAver(VL, VR, AuxL, AuxR, mL, mR)
    FT = eltype(VL)
    qHat, uAv, vAv, wAv, m1, m2, m3 = contravariant_velocity(VL, VR, mL, mR)
    RhoM = density_mean(Mean, VL[1], VR[1], FT(Gamma))
    pM = FT(0.5) * ((AuxL[pPos] + AuxR[pPos]) + RhoM * (AuxR[GPPos] - AuxL[GPPos]))
    fRho = RhoM * qHat
    fRhoTh = stolarsky_mean(VL[1] * VL[5], VR[1] * VR[5], FT(Gamma)) * qHat
    return flux_vector(fRho, fRhoTh, pM, uAv, vAv, wAv, m1, m2, m3)
  end
  return FluxNonLinAver
end

# EC: f_rho = rhoM {{v}}, f_rhotheta = f_rho / {{1/theta}}_log
Base.@kwdef struct PotTempEC <: AverageFlux end

function (::PotTempEC)(::Val{pPos}, ::Val{GPPos}, Mean::DensityMean, Gamma) where {pPos, GPPos}
  @inline function FluxNonLinAver(VL, VR, AuxL, AuxR, mL, mR)
    FT = eltype(VL)
    qHat, uAv, vAv, wAv, m1, m2, m3 = contravariant_velocity(VL, VR, mL, mR)
    RhoM = density_mean(Mean, VL[1], VR[1], FT(Gamma))
    pM = FT(0.5) * ((AuxL[pPos] + AuxR[pPos]) + RhoM * (AuxR[GPPos] - AuxL[GPPos]))
    fRho = RhoM * qHat
    fRhoTh = fRho / ln_mean(FT(1) / VL[5], FT(1) / VR[5])
    return flux_vector(fRho, fRhoTh, pM, uAv, vAv, wAv, m1, m2, m3)
  end
  return FluxNonLinAver
end

# ETEC: f_rhotheta = {{rho theta}}_gamma {{v}}, f_rho = f_rhotheta {{1/theta}}_log,
# the density mean of the gravity term is fixed to rhoM = f_rho / {{v}}
Base.@kwdef struct PotTempETEC <: AverageFlux end

function (::PotTempETEC)(::Val{pPos}, ::Val{GPPos}, Gamma) where {pPos, GPPos}
  @inline function FluxNonLinAver(VL, VR, AuxL, AuxR, mL, mR)
    FT = eltype(VL)
    qHat, uAv, vAv, wAv, m1, m2, m3 = contravariant_velocity(VL, VR, mL, mR)
    RhoM = stolarsky_mean(VL[1] * VL[5], VR[1] * VR[5], FT(Gamma)) *
      ln_mean(FT(1) / VL[5], FT(1) / VR[5])
    pM = FT(0.5) * ((AuxL[pPos] + AuxR[pPos]) + RhoM * (AuxR[GPPos] - AuxL[GPPos]))
    fRho = RhoM * qHat
    fRhoTh = stolarsky_mean(VL[1] * VL[5], VR[1] * VR[5], FT(Gamma)) * qHat
    return flux_vector(fRho, fRhoTh, pM, uAv, vAv, wAv, m1, m2, m3)
  end
  return FluxNonLinAver
end
