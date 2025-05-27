abstract type NonConservativeFlux end

Base.@kwdef struct BuoyancyFlux <: NonConservativeFlux end

function (::BuoyancyFlux)(RhoPos,GPPos)
  @inline function Flux(VL,VR,AuxL,AuxR)
    GPL = AuxL[GPPos]
    GPR = AuxR[GPPos]
    RhoL = VL[RhoPos]
    RhoR = VR[RhoPos]
    BF = eltype(VL)(0.5) * (RhoL + RhoR) * (GPR - GPL)
  end
  return Flux
end

abstract type Flux end

Base.@kwdef struct EulerFlux <: Flux end

function (::EulerFlux)(RhoPos,wPos,ThPos,pPos)
  @inline function Flux(flux,V,Aux)
    p = Aux[pPos]
    w = V[wPos] / V[RhoPos]
    Th = V[ThPos] / V[RhoPos]

    flux[RhoPos] = V[uPos]
    flux[wPos] = V[uPos] * w + p
    flux[ThPos] = V[uPos] * Th

  end
  return Flux
end

Base.@kwdef struct AccousticFlux <: Flux end

function (::AccousticFlux)(pPos,wPos,cS)
  @inline function Flux(flux,V,Aux)
    p = V[pPos]
    w = V[wPos]

    flux[pPos] = cS^2 * w
    flux[wPos] = p

  end
  return Flux
end

abstract type AverageFluxV end

Base.@kwdef struct KennedyGruberGravV <: AverageFluxV end

function (::KennedyGruberGravV)(RhoPos,wPos,ThPos,pPos,GPPos)
  @inline function FluxNonLinAver!(flux,VL,VR,AuxL,AuxR,m_L,m_R)
    FT = eltype(flux)
    pL = AuxL[pPos]
    pR = AuxR[pPos]
    GPL = AuxL[GPPos]
    GPR = AuxR[GPPos]
    RhoL = VL[RhoPos]
    RhoR = VR[RhoPos]
    wL = VL[wPos] / RhoL
    ThL = VL[ThPos] / RhoL
    wR = VR[wPos] / RhoR
    ThR = VR[ThPos] / RhoR

    pAv = FT(0.5) * ((pL + pR) + FT(0.5) * (RhoL + RhoR) * (GPR - GPL))
    wAv = FT(0.5) * (wL + wR)
    RhoAv = FT(0.5) * (RhoL + RhoR)
    ThAv = FT(0.5) * (ThL + ThR)
    mAv3 = FT(0.5) * (m_L + m_R)
    qHat = mAv3 * wAv
    flux[1] = RhoAv * qHat
    flux[2] = flux[1] * wAv + mAv3 * pAv
    flux[3] = flux[1] * ThAv

  end
  return FluxNonLinAver!
end

Base.@kwdef struct KennedyGruberAccousticV <: AverageFluxV end

function (::KennedyGruberAccousticV)(pPos,wPos,cS)
  @inline function FluxNonLinAver!(flux,VL,VR,AuxL,AuxR,m_L,m_R)
    FT = eltype(flux)
    pL = VL[pPos]
    pR = VR[pPos]
    wL = VL[wPos]
    wR = VR[wPos]

    pAv = FT(0.5) * (pL + pR) 
    wAv = FT(0.5) * (wL + wR)
    mAv3 = FT(0.5) * (m_L + m_R)
    qHat = mAv3 * wAv
    flux[1] = cS^2 * qHat
    flux[2] = mAv3 * pAv
  end
  return FluxNonLinAver!
end


abstract type RiemannSolverV end

Base.@kwdef struct RiemannLMARSV <: RiemannSolverV end

function (::RiemannLMARSV)(Param,Phys,RhoPos,wPos,ThPos,pPos)
  @inline function RiemannByLMARSNonLin!(F,VLL,VRR,AuxL,AuxR)

    FT = eltype(F)
    cS = Param.cS
    pLL = AuxL[pPos]
    pRR = AuxR[pPos]
    RhoM = FT(0.5) * (VLL[RhoPos] + VRR[RhoPos])
    vLL = VLL[wPos] / VLL[RhoPos]
    vRR = VRR[wPos] / VRR[RhoPos]
    pM = FT(0.5) * (pLL + pRR) - FT(0.5) * cS * RhoM * (vRR - vLL)
    vM = FT(0.5) * (vRR + vLL) - FT(1.0) /(FT(2.0) * cS) * (pRR - pLL) / RhoM
    if vM > FT(0)
      F[RhoPos] = vM * VLL[RhoPos]
      F[wPos] = vM * VLL[wPos] + pM
      F[ThPos] = vM * VLL[ThPos]
    else
      F[RhoPos] = vM * VRR[RhoPos]
      F[wPos] = vM * VRR[wPos] + pM
      F[ThPos] = vM * VRR[ThPos]
    end
  end
  return RiemannByLMARSNonLin!
end

Base.@kwdef struct RiemannAccousticV <: RiemannSolverV end

function (::RiemannAccousticV)(Param,Phys,pPos,wPos,cS)
  @inline function RiemannSolver(F,VLL,VRR,AuxL,AuxR)

    pLL = VLL[pPos]
    pRR = VRR[pPos]
    vLL = VLL[wPos]
    vRR = VRR[wPos]
    pM = 0.5 * (pLL + pRR) - 0.5 * cS * (vRR - vLL)
    vM = 0.5 * (vRR + vLL) - 0.5 / cS * (pRR - pLL) 
    F[pPos] = cS^2 * vM
    F[wPos] = pM
  end
  return RiemannSolver
end

