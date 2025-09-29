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

Base.@kwdef struct BuoyancyFluxSlow <: NonConservativeFlux end

function (::BuoyancyFluxSlow)(RhoPos,GPPos)
  @inline function Flux(VL,VR,AuxL,AuxR)
    GPL = AuxL[GPPos]
    GPR = AuxR[GPPos]
    RhoL = VL[RhoPos]
    RhoR = VR[RhoPos]
    BF = 0.0
  end
  return Flux
end
abstract type Flux end

Base.@kwdef struct EulerFlux <: Flux end

function (::EulerFlux)(RhoPos,uPos,vPos,wPos,ThPos,pPos)
  @inline function Flux(flux,V,Aux)
    p = Aux[pPos]
    u = V[uPos] / V[RhoPos]
    v = V[vPos] / V[RhoPos]
    w = V[wPos] / V[RhoPos]
    Th = V[ThPos] / V[RhoPos]

    flux[1,RhoPos] = V[uPos]
    flux[1,uPos] = V[uPos] * u + p
    flux[1,vPos] = V[uPos] * v
    flux[1,wPos] = V[uPos] * w
    flux[1,ThPos] = V[uPos] * Th

    flux[2,RhoPos] = V[vPos]
    flux[2,uPos] = V[vPos] * u
    flux[2,vPos] = V[vPos] * v + p
    flux[2,wPos] = V[vPos] * w
    flux[2,ThPos] = V[vPos] * Th


    flux[3,RhoPos] = V[wPos]
    flux[3,uPos] = V[wPos] * u
    flux[3,vPos] = V[wPos] * v
    flux[3,wPos] = V[wPos] * w + p
    flux[3,ThPos] = V[wPos] * Th
  end
  return Flux
end

Base.@kwdef struct LinearizedEulerFlux <: Flux end

function (::LinearizedEulerFlux)(RhoPos,uPos,vPos,wPos,RhoThPos,dpdRhoThPos,ThPos)
  @inline function Flux(flux,V,Aux)
    p = Aux[dpdRhoThPos]*V[RhoThPos]
    Th = Aux[ThPos]

    flux[1,RhoPos] = V[uPos]
    flux[1,uPos] = p
    flux[1,vPos] = eltype(flux)(0)
    flux[1,wPos] = eltype(flux)(0)
    flux[1,RhoThPos] = V[uPos] * Th

    flux[2,RhoPos] = V[vPos]
    flux[2,uPos] = eltype(flux)(0)
    flux[2,vPos] = p
    flux[2,wPos] = eltype(flux)(0)
    flux[2,RhoThPos] = V[vPos] * Th


    flux[3,RhoPos] = V[wPos]
    flux[3,uPos] = eltype(flux)(0)
    flux[3,vPos] = eltype(flux)(0)
    flux[3,wPos] = p
    flux[3,RhoThPos] = V[wPos] * Th
  end
  return Flux
end

Base.@kwdef struct LinearBoussinesqFlux <: Flux end

function (::LinearBoussinesqFlux)(Param,pPos,uPos,vPos,wPos,bPos)
  @inline function Flux(flux,V,Aux)

    p = V[pPos]
    u = V[uPos]
    v = V[vPos]
    w = V[wPos]
    b = V[bPos]

    flux[1,pPos] = -Param.U * p - Param.cS^2 * u
    flux[1,uPos] = -Param.U * u - p
    flux[1,vPos] = eltype(flux)(0)
    flux[1,wPos] = -Param.U * w 
    flux[1,bPos] = eltype(flux)(0)

    flux[2,pPos] = -Param.cS^2 * v
    flux[2,uPos] = eltype(flux)(0)
    flux[2,vPos] = -p
    flux[2,wPos] = eltype(flux)(0)
    flux[2,bPos] = eltype(flux)(0)

    flux[3,pPos] = -Param.cS^2 * w
    flux[3,uPos] = -p
    flux[3,vPos] = eltype(flux)(0)
    flux[3,wPos] = eltype(flux)(0)
    flux[3,bPos] = eltype(flux)(0)

  end
  return Flux
end

abstract type AverageFlux end

Base.@kwdef struct KennedyGruber <: AverageFlux end

function (::KennedyGruber)(RhoPos,uPos,vPos,wPos,ThPos,pPos)
  @inline function FluxNonLinAver!(flux,VL,VR,AuxL,AuxR,m_L,m_R)
    FT = eltype(flux)
    pL = AuxL[pPos]
    pR = AuxR[pPos]
    RhoL = VL[RhoPos]
    RhoR = VR[RhoPos]
    uL = VL[uPos] / RhoL
    vL = VL[vPos] / RhoL
    wL = VL[wPos] / RhoL
    ThL = VL[ThPos] / RhoL
    uR = VR[uPos] / RhoR
    vR = VR[vPos] / RhoR
    wR = VR[wPos] / RhoR
    ThR = VR[ThPos] / RhoR

    pAv = FT(0.5) * (pL + pR)
    uAv = FT(0.5) * (uL + uR)
    vAv = FT(0.5) * (vL + vR)
    wAv = FT(0.5) * (wL + wR)
    RhoAv = FT(0.5) * (RhoL + RhoR)
    ThAv = FT(0.5) * (ThL + ThR)
    mAv1 = FT(0.5) * (m_L[1] + m_R[1])
    mAv2 = FT(0.5) * (m_L[2] + m_R[2])
    mAv3 = FT(0.5) * (m_L[3] + m_R[3])
    qHat = mAv1 * uAv + mAv2 * vAv + mAv3 * wAv
    flux[1] = RhoAv * qHat
    flux[2] = flux[1] * uAv + mAv1 * pAv
    flux[3] = flux[1] * vAv + mAv2 * pAv
    flux[4] = flux[1] * wAv + mAv3 * pAv
    flux[5] = flux[1] * ThAv
  end
  return FluxNonLinAver!
end

Base.@kwdef struct ArtianoGrav <: AverageFlux end

function (::ArtianoGrav)(RhoPos,uPos,vPos,wPos,ThPos,pPos,GPPos,Phys)
  @inline function FluxNonLinAver!(flux,VL,VR,AuxL,AuxR,m_L,m_R)
    FT = eltype(flux)

    RhoAv = FT(0.5) * (VL[RhoPos] + VR[RhoPos])
    pAv = FT(0.5) * ((AuxL[pPos] + AuxR[pPos]) + 
      RhoAv * (AuxR[GPPos] - AuxL[GPPos]))
    uAv = FT(0.5) * (VL[uPos] + VR[uPos])
    vAv = FT(0.5) * (VL[vPos] + VR[vPos])
    wAv = FT(0.5) * (VL[wPos] + VR[wPos])
    RhoThAv = stolarsky_mean(VL[ThPos], VR[ThPos], Phys.Gamma)
    RhoAv = FT(0.5) * (VL[RhoPos] + VR[RhoPos])
    ThAv = FT(0.5) * (VL[ThPos] + VR[ThPos])
    mAv1 = FT(0.5) * (m_L[1] + m_R[1])
    mAv2 = FT(0.5) * (m_L[2] + m_R[2])
    mAv3 = FT(0.5) * (m_L[3] + m_R[3])
    qHat = mAv1 * uAv + mAv2 * vAv + mAv3 * wAv
    flux[1] = RhoAv * qHat
    flux[2] = flux[1] * uAv + mAv1 * pAv
    flux[3] = flux[1] * vAv + mAv2 * pAv
    flux[4] = flux[1] * wAv + mAv3 * pAv
    flux[5] = RhoThAv * qHat

  end
  return FluxNonLinAver!
end

Base.@kwdef struct ArtianoExGrav <: AverageFlux end

function (::ArtianoExGrav)(RhoPos,uPos,vPos,wPos,ThPos,pPos,GPPos,Phys)
  @inline function FluxNonLinAver!(flux,VL,VR,AuxL,AuxR,m_L,m_R)
    FT = eltype(flux)

    RhoAv = FT(0.5) * (VL[RhoPos] + VR[RhoPos])
    ThAv = FT(0.5) * (VL[ThPos] + VR[ThPos])
    RhoThAv = FT(0.5) * (VL[RhoPos] * VL[ThPos] + VR[RhoPos] * VR[ThPos])
    pAv = FT(0.5) * (Phys.Cpd * RhoThAv * (AuxR[3] - AuxL[3]) +
      RhoAv * (AuxR[GPPos] - AuxL[GPPos]))
    uAv = FT(0.5) * (VL[uPos] + VR[uPos])
    vAv = FT(0.5) * (VL[vPos] + VR[vPos])
    wAv = FT(0.5) * (VL[wPos] + VR[wPos])
    mAv1 = FT(0.5) * (m_L[1] + m_R[1])
    mAv2 = FT(0.5) * (m_L[2] + m_R[2])
    mAv3 = FT(0.5) * (m_L[3] + m_R[3])
    qHat = mAv1 * uAv + mAv2 * vAv + mAv3 * wAv
    flux[1] = RhoAv * qHat
    flux[2] = flux[1] * uAv + mAv1 * pAv
    flux[3] = flux[1] * vAv + mAv2 * pAv
    flux[4] = flux[1] * wAv + mAv3 * pAv
    flux[5] = flux[1] * ThAv

  end
  return FluxNonLinAver!
end

Base.@kwdef struct KennedyGruberExPGrav <: AverageFlux end

function (::KennedyGruberExPGrav)(RhoPos,uPos,vPos,wPos,ThPos,pPos,GPPos,Phys)
  @inline function FluxNonLinAver!(flux,VL,VR,AuxL,AuxR,m_L,m_R)
    FT = eltype(flux)

    RhoAv = FT(0.5) * (VL[RhoPos] + VR[RhoPos])
    ThAv = FT(0.5) * (VL[ThPos] + VR[ThPos])
    pAv = FT(0.5) * (Phys.Cpd * VL[RhoPos] * VL[ThPos] * (AuxL[3] + AuxR[3]) +
      RhoAv * (AuxR[GPPos] - AuxL[GPPos]))
    uAv = FT(0.5) * (VL[uPos] + VR[uPos])
    vAv = FT(0.5) * (VL[vPos] + VR[vPos])
    wAv = FT(0.5) * (VL[wPos] + VR[wPos])
    mAv1 = FT(0.5) * (m_L[1] + m_R[1])
    mAv2 = FT(0.5) * (m_L[2] + m_R[2])
    mAv3 = FT(0.5) * (m_L[3] + m_R[3])
    qHat = mAv1 * uAv + mAv2 * vAv + mAv3 * wAv
    flux[1] = RhoAv * qHat
    flux[2] = flux[1] * uAv + mAv1 * pAv
    flux[3] = flux[1] * vAv + mAv2 * pAv
    flux[4] = flux[1] * wAv + mAv3 * pAv
    flux[5] = flux[1] * ThAv

  end
  return FluxNonLinAver!
end


Base.@kwdef struct ArtianoExner <: AverageFlux end

function (::ArtianoExner)(RhoPos,uPos,vPos,wPos,ThPos,pPos,GPPos,Phys)
  @inline function FluxNonLinAver!(flux,VL,VR,AuxL,AuxR,m_L,m_R)
    FT = eltype(flux)
    ExpL = AuxL[3]
    ExpR = AuxR[3]
    RhoThetaL = VL[ThPos]*VL[RhoPos]
    RhoThetaR = VR[ThPos]*VR[RhoPos]
    uAv = FT(0.5) * (VL[uPos] + VR[uPos])
    vAv = FT(0.5) * (VL[vPos] + VR[vPos])
    wAv = FT(0.5) * (VL[wPos] + VR[wPos])
    RhoAv = FT(0.5) * (VL[RhoPos] + VR[RhoPos])
    ThAv = FT(0.5) * (VL[ThPos] + VR[ThPos])
    pAv = FT(0.5) * (Phys.Cpd*RhoAv*ThAv*(ExpR - ExpL) + FT(0.5) * (VL[RhoPos] + VR[RhoPos])* (AuxR[GPPos] - AuxL[GPPos]))
    #pAv = FT(0.5) * ((AuxL[pPos] + AuxR[pPos]) + FT(0.5) * (VL[RhoPos] + VR[RhoPos]) * (AuxR[GPPos] - AuxL[GPPos]))
    mAv1 = FT(0.5) * (m_L[1] + m_R[1])
    mAv2 = FT(0.5) * (m_L[2] + m_R[2])
    mAv3 = FT(0.5) * (m_L[3] + m_R[3])
    qHat = mAv1 * uAv + mAv2 * vAv + mAv3 * wAv
    flux[1] = RhoAv * qHat
    flux[2] = flux[1] * uAv + mAv1 * pAv
    flux[3] = flux[1] * vAv + mAv2 * pAv
    flux[4] = flux[1] * wAv + mAv3 * pAv
    flux[5] = flux[1] * ThAv	
  end
  return FluxNonLinAver!
end


Base.@kwdef struct KennedyGruberGrav <: AverageFlux end

function (::KennedyGruberGrav)(RhoPos,uPos,vPos,wPos,ThPos,pPos,GPPos)
  @inline function FluxNonLinAver!(flux,VL,VR,AuxL,AuxR,m_L,m_R)
    FT = eltype(flux)
    pAv = FT(0.5) * ((AuxL[pPos] + AuxR[pPos]) +
			FT(0.5) * (VL[RhoPos] + VR[RhoPos]) * (AuxR[GPPos] - AuxL[GPPos]))
    uAv = FT(0.5) * (VL[uPos] + VR[uPos])
    vAv = FT(0.5) * (VL[vPos] + VR[vPos])
    wAv = FT(0.5) * (VL[wPos] + VR[wPos])
    RhoAv = FT(0.5) * (VL[RhoPos] + VR[RhoPos])
    ThAv = FT(0.5) * (VL[ThPos] + VR[ThPos])
    mAv1 = FT(0.5) * (m_L[1] + m_R[1])
    mAv2 = FT(0.5) * (m_L[2] + m_R[2])
    mAv3 = FT(0.5) * (m_L[3] + m_R[3])
    qHat = mAv1 * uAv + mAv2 * vAv + mAv3 * wAv
    flux[1] = RhoAv * qHat
    flux[2] = flux[1] * uAv + mAv1 * pAv
    flux[3] = flux[1] * vAv + mAv2 * pAv
    flux[4] = flux[1] * wAv + mAv3 * pAv
    flux[5] = flux[1] * ThAv

  end
  return FluxNonLinAver!
end

Base.@kwdef struct KennedyGruberGravFast <: AverageFlux end

function (::KennedyGruberGravFast)(RhoPos,uPos,vPos,wPos,ThPos,pPos,GPPos)
  @inline function FluxNonLinAver!(flux,VL,VR,AuxL,AuxR,m_L,m_R)
    FT = eltype(flux)
    pL = AuxL[pPos]
    pR = AuxR[pPos]
    GPL = AuxL[GPPos]
    GPR = AuxR[GPPos]
    RhoL = VL[RhoPos]
    RhoR = VR[RhoPos]
    uL = VL[uPos] / RhoL
    vL = VL[vPos] / RhoL
    wL = VL[wPos] / RhoL
    ThL = VL[ThPos] / RhoL
    uR = VR[uPos] / RhoR
    vR = VR[vPos] / RhoR
    wR = VR[wPos] / RhoR
    ThR = VR[ThPos] / RhoR

    pAv = FT(0.5) * ((pL + pR) + FT(0.5) * (RhoL + RhoR) * (GPR - GPL))
    uAv = FT(0.5) * (uL + uR)
    vAv = FT(0.5) * (vL + vR)
    wAv = FT(0.5) * (wL + wR)
    RhoAv = FT(0.5) * (RhoL + RhoR)
    ThAv = FT(0.5) * (ThL + ThR)
    mAv1 = FT(0.5) * (m_L[1] + m_R[1])
    mAv2 = FT(0.5) * (m_L[2] + m_R[2])
    mAv3 = FT(0.5) * (m_L[3] + m_R[3])
    qHat = mAv1 * uAv + mAv2 * vAv + mAv3 * wAv
    flux[1] = RhoAv * qHat
    flux[2] = mAv1 * pAv
    flux[3] = mAv2 * pAv
    flux[4] = mAv3 * pAv
    flux[5] = flux[1] * ThAv

  end
  return FluxNonLinAver!
end


Base.@kwdef struct KennedyGruberGravSlow <: AverageFlux end

function (::KennedyGruberGravSlow)(RhoPos,uPos,vPos,wPos,ThPos,pPos,GPPos)
  @inline function FluxNonLinAver!(flux,VL,VR,AuxL,AuxR,m_L,m_R)
    FT = eltype(flux)
    RhoL = VL[RhoPos]
    RhoR = VR[RhoPos]
    uL = VL[uPos] / RhoL
    vL = VL[vPos] / RhoL
    wL = VL[wPos] / RhoL
    uR = VR[uPos] / RhoR
    vR = VR[vPos] / RhoR
    wR = VR[wPos] / RhoR

    uAv = FT(0.5) * (uL + uR)
    vAv = FT(0.5) * (vL + vR)
    wAv = FT(0.5) * (wL + wR)
    RhoAv = FT(0.5) * (RhoL + RhoR)
    mAv1 = FT(0.5) * (m_L[1] + m_R[1])
    mAv2 = FT(0.5) * (m_L[2] + m_R[2])
    mAv3 = FT(0.5) * (m_L[3] + m_R[3])
    qHat = mAv1 * uAv + mAv2 * vAv + mAv3 * wAv
    temp = RhoAv * qHat
    flux[2] = temp * uAv
    flux[3] = temp * vAv
    flux[4] = temp * wAv
    flux[5] = 0.0
    flux[1] = 0.0

  end
  return FluxNonLinAver!
end

abstract type RiemannSolver end

Base.@kwdef struct RiemannLMARS <: RiemannSolver end
Base.@kwdef struct RiemannExnerLMARS <: RiemannSolver end
Base.@kwdef struct ArtianoEnergyStable <: RiemannSolver end
Base.@kwdef struct RiemannExLMARS <: RiemannSolver end
Base.@kwdef struct RiemannExPLMARS <: RiemannSolver end
Base.@kwdef struct RiemannLMARSFast <: RiemannSolver end
Base.@kwdef struct RiemannLMARSSlow <: RiemannSolver end

function (::RiemannLMARS)(Param,Phys,hPos,uPos,vPos,wPos,pPos)
  @inline function RiemannByLMARSNonLin!(F,VLL,VRR,AuxL,AuxR)

    FT = eltype(F)
    cS = Param.cS
    pLL = AuxL[pPos]
    pRR = AuxR[pPos]
    hM = FT(0.5) * (VLL[hPos] + VRR[hPos])
    vLL = VLL[uPos] / VLL[hPos]
    vRR = VRR[uPos] / VRR[hPos]
    pM = FT(0.5) * (pLL + pRR) - FT(0.5) * cS * hM * (vRR - vLL)
    vM = FT(0.5) * (vRR + vLL) - FT(1.0) /(FT(2.0) * cS) * (pRR - pLL) / hM
    if vM > FT(0)
      F[hPos] = vM * VLL[hPos]
      F[uPos] = vM * VLL[uPos] + pM
      F[vPos] = vM * VLL[vPos]
      F[wPos] = vM * VLL[wPos]
    else
      F[hPos] = vM * VRR[hPos]
      F[uPos] = vM * VRR[uPos] + pM
      F[vPos] = vM * VRR[vPos]
      F[wPos] = vM * VRR[wPos]
    end
  end
  return RiemannByLMARSNonLin!
end

function (::RiemannLMARS)(Param,Phys,RhoPos,uPos,vPos,wPos,ThPos,pPos)
  @inline function RiemannByLMARSNonLin!(F,VLL,VRR,AuxL,AuxR,Normal)
    FT = eltype(F)
    cS = Param.cS
    pLL = AuxL[pPos]
    pRR = AuxR[pPos]
    RhoM = FT(0.5) * (VLL[RhoPos] + VRR[RhoPos])
    vLL = (VLL[uPos] * Normal[1] + VLL[vPos] * Normal[2] + VLL[wPos] * Normal[3]) / VLL[RhoPos]
    vRR = (VRR[uPos] * Normal[1] + VRR[vPos] * Normal[2] + VRR[wPos] * Normal[3]) / VRR[RhoPos]
    pM = FT(0.5) * (pLL + pRR) - FT(0.5) * cS * RhoM * (vRR - vLL)
    vM = FT(0.5) * (vRR + vLL) - FT(1.0) /(FT(2.0) * cS) * (pRR - pLL) / RhoM
    if vM > FT(0)
      F[RhoPos] = vM * VLL[RhoPos]
      F[uPos] = vM * VLL[uPos] + Normal[1] * pM
      F[vPos] = vM * VLL[vPos] + Normal[2] * pM
      F[wPos] = vM * VLL[wPos] + Normal[3] * pM
      F[ThPos] = vM * VLL[ThPos]
    else
      F[RhoPos] = vM * VRR[RhoPos]
      F[uPos] = vM * VRR[uPos] + Normal[1] * pM
      F[vPos] = vM * VRR[vPos] + Normal[2] * pM
      F[wPos] = vM * VRR[wPos] + Normal[3] * pM
      F[ThPos] = vM * VRR[ThPos]
    end
  end
  return RiemannByLMARSNonLin!
end

function (::ArtianoEnergyStable)(Param,Phys,RhoPos,uPos,vPos,wPos,ThPos,pPos)
  @inline function RiemannByLMARSNonLin!(FL,FR,VLL,VRR,AuxL,AuxR,Normal)
FT = eltype(FL)
    cS = Param.cS
    pLL = AuxL[pPos]
    pRR = AuxR[pPos]

    uAv = (VLL[uPos]/VLL[RhoPos] + VRR[uPos]/VRR[RhoPos] ) * FT(0.5)		
    vAv = (VLL[vPos]/VLL[RhoPos] + VRR[vPos]/VRR[RhoPos] ) * FT(0.5)		
    wAv = (VLL[wPos]/VLL[RhoPos] + VRR[wPos]/VRR[RhoPos] ) * FT(0.5)		

    RhoM = FT(0.5) * (VLL[RhoPos] + VRR[RhoPos])
    ThM = FT(0.5) * (VLL[ThPos]/VLL[RhoPos] + VRR[ThPos]/VRR[RhoPos])
    vLL = (VLL[uPos] * Normal[1] + VLL[vPos] * Normal[2] + VLL[wPos] * Normal[3]) / VLL[RhoPos]
    vRR = (VRR[uPos] * Normal[1] + VRR[vPos] * Normal[2] + VRR[wPos] * Normal[3]) / VRR[RhoPos]
    norm_ = sqrt(Normal[1]^2 + Normal[2]^2 + Normal[3]^2) 
    vM = FT(0.5) * (vRR + vLL)
    pM = -FT(0.5) * Phys.Cpd * RhoM * ThM * (AuxR[3] - AuxL[3])
    diss = FT(0.5) * cS * RhoM * (vRR - vLL)/norm_ 
    if vM > FT(0)
    ThUpwind = VLL[ThPos]/VLL[RhoPos]
    else 
    ThUpwind = VRR[ThPos]/VRR[RhoPos]			
    end
	 Cad = FT(0.5) * abs((vRR + vLL))	

      FL[RhoPos] = RhoM * vM 
      FL[uPos] = FL[RhoPos] * uAv - diss * Normal[1] - FT(0.5) * RhoM * Cad * (VRR[uPos]/VRR[RhoPos] - VLL[uPos]/VLL[RhoPos])
      FL[vPos] = FL[RhoPos] * vAv - diss * Normal[2] - FT(0.5) * RhoM * Cad * (VRR[vPos]/VRR[RhoPos] - VLL[vPos]/VLL[RhoPos])
      FL[wPos] = FL[RhoPos] * wAv - diss * Normal[3] - FT(0.5) * RhoM * Cad * (VRR[wPos]/VRR[RhoPos] - VLL[wPos]/VLL[RhoPos])
      FL[ThPos] = FL[RhoPos] * ThUpwind
     
		FR[RhoPos] = FL[RhoPos]
		FR[uPos] = FL[uPos] - pM * Normal[1]
		FR[vPos] = FL[vPos] - pM * Normal[2]
		FR[wPos] = FL[wPos] - pM * Normal[3]
		FR[ThPos] = FL[ThPos]
		FL[RhoPos] = FL[RhoPos]
		FL[uPos] = FL[uPos] + pM * Normal[1]
		FL[vPos] = FL[vPos] + pM * Normal[2]
		FL[wPos] = FL[wPos] + pM * Normal[3]
		FL[ThPos] = FL[ThPos]

  end
  return RiemannByLMARSNonLin!
end


function (::RiemannExnerLMARS)(Param,Phys,RhoPos,uPos,vPos,wPos,ThPos,pPos)
  @inline function RiemannByLMARSNonLin!(FL,FR,VLL,VRR,AuxL,AuxR,Normal)
    FT = eltype(FL)
    cS = Param.cS
    pLL = AuxL[pPos]
    pRR = AuxR[pPos]
    RhoM = FT(0.5) * (VLL[RhoPos] + VRR[RhoPos])
    ThM = FT(0.5) * (VLL[ThPos]/VLL[RhoPos] + VRR[ThPos]/VRR[RhoPos])
    vLL = (VLL[uPos] * Normal[1] + VLL[vPos] * Normal[2] + VLL[wPos] * Normal[3]) / VLL[RhoPos]
    vRR = (VRR[uPos] * Normal[1] + VRR[vPos] * Normal[2] + VRR[wPos] * Normal[3]) / VRR[RhoPos]
    dp = -FT(0.5) * cS * RhoM * (vRR - vLL)
    vM = FT(0.5) * (vRR + vLL) - FT(1.0) /(FT(2.0) * cS) * (pRR - pLL) / RhoM
    pM = -FT(0.5) * Phys.Cpd * RhoM * ThM * (AuxR[3] - AuxL[3])
    #pM = FT(0.5) * (pLL + pRR)
    if vM > FT(0)
      FL[RhoPos] = vM * VLL[RhoPos]
      FL[uPos] = vM * VLL[uPos] + Normal[1] * dp
      FL[vPos] = vM * VLL[vPos] + Normal[2] * dp
      FL[wPos] = vM * VLL[wPos] + Normal[3] * dp
      FL[ThPos] = vM * VLL[ThPos]
    else
      FL[RhoPos] = vM * VRR[RhoPos]
      FL[uPos] = vM * VRR[uPos] + Normal[1] * dp
      FL[vPos] = vM * VRR[vPos] + Normal[2] * dp
      FL[wPos] = vM * VRR[wPos] + Normal[3] * dp
      FL[ThPos] = vM * VRR[ThPos]
    end

    FR[RhoPos] = FL[RhoPos] 
    FR[uPos] = FL[uPos] - pM * Normal[1]
    FR[vPos] = FL[vPos] - pM * Normal[2]
    FR[wPos] = FL[wPos] - pM * Normal[3]
    FR[ThPos] = FL[ThPos]
    
    #FL[RhoPos] = FL[RhoPos] 
    FL[uPos] = FL[uPos] + pM * Normal[1]
    FL[vPos] = FL[vPos] + pM * Normal[2]
    FL[wPos] = FL[wPos] + pM * Normal[3]
    #FL[ThPos] = FL[ThPos]	
  end
  return RiemannByLMARSNonLin!
end

function (::RiemannExLMARS)(Param,Phys,RhoPos,uPos,vPos,wPos,ThPos,pPos)
  @inline function RiemannByLMARSNonLin!(F,VLL,VRR,AuxL,AuxR,Normal)

    FT = eltype(F)
    cS = Param.cS
    pLL = AuxL[pPos]
    pRR = AuxR[pPos]
    ExpLL = AuxL[3]
    ExpRR = AuxR[3]
    RhoM = FT(0.5) * (VLL[RhoPos] + VRR[RhoPos])
    RhoThM = FT(0.5) * (VLL[ThPos] + VRR[ThPos])
    vLL = (VLL[uPos] * Normal[1] + VLL[vPos] * Normal[2] + VLL[wPos] * Normal[3]) / VLL[RhoPos]
    vRR = (VRR[uPos] * Normal[1] + VRR[vPos] * Normal[2] + VRR[wPos] * Normal[3]) / VRR[RhoPos]
#   pM = FT(0.5) * (pLL + pRR) - FT(0.5) * cS * RhoM * (vRR - vLL)
    vM = FT(0.5) * (vRR + vLL) - FT(1.0) /(FT(2.0) * cS) * (pRR - pLL) / RhoM
    pM = -FT(0.5) * cS * RhoM * (vRR - vLL)
    pM += 0.5 * Phys.Cpd * RhoThM * (ExpRR - ExpLL)  
    if vM > FT(0)
      F[RhoPos] = vM * VLL[RhoPos]
      F[uPos] = vM * VLL[uPos] + Normal[1] * pM
      F[vPos] = vM * VLL[vPos] + Normal[2] * pM
      F[wPos] = vM * VLL[wPos] + Normal[3] * pM
      F[ThPos] = vM * VLL[ThPos]
    else
      F[RhoPos] = vM * VRR[RhoPos]
      F[uPos] = vM * VRR[uPos] + Normal[1] * pM
      F[vPos] = vM * VRR[vPos] + Normal[2] * pM
      F[wPos] = vM * VRR[wPos] + Normal[3] * pM
      F[ThPos] = vM * VRR[ThPos]
    end
  end
  return RiemannByLMARSNonLin!
end

function (::RiemannExPLMARS)(Param,Phys,RhoPos,uPos,vPos,wPos,ThPos,pPos)
  @inline function RiemannByLMARSNonLin!(F,VLL,VRR,AuxL,AuxR,Normal)

    FT = eltype(F)
    cS = Param.cS
    pLL = AuxL[pPos]
    pRR = AuxR[pPos]
    ExpLL = AuxL[3]
    ExpRR = AuxR[3]
    RhoM = FT(0.5) * (VLL[RhoPos] + VRR[RhoPos])
    RhoThM = FT(0.5) * (VLL[ThPos] + VRR[ThPos])
    vLL = (VLL[uPos] * Normal[1] + VLL[vPos] * Normal[2] + VLL[wPos] * Normal[3]) / VLL[RhoPos]
    vRR = (VRR[uPos] * Normal[1] + VRR[vPos] * Normal[2] + VRR[wPos] * Normal[3]) / VRR[RhoPos]
    pM = FT(0.5) * Phys.Cpd * (VLL[ThPos] * ExpLL + VRR[ThPos]* ExpRR) - FT(0.5) * cS * RhoM * (vRR - vLL)
    vM = FT(0.5) * (vRR + vLL) - FT(1.0) /(FT(2.0) * cS) * (pRR - pLL) / RhoM
    if vM > FT(0)
      F[RhoPos] = vM * VLL[RhoPos]
      F[uPos] = vM * VLL[uPos] + Normal[1] * pM
      F[vPos] = vM * VLL[vPos] + Normal[2] * pM
      F[wPos] = vM * VLL[wPos] + Normal[3] * pM
      F[ThPos] = vM * VLL[ThPos]
    else
      F[RhoPos] = vM * VRR[RhoPos]
      F[uPos] = vM * VRR[uPos] + Normal[1] * pM
      F[vPos] = vM * VRR[vPos] + Normal[2] * pM
      F[wPos] = vM * VRR[wPos] + Normal[3] * pM
      F[ThPos] = vM * VRR[ThPos]
    end
  end
  return RiemannByLMARSNonLin!
end


function (::RiemannLMARSFast)(Param,Phys,RhoPos,uPos,vPos,wPos,ThPos,pPos)
  @inline function RiemannByLMARSNonLin!(F,VLL,VRR,AuxL,AuxR,Normal)

    FT = eltype(F)
    cS = Param.cS
    pLL = AuxL[dpdRhoThPos] * VLL[RhoThPos]
    pRR = AuxR[dpdRhoThPos] * VRR[RhoThPos]
    RhoM = FT(0.5) * (VLL[RhoPos] + VRR[RhoPos])
    vLL = (VLL[uPos] * Normal[1] + VLL[vPos] * Normal[2] + VLL[wPos] * Normal[3]) / VLL[RhoPos]
    vRR = (VRR[uPos] * Normal[1] + VRR[vPos] * Normal[2] + VRR[wPos] * Normal[3]) / VRR[RhoPos]
    pM = FT(0.5) * (pLL + pRR) - FT(0.5) * cS * RhoM * (vRR - vLL)
    vM = FT(0.5) * (vRR + vLL) - FT(1.0) /(FT(2.0) * cS) * (pRR - pLL) / RhoM
    F[uPos] = Normal[1] * pM
    F[vPos] = Normal[2] * pM
    F[wPos] = Normal[3] * pM
    if vM > FT(0)
      F[RhoPos] = vM * VLL[RhoPos]
      F[RhoThPos] = F[RhoPos] * AuxL[ThPos]
    else
      F[RhoPos] = vM * VRR[RhoPos]
      F[RhoThPos] = F[RhoPos] * AuxR[ThPos]
    end
  end
  return RiemannByLMARSLin!
end
Base.@kwdef struct RiemannBoussinesqLMARS <: RiemannSolver end

function (::RiemannBoussinesqLMARS)(Param,pPos,uPos,vPos,wPos,bPos)
  @inline function RiemannByLMARS(F,VLL,VRR,AuxL,AuxR,Normal)

    FT = eltype(F)
    cS = Param.cS
    pLL = VLL[pPos]
    pRR = VRR[pPos]
    U = Normal[1] * Param.U
    vLL = VLL[uPos] * Normal[1] + VLL[vPos] * Normal[2] + VLL[wPos] * Normal[3] 
    vRR = VRR[uPos] * Normal[1] + VRR[vPos] * Normal[2] + VRR[wPos] * Normal[3]
    pM = FT(0.5) * (pLL + pRR) - FT(0.5) * cS * (vRR - vLL)
    vM = FT(0.5) * (vRR + vLL) - FT(1.0) /(FT(2.0) * cS) * (pRR - pLL) 
    if Param.U > 0
      F[pPos] = U * VLL[pPos] + Param.cS^2*vM 
      F[uPos] = U * VLL[uPos] + Normal[1] * pM
      F[vPos] = U * VLL[vPos] + Normal[2] * pM
      F[wPos] = U * VLL[wPos] + Normal[3] * pM
      F[bPos] = U * VLL[bPos]
    else
      F[pPos] = U * VRR[pPos] + Param.cS^2*vM 
      F[uPos] = U * VRR[uPos] + Normal[1] * pM 
      F[vPos] = U * VRR[vPos] + Normal[2] * pM
      F[wPos] = U * VRR[wPos] + Normal[3] * pM
      F[bPos] = U * VRR[bPos]
    end
  end
  return RiemannByLMARS
end

@inline function ln_mean(x, y)
  FT = eltype(x)  
  if abs(x-y) <= max(x,y) * eps(FT)
    return FT(0.5) * (x + y)
  else
    return (x - y) / log(x / y)  
  end
end  

@inline function stolarsky_mean(x, y, gamma) 
  FT = eltype(x)
  epsilon_f2 = convert(FT, 1.0e-4)
  f2 = (x * (x - 2 * y) + y * y) / (x * (x + 2 * y) + y * y) # f2 = f^2
   if f2 < epsilon_f2
     # convenience coefficients
     c1 = FT(1 / 3) * (gamma - FT(2))
     c2 = FT(-1 / 15) * (gamma + FT(1)) * (gamma - FT(3)) * c1
     c3 = FT(-1 / 21) * (2 * gamma * (gamma - FT(2)) - FT(9)) * c2
    return FT(0.5) * (x + y) * @evalpoly(f2, FT(1), c1, c2, c3)
  else
    return (gamma - FT(1)) / gamma * (y^gamma - x^gamma) /
     (y^(gamma - FT(1)) - x^(gamma - FT(1)))
  end
end






