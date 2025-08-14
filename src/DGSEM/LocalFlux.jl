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

Base.@kwdef struct KennedyGruberGrav <: AverageFlux end

function (::KennedyGruberGrav)(RhoPos,uPos,vPos,wPos,ThPos,pPos,GPPos)
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
    flux[2] = flux[1] * uAv + mAv1 * pAv
    flux[3] = flux[1] * vAv + mAv2 * pAv
    flux[4] = flux[1] * wAv + mAv3 * pAv
    flux[5] = flux[1] * ThAv

  end
  return FluxNonLinAver!
end

Base.@kwdef struct KennedyGruberGravMod <: AverageFlux end

function (::KennedyGruberGravMod)(RhoPos,uPos,vPos,wPos,ThPos,pPos,GPPos)
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
    
    #pAv = FT(0.5) * ((pL + pR) + FT(0.5) * (RhoL + RhoR) * (GPR - GPL))
    uAv = FT(0.5) * (uL + uR)
    vAv = FT(0.5) * (vL + vR)
    wAv = FT(0.5) * (wL + wR)
    RhoAv = FT(0.5) * (RhoL + RhoR)
    ThAv = FT(0.5) * (ThL + ThR)

    mAv1 = FT(0.5) * (m_L[1] + m_R[1])
    mAv2 = FT(0.5) * (m_L[2] + m_R[2])
    mAv3 = FT(0.5) * (m_L[3] + m_R[3])
    qHat = mAv1 * uAv + mAv2 * vAv + mAv3 * wAv
    gamma = FT(1004.0)/FT(717.0)
    p0 = FT(100000.0)
    R = FT(287.0)
    dpdrhoL = gamma * pL / (RhoL * ThL)
    dpdrhoR = gamma * pR / (RhoR * ThR)
    dpdrhoAv = FT(0.5) * (dpdrhoL + dpdrhoR)
    grav_term = FT(0.5) * RhoAv * (GPR - GPL)
    dpdrho = FT(0.5) * dpdrhoL * (RhoR * ThR - RhoL * ThL)
    # @show dpdrhoL * RhoR * ThR
    # @show pL
    # @show dpdrhoL
    # @show gamma * R * (R * RhoL * ThL/p0)^(gamma-1.0) * RhoR * ThR
    pAv = grav_term + dpdrho
    flux[1] = RhoAv * qHat
    flux[2] = flux[1] * uAv + mAv1 * pAv
    flux[3] = flux[1] * vAv + mAv2 * pAv
    flux[4] = flux[1] * wAv + mAv3 * pAv
    flux[5] = flux[1] * ThAv
  end
  return FluxNonLinAver!
end

abstract type RiemannSolver end

Base.@kwdef struct RiemannLMARS <: RiemannSolver end
Base.@kwdef struct RiemannLMARSMod <: RiemannSolver end

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

function (::RiemannLMARSMod)(Param,Phys,RhoPos,uPos,vPos,wPos,ThPos,pPos)
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

    RhoL = VLL[RhoPos]
    RhoR = VRR[RhoPos]
    ThL = VLL[ThPos] / RhoL
    ThR = VRR[ThPos] / RhoR
    dpdrhoL = 1.4002789f0 * pLL / (VLL[RhoPos] * VLL[ThPos])
    dpdrhoR = 1.4002789f0 * pRR / (VRR[RhoPos] * VRR[ThPos])
    if FT(0.5) * (vLL + vRR) > FT(0)
    pM = FT(0.5) * dpdrhoL * (RhoR * ThR - RhoL * ThL)
    else
    pM = FT(0.5) * dpdrhoR * (RhoR * ThR - RhoL * ThL)
    end
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






