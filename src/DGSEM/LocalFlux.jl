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

    flux[1,RhoPos] = -V[uPos]
    flux[1,uPos] = -V[uPos] * u - p
    flux[1,vPos] = -V[uPos] * v
    flux[1,wPos] = -V[uPos] * w
    flux[1,ThPos] = -V[uPos] * Th

    flux[2,RhoPos] = -V[vPos]
    flux[2,uPos] = -V[vPos] * u
    flux[2,vPos] = -V[vPos] * v - p
    flux[2,wPos] = -V[vPos] * w
    flux[2,ThPos] = -V[vPos] * Th


    flux[3,RhoPos] = -V[wPos]
    flux[3,uPos] = -V[wPos] * u
    flux[3,vPos] = -V[wPos] * v
    flux[3,wPos] = -V[wPos] * w - p
    flux[3,ThPos] = -V[wPos] * Th
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


abstract type RiemannSolver end

Base.@kwdef struct RiemannLMARS <: RiemannSolver end

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







#=
function flux_nonconservative_waruszewski(u_ll, u_rr, normal_direction::AbstractVector,
                                          equations::CompressibleEulerEquationsWithGravity2D)
    rho_ll, _, _, _, phi_ll = u_ll
    rho_rr, _, _, _, phi_rr = u_rr

    # We omit the 0.5 in the density average since Trixi.jl always multiplies the non-conservative flux with 0.5
    noncons = ln_mean(rho_ll, rho_rr) * (phi_rr - phi_ll)
    # noncons = 0.5 * (rho_ll + rho_rr) * (phi_rr - phi_ll)

    f0 = zero(eltype(u_ll))
    return SVector(f0, noncons * normal_direction[1], noncons * normal_direction[2],
                   f0, f0)
end
=#
