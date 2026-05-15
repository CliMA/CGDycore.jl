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


function (::KennedyGruberGrav)(RhoPos, uPos, vPos, wPos, ThPos, pPos, GPPos, ::Grids.Quad)

  @inline function FluxNonLinAver!(flux,
      VLoc, AuxLoc, dXdxILoc,
      K1, Iz1, iD1,   # left state indices  (localidx into @localmem)
      K2, Iz2, iD2,   # right state indices
      ::Val{dir}) where {dir}

    FT = eltype(flux)

    # ------ read left state directly from shared memory ------
    RhoL  = VLoc[K1, Iz1, iD1, RhoPos]
    uL    = VLoc[K1, Iz1, iD1, uPos]
    vL    = VLoc[K1, Iz1, iD1, vPos]
    wL    = VLoc[K1, Iz1, iD1, wPos]
    ThL   = VLoc[K1, Iz1, iD1, ThPos]
    pL    = AuxLoc[K1, Iz1, iD1, pPos]
    GPL   = AuxLoc[K1, Iz1, iD1, GPPos]

    # ------ read right state directly from shared memory -----
    RhoR  = VLoc[K2, Iz2, iD2, RhoPos]
    uR    = VLoc[K2, Iz2, iD2, uPos]
    vR    = VLoc[K2, Iz2, iD2, vPos]
    wR    = VLoc[K2, Iz2, iD2, wPos]
    ThR   = VLoc[K2, Iz2, iD2, ThPos]
    pR    = AuxLoc[K2, Iz2, iD2, pPos]
    GPR   = AuxLoc[K2, Iz2, iD2, GPPos]

    # ------ read metric (dXdxI row) from shared memory -------
    # dXdxILoc layout: (dir, j, K, Iz, iD)  — pass dir as Val for unrolling
    m_L1  = dXdxILoc[dir, 1, K1, Iz1, iD1]
    m_L2  = dXdxILoc[dir, 2, K1, Iz1, iD1]
    m_L3  = dXdxILoc[dir, 3, K1, Iz1, iD1]
    m_R1  = dXdxILoc[dir, 1, K2, Iz2, iD2]
    m_R2  = dXdxILoc[dir, 2, K2, Iz2, iD2]
    m_R3  = dXdxILoc[dir, 3, K2, Iz2, iD2]

    # ------ Kennedy-Gruber averages --------------------------
    RhoAv = FT(0.5) * (RhoL + RhoR)
    pAv   = FT(0.5) * ((pL + pR) + RhoAv * (GPR - GPL))
    uAv   = FT(0.5) * (uL / RhoL + uR / RhoR)
    vAv   = FT(0.5) * (vL / RhoL + vR / RhoR)
    wAv   = FT(0.5) * (wL / RhoL + wR / RhoR)
    ThAv  = FT(0.5) * (ThL / RhoL + ThR / RhoR)
    mAv1  = FT(0.5) * (m_L1 + m_R1)
    mAv2  = FT(0.5) * (m_L2 + m_R2)
    mAv3  = FT(0.5) * (m_L3 + m_R3)

    qHat  = mAv1 * uAv + mAv2 * vAv + mAv3 * wAv

    flux[1] = RhoAv * qHat
    flux[2] = flux[1] * uAv + mAv1 * pAv
    flux[3] = flux[1] * vAv + mAv2 * pAv
    flux[4] = flux[1] * wAv + mAv3 * pAv
    flux[5] = flux[1] * ThAv
  end

  return FluxNonLinAver!
end

 
function (::KennedyGruberGrav)(RhoPos, uPos, vPos, wPos, ThPos, pPos, GPPos, ::Grids.Tri)

  @inline function FluxNonLinAver!(flux,
      VLoc, AuxLoc, dXdxILoc,
      ID1, iz1,
      ID2, iz2,
      ::Val{dir}) where {dir}   # dir = 1 for fTilde, 2 for gTilde

    FT = eltype(flux)

    RhoL = VLoc[ID1, iz1, RhoPos]
    uL   = VLoc[ID1, iz1, uPos]
    vL   = VLoc[ID1, iz1, vPos]
    wL   = VLoc[ID1, iz1, wPos]
    ThL  = VLoc[ID1, iz1, ThPos]
    pL   = AuxLoc[ID1, iz1, pPos]
    GPL  = AuxLoc[ID1, iz1, GPPos]

    RhoR = VLoc[ID2, iz2, RhoPos]
    uR   = VLoc[ID2, iz2, uPos]
    vR   = VLoc[ID2, iz2, vPos]
    wR   = VLoc[ID2, iz2, wPos]
    ThR  = VLoc[ID2, iz2, ThPos]
    pR   = AuxLoc[ID2, iz2, pPos]
    GPR  = AuxLoc[ID2, iz2, GPPos]

    m_L1 = dXdxILoc[dir, 1, ID1, iz1]
    m_L2 = dXdxILoc[dir, 2, ID1, iz1]
    m_L3 = dXdxILoc[dir, 3, ID1, iz1]
    m_R1 = dXdxILoc[dir, 1, ID2, iz2]
    m_R2 = dXdxILoc[dir, 2, ID2, iz2]
    m_R3 = dXdxILoc[dir, 3, ID2, iz2]

    RhoAv = FT(0.5) * (RhoL + RhoR)
    pAv   = FT(0.5) * ((pL + pR) + RhoAv * (GPR - GPL))
    uAv   = FT(0.5) * (uL / RhoL + uR / RhoR)
    vAv   = FT(0.5) * (vL / RhoL + vR / RhoR)
    wAv   = FT(0.5) * (wL / RhoL + wR / RhoR)
    ThAv  = FT(0.5) * (ThL / RhoL + ThR / RhoR)
    mAv1  = FT(0.5) * (m_L1 + m_R1)
    mAv2  = FT(0.5) * (m_L2 + m_R2)
    mAv3  = FT(0.5) * (m_L3 + m_R3)
    qHat  = mAv1 * uAv + mAv2 * vAv + mAv3 * wAv

    flux[1] = RhoAv * qHat
    flux[2] = flux[1] * uAv + mAv1 * pAv
    flux[3] = flux[1] * vAv + mAv2 * pAv
    flux[4] = flux[1] * wAv + mAv3 * pAv
    flux[5] = flux[1] * ThAv
  end

  return FluxNonLinAver!
end

Base.@kwdef struct KennedyGruberIEGrav <: AverageFlux end

function (::KennedyGruberIEGrav)(RhoPos,uPos,vPos,wPos,IEPos,pPos,GPPos)
  @inline function FluxNonLinAver!(flux,VL,VR,AuxL,AuxR,m_L,m_R)
    FT = eltype(flux)
    pAv = FT(0.5) * ((AuxL[pPos] + AuxR[pPos]) +
			FT(0.5) * (VL[RhoPos] + VR[RhoPos]) * (AuxR[GPPos] - AuxL[GPPos]))
    uAv = FT(0.5) * (VL[uPos] + VR[uPos])
    vAv = FT(0.5) * (VL[vPos] + VR[vPos])
    wAv = FT(0.5) * (VL[wPos] + VR[wPos])
    RhoAv = FT(0.5) * (VL[RhoPos] + VR[RhoPos])
    IEAv = FT(0.5) * (VL[IEPos] + VR[IEPos])
    mAv1 = FT(0.5) * (m_L[1] + m_R[1])
    mAv2 = FT(0.5) * (m_L[2] + m_R[2])
    mAv3 = FT(0.5) * (m_L[3] + m_R[3])
    qRHat = VR[uPos] * mAv1 + VR[vPos] * mAv2 + VR[wPos] * mAv3
    qLHat = VL[uPos] * mAv1 + VL[vPos] * mAv2 + VL[wPos] * mAv3
    qHat = FT(0.5) * (qLHat + qRHat)
    flux[1] = RhoAv * qHat
    flux[2] = flux[1] * uAv + mAv1 * pAv
    flux[3] = flux[1] * vAv + mAv2 * pAv
    flux[4] = flux[1] * wAv + mAv3 * pAv
    flux[5] = flux[1] * IEAv + 
      FT(0.25) * (AuxL[pPos] - AuxR[pPos]) * (qLHat-qRHat) + AuxL[pPos] * qHat

  end
  return FluxNonLinAver!
end

Base.@kwdef struct KennedyGruberGravFast <: AverageFlux end

function (::KennedyGruberGravFast)(RhoPos,uPos,vPos,wPos,ThPos,pPos,GPPos)
  @inline function FluxNonLinAver!(flux,VL,VR,AuxL,AuxR,m_L,m_R)
    FT = eltype(flux)
    RhoAv = FT(0.5) * (VL[RhoPos] + VR[RhoPos])
    pL = AuxL[pPos]
    pR = AuxR[pPos]
#   pAv = FT(0.5) * ((AuxL[pPos] + AuxR[pPos]) + RhoAv * (AuxR[GPPos] - AuxL[GPPos]))
    pAv = FT(0.5) * ((pL + pR) + RhoAv * (AuxR[GPPos] - AuxL[GPPos]))
    uAv = FT(0.5) * (VL[uPos] / VL[RhoPos] + VR[uPos] / VR[RhoPos])
    vAv = FT(0.5) * (VL[vPos] / VL[RhoPos] + VR[vPos] / VR[RhoPos])
    wAv = FT(0.5) * (VL[wPos] / VL[RhoPos] + VR[wPos] / VR[RhoPos])
    ThAv = FT(0.5) * (VL[ThPos] / VL[RhoPos] + VR[ThPos] / VR[RhoPos])
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

function (::KennedyGruberGravFast)(RhoPos,uPos,vPos,wPos,ThPos,pRhoThPos,ThAuxPos,GPPos)
  @inline function FluxNonLinAverSemi!(flux,VL,VR,AuxL,AuxR,m_L,m_R)
    FT = eltype(flux)
    RhoAv = FT(0.5) * (VL[RhoPos] + VR[RhoPos])
    pL = AuxL[pRhoThPos] * VL[ThPos]
    pR = AuxR[pRhoThPos] * VR[ThPos]
    pAv = FT(0.5) * ((pL + pR) + RhoAv * (AuxR[GPPos] - AuxL[GPPos]))
    uAv = FT(0.5) * (VL[uPos] / VL[RhoPos] + VR[uPos] / VR[RhoPos])
    vAv = FT(0.5) * (VL[vPos] / VL[RhoPos] + VR[vPos] / VR[RhoPos])
    wAv = FT(0.5) * (VL[wPos] / VL[RhoPos] + VR[wPos] / VR[RhoPos])
    ThAv = FT(0.5) * (AuxL[ThAuxPos] + AuxR[ThAuxPos])
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
  return FluxNonLinAverSemi!
end

Base.@kwdef struct KennedyGruberGravLinFast <: AverageFlux end

function (::KennedyGruberGravLinFast)(RhoPos,uPos,vPos,wPos,ThPos,pRhoThPos,ThAuxPos,GPPos,::Grids.Tri)
  @inline function FluxNonLinAverSemi!(flux,
      VLoc, AuxLoc, dXdxILoc,
      ID1, iz1,
      ID2, iz2,
      ::Val{dir}) where {dir}   # dir = 1 for fTilde, 2 for gTilde
    FT = eltype(flux)
    RhoL = VLoc[ID1, iz1, RhoPos]
    uL   = VLoc[ID1, iz1, uPos]
    vL   = VLoc[ID1, iz1, vPos]
    wL   = VLoc[ID1, iz1, wPos]
    ThL  = VLoc[ID1, iz1, ThPos]
    pLTh   = AuxLoc[ID1, iz1, pRhoThpos]
    GPL  = AuxLoc[ID1, iz1, GPPos]

    RhoR = VLoc[ID2, iz2, RhoPos]
    uR   = VLoc[ID2, iz2, uPos]
    vR   = VLoc[ID2, iz2, vPos]
    wR   = VLoc[ID2, iz2, wPos]
    ThR  = VLoc[ID2, iz2, ThPos]
    pRTh   = AuxLoc[ID2, iz2, pRhoThPos]
    GPR  = AuxLoc[ID2, iz2, GPPos]

    m_L1 = dXdxILoc[dir, 1, ID1, iz1]
    m_L2 = dXdxILoc[dir, 2, ID1, iz1]
    m_L3 = dXdxILoc[dir, 3, ID1, iz1]
    m_R1 = dXdxILoc[dir, 1, ID2, iz2]
    m_R2 = dXdxILoc[dir, 2, ID2, iz2]
    m_R3 = dXdxILoc[dir, 3, ID2, iz2]
    RhoAv = FT(0.5) * (RhoL + RhoR)
    pL = pLTh* ThL
    pR = pRTh* ThR
    pAv = FT(0.5) * ((pL + pR) + RhoAv * (GPR - GPL))
    uAv = FT(0.5) * (uL + uR)
    vAv = FT(0.5) * (vL + vR)
    wAv = FT(0.5) * (wL + wR)
    ThAv = FT(0.5) * (AuxL[ID1, iz1, ThAuxPos] + AuxR[ID2, iz2, ThAuxPos])
    mAv1 = FT(0.5) * (m_L1 + m_R1)
    mAv2 = FT(0.5) * (m_L2 + m_R2)
    mAv3 = FT(0.5) * (m_L3 + m_R3)
    qHat = mAv1 * uAv + mAv2 * vAv + mAv3 * wAv
    flux[1] = qHat
    flux[2] = mAv1 * pAv
    flux[3] = mAv2 * pAv
    flux[4] = mAv3 * pAv
    flux[5] = flux[1] * ThAv

  end
  return FluxNonLinAverSemi!
end


function (::KennedyGruberGravLinFast)(RhoPos,uPos,vPos,wPos,ThPos,pRhoThPos,ThAuxPos,GPPos,::Grids.Quad)
  @inline function FluxNonLinAverSemi!(flux,
      VLoc, AuxLoc, dXdxILoc,
      K1, Iz1, iD1,   # left state indices  (localidx into @localmem)
      K2, Iz2, iD2,   # right state indices
      ::Val{dir}) where {dir}
    FT = eltype(flux)
    RhoL  = VLoc[K1, Iz1, iD1, RhoPos]
    uL    = VLoc[K1, Iz1, iD1, uPos]
    vL    = VLoc[K1, Iz1, iD1, vPos]
    wL    = VLoc[K1, Iz1, iD1, wPos]
    ThL   = VLoc[K1, Iz1, iD1, ThPos]
    pLTh   = AuxLoc[K1, Iz1, iD1, pRhoThPos]
    GPL  = AuxLoc[K1, Iz1, iD1, GPPos]

    RhoR  = VLoc[K2, Iz2, iD2, RhoPos]
    uR    = VLoc[K2, Iz2, iD2, uPos]
    vR    = VLoc[K2, Iz2, iD2, vPos]
    wR    = VLoc[K2, Iz2, iD2, wPos]
    ThR   = VLoc[K2, Iz2, iD2, ThPos]
    pRTh   = AuxLoc[K2, Iz2, iD2, pRhoThPos]
    GPR  = AuxLoc[K2, Iz2, iD2, GPPos]

    m_L1  = dXdxILoc[dir, 1, K1, Iz1, iD1]
    m_L2  = dXdxILoc[dir, 2, K1, Iz1, iD1]
    m_L3  = dXdxILoc[dir, 3, K1, Iz1, iD1]
    m_R1  = dXdxILoc[dir, 1, K2, Iz2, iD2]
    m_R2  = dXdxILoc[dir, 2, K2, Iz2, iD2]
    m_R3  = dXdxILoc[dir, 3, K2, Iz2, iD2]
    RhoAv = FT(0.5) * (RhoL + RhoR)
    pL = pLTh* ThL
    pR = pRTh* ThR
    pAv = FT(0.5) * ((pL + pR) + RhoAv * (GPR - GPL))
    uAv = FT(0.5) * (uL + uR)
    vAv = FT(0.5) * (vL + vR)
    wAv = FT(0.5) * (wL + wR)
    ThAv = FT(0.5) * (AuxLoc[K1, Iz1, iD1, ThAuxPos] + AuxLoc[K2, Iz2, iD2, ThAuxPos])
    mAv1 = FT(0.5) * (m_L1 + m_R1)
    mAv2 = FT(0.5) * (m_L2 + m_R2)
    mAv3 = FT(0.5) * (m_L3 + m_R3)
    qHat = mAv1 * uAv + mAv2 * vAv + mAv3 * wAv
    flux[1] = qHat
    flux[2] = mAv1 * pAv
    flux[3] = mAv2 * pAv
    flux[4] = mAv3 * pAv
    flux[5] = flux[1] * ThAv

  end
  return FluxNonLinAverSemi!
end

Base.@kwdef struct KennedyGruberGravLinFast1 <: AverageFlux end

function (::KennedyGruberGravLinFast1)(RhoPos,uPos,vPos,wPos,ThPos,pRhoThPos,ThAuxPos,GPPos,::Grids.Quad)
  @inline function FluxNonLinAverSemi!(flux,
      VLoc, AuxLoc, dXdxILoc,
      K1, Iz1, iD1,   # left state indices  (localidx into @localmem)
      K2, Iz2, iD2,   # right state indices
      ::Val{dir}) where {dir}
    FT = eltype(flux)
    RhoL  = VLoc[K1, Iz1, iD1, RhoPos]
    uL    = VLoc[K1, Iz1, iD1, uPos]
    vL    = VLoc[K1, Iz1, iD1, vPos]
    wL    = VLoc[K1, Iz1, iD1, wPos]
    RhoThL   = VLoc[K1, Iz1, iD1, ThPos]
    pLTh   = AuxLoc[K1, Iz1, iD1, pRhoThPos]
    GPL  = AuxLoc[K1, Iz1, iD1, GPPos]

    RhoR  = VLoc[K2, Iz2, iD2, RhoPos]
    uR    = VLoc[K2, Iz2, iD2, uPos]
    vR    = VLoc[K2, Iz2, iD2, vPos]
    wR    = VLoc[K2, Iz2, iD2, wPos]
    RhoThR   = VLoc[K2, Iz2, iD2, ThPos]
    pRTh   = AuxLoc[K2, Iz2, iD2, pRhoThPos]
    GPR  = AuxLoc[K2, Iz2, iD2, GPPos]

    m_L1  = dXdxILoc[dir, 1, K1, Iz1, iD1]
    m_L2  = dXdxILoc[dir, 2, K1, Iz1, iD1]
    m_L3  = dXdxILoc[dir, 3, K1, Iz1, iD1]
    m_R1  = dXdxILoc[dir, 1, K2, Iz2, iD2]
    m_R2  = dXdxILoc[dir, 2, K2, Iz2, iD2]
    m_R3  = dXdxILoc[dir, 3, K2, Iz2, iD2]
    pL = pLTh * RhoThL
    pR = pRTh * RhoThR
    pAvL = pL  + FT(0.5) * RhoL * (GPR - GPL)
    pAvR = pR  + FT(0.5) * RhoR * (GPR - GPL)
    ThL = AuxLoc[K1, Iz1, iD1, ThAuxPos]
    ThR = AuxLoc[K2, Iz2, iD2, ThAuxPos]
    qHatL = m_L1 * uL + m_L2 * vL + m_L3 * wL
    qHatR = m_R1 * uR + m_R2 * vR + m_R3 * wR
    flux[1] = FT(0.5) * (qHatL + qHatR)
    flux[2] = FT(0.5) * (m_L1 * pAvL + m_R1 * pAvR)
    flux[3] = FT(0.5) * (m_L2 * pAvL + m_R2 * pAvR)
    flux[4] = FT(0.5) * (m_L3 * pAvL + m_R3 * pAvR)
    flux[5] = FT(0.5) * (qHatL * ThL + qHatR * ThR)

  end
  return FluxNonLinAverSemi!
end

Base.@kwdef struct KennedyGruberGravSlow <: AverageFlux end

function (::KennedyGruberGravSlow)(RhoPos,uPos,vPos,wPos,ThPos,pPos,GPPos)
  @inline function FluxNonLinAver!(flux,VL,VR,AuxL,AuxR,m_L,m_R)
    FT = eltype(flux)
    uAv = FT(0.5) * (VL[uPos] / VL[RhoPos] + VR[uPos] / VR[RhoPos])
    vAv = FT(0.5) * (VL[vPos] / VL[RhoPos] + VR[vPos] / VR[RhoPos])
    wAv = FT(0.5) * (VL[wPos] / VL[RhoPos] + VR[wPos] / VR[RhoPos])
    RhoAv = FT(0.5) * (VL[RhoPos] + VR[RhoPos])
    mAv1 = FT(0.5) * (m_L[1] + m_R[1])
    mAv2 = FT(0.5) * (m_L[2] + m_R[2])
    mAv3 = FT(0.5) * (m_L[3] + m_R[3])
    qHat = mAv1 * uAv + mAv2 * vAv + mAv3 * wAv
    flux[1] = FT(0)
    flux[2] = RhoAv * qHat * uAv
    flux[3] = RhoAv * qHat * vAv
    flux[4] = RhoAv * qHat * wAv
    flux[5] = FT(0)
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
Base.@kwdef struct RiemannLMARSLinFast <: RiemannSolver end
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
  @inline function RiemannByLMARSNonLin!(F,VLL,VRR,AuxL,AuxR,n1,n2,n3)
    FT = eltype(F)
    cS = Param.cS
    pLL = AuxL[pPos]
    pRR = AuxR[pPos]
    RhoM = FT(0.5) * (VLL[RhoPos] + VRR[RhoPos])
    vLL = (VLL[uPos] * n1 + VLL[vPos] * n2 + VLL[wPos] * n3) / VLL[RhoPos]
    vRR = (VRR[uPos] * n1 + VRR[vPos] * n2 + VRR[wPos] * n3) / VRR[RhoPos]
    pM = FT(0.5) * (pLL + pRR) - FT(0.5) * cS * RhoM * (vRR - vLL)
    vM = FT(0.5) * (vRR + vLL) - FT(1.0) /(FT(2.0) * cS) * (pRR - pLL) / RhoM
    if vM > FT(0)
      F[RhoPos] = vM * VLL[RhoPos]
      F[uPos] = vM * VLL[uPos] + n1 * pM
      F[vPos] = vM * VLL[vPos] + n2 * pM
      F[wPos] = vM * VLL[wPos] + n3 * pM
      F[ThPos] = vM * VLL[ThPos]
    else
      F[RhoPos] = vM * VRR[RhoPos]
      F[uPos] = vM * VRR[uPos] + n1 * pM
      F[vPos] = vM * VRR[vPos] + n2 * pM
      F[wPos] = vM * VRR[wPos] + n3 * pM
      F[ThPos] = vM * VRR[ThPos]
    end
  end
  return RiemannByLMARSNonLin!
end

function (::RiemannLMARSSlow)(Param,Phys,RhoPos,uPos,vPos,wPos,ThPos,pPos)
  @inline function RiemannByLMARSNonLin!(F,VLL,VRR,AuxL,AuxR,n1,n2,n3)
    FT = eltype(F)
    cS = Param.cS
    pLL = AuxL[pPos]
    pRR = AuxR[pPos]
    RhoM = FT(0.5) * (VLL[RhoPos] + VRR[RhoPos])
    vLL = (VLL[uPos] * n1 + VLL[vPos] * n2 + VLL[wPos] * n3) / VLL[RhoPos]
    vRR = (VRR[uPos] * n1 + VRR[vPos] * n2 + VRR[wPos] * n3) / VRR[RhoPos]
    vM = FT(0.5) * (vRR + vLL) - FT(1.0) /(FT(2.0) * cS) * (pRR - pLL) / RhoM
    if vM > FT(0)
      F[RhoPos] = FT(0)
      F[uPos] = vM * VLL[uPos]
      F[vPos] = vM * VLL[vPos]
      F[wPos] = vM * VLL[wPos]
      F[ThPos] = FT(0)
    else
      F[RhoPos] = FT(0)
      F[uPos] = vM * VRR[uPos]
      F[vPos] = vM * VRR[vPos]
      F[wPos] = vM * VRR[wPos]
      F[ThPos] = FT(0)
    end
  end
  return RiemannByLMARSNonLin!
end

function (::ArtianoEnergyStable)(Param,Phys,RhoPos,uPos,vPos,wPos,ThPos,pPos)
  @inline function RiemannByLMARSNonLin!(FL,FR,VLL,VRR,AuxL,AuxR,n1,n2,n3)
    FT = eltype(FL)
    cS = Param.cS
    pLL = AuxL[pPos]
    pRR = AuxR[pPos]

    uAv = (VLL[uPos]/VLL[RhoPos] + VRR[uPos]/VRR[RhoPos] ) * FT(0.5)		
    vAv = (VLL[vPos]/VLL[RhoPos] + VRR[vPos]/VRR[RhoPos] ) * FT(0.5)		
    wAv = (VLL[wPos]/VLL[RhoPos] + VRR[wPos]/VRR[RhoPos] ) * FT(0.5)		

    RhoM = FT(0.5) * (VLL[RhoPos] + VRR[RhoPos])
    ThM = FT(0.5) * (VLL[ThPos]/VLL[RhoPos] + VRR[ThPos]/VRR[RhoPos])
    vLL = (VLL[uPos] * n1 + VLL[vPos] * n2 + VLL[wPos] * n3) / VLL[RhoPos]
    vRR = (VRR[uPos] * n1 + VRR[vPos] * n2 + VRR[wPos] * n3) / VRR[RhoPos]
    norm_ = sqrt(n1^2 + n2^2 + n3^2) 
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
    FL[uPos] = FL[RhoPos] * uAv - diss * n1 - FT(0.5) * RhoM * Cad * (VRR[uPos]/VRR[RhoPos] - VLL[uPos]/VLL[RhoPos])
    FL[vPos] = FL[RhoPos] * vAv - diss * n2 - FT(0.5) * RhoM * Cad * (VRR[vPos]/VRR[RhoPos] - VLL[vPos]/VLL[RhoPos])
    FL[wPos] = FL[RhoPos] * wAv - diss * n3 - FT(0.5) * RhoM * Cad * (VRR[wPos]/VRR[RhoPos] - VLL[wPos]/VLL[RhoPos])
    FL[ThPos] = FL[RhoPos] * ThUpwind
       
    FR[RhoPos] = FL[RhoPos]
    FR[uPos] = FL[uPos] - pM * n1
    FR[vPos] = FL[vPos] - pM * n2
    FR[wPos] = FL[wPos] - pM * n3
    FR[ThPos] = FL[ThPos]
    FL[RhoPos] = FL[RhoPos]
    FL[uPos] = FL[uPos] + pM * n1
    FL[vPos] = FL[vPos] + pM * n2
    FL[wPos] = FL[wPos] + pM * n3
    FL[ThPos] = FL[ThPos]

  end
  return RiemannByLMARSNonLin!
end


function (::RiemannExnerLMARS)(Param,Phys,RhoPos,uPos,vPos,wPos,ThPos,pPos)
  @inline function RiemannByLMARSNonLin!(FL,FR,VLL,VRR,AuxL,AuxR,n1,n2,n3)
    FT = eltype(FL)
    cS = Param.cS
    pLL = AuxL[pPos]
    pRR = AuxR[pPos]
    RhoM = FT(0.5) * (VLL[RhoPos] + VRR[RhoPos])
    ThM = FT(0.5) * (VLL[ThPos]/VLL[RhoPos] + VRR[ThPos]/VRR[RhoPos])
    vLL = (VLL[uPos] * n1 + VLL[vPos] * n2 + VLL[wPos] * n3) / VLL[RhoPos]
    vRR = (VRR[uPos] * n1 + VRR[vPos] * n2 + VRR[wPos] * n3) / VRR[RhoPos]
    dp = -FT(0.5) * cS * RhoM * (vRR - vLL)
    vM = FT(0.5) * (vRR + vLL) - FT(1.0) /(FT(2.0) * cS) * (pRR - pLL) / RhoM
    pM = -FT(0.5) * Phys.Cpd * RhoM * ThM * (AuxR[3] - AuxL[3])
    if vM > FT(0)
      FL[RhoPos] = vM * VLL[RhoPos]
      FL[uPos] = vM * VLL[uPos] + n1 * dp
      FL[vPos] = vM * VLL[vPos] + n2 * dp
      FL[wPos] = vM * VLL[wPos] + n3 * dp
      FL[ThPos] = vM * VLL[ThPos]
    else
      FL[RhoPos] = vM * VRR[RhoPos]
      FL[uPos] = vM * VRR[uPos] + n1 * dp
      FL[vPos] = vM * VRR[vPos] + n2 * dp
      FL[wPos] = vM * VRR[wPos] + n3 * dp
      FL[ThPos] = vM * VRR[ThPos]
    end

    FR[RhoPos] = FL[RhoPos]
    FR[uPos] = FL[uPos] - pM * n1
    FR[vPos] = FL[vPos] - pM * n2
    FR[wPos] = FL[wPos] - pM * n3
    FR[ThPos] = FL[ThPos]

    #FL[RhoPos] = FL[RhoPos]
    FL[uPos] = FL[uPos] + pM * n1
    FL[vPos] = FL[vPos] + pM * n2
    FL[wPos] = FL[wPos] + pM * n3
    #FL[ThPos] = FL[ThPos]

  end
  return RiemannByLMARSNonLin!
end

function (::RiemannExLMARS)(Param,Phys,RhoPos,uPos,vPos,wPos,ThPos,pPos)
  @inline function RiemannByLMARSNonLin!(F,VLL,VRR,AuxL,AuxR,n1,n2,n3)

    FT = eltype(F)
    cS = Param.cS
    pLL = AuxL[pPos]
    pRR = AuxR[pPos]
    ExpLL = AuxL[3]
    ExpRR = AuxR[3]
    RhoM = FT(0.5) * (VLL[RhoPos] + VRR[RhoPos])
    RhoThM = FT(0.5) * (VLL[ThPos] + VRR[ThPos])
    vLL = (VLL[uPos] * n1 + VLL[vPos] * n2 + VLL[wPos] * n3) / VLL[RhoPos]
    vRR = (VRR[uPos] * n1 + VRR[vPos] * n2 + VRR[wPos] * n3) / VRR[RhoPos]
#   pM = FT(0.5) * (pLL + pRR) - FT(0.5) * cS * RhoM * (vRR - vLL)
    vM = FT(0.5) * (vRR + vLL) - FT(1.0) /(FT(2.0) * cS) * (pRR - pLL) / RhoM
    pM = -FT(0.5) * cS * RhoM * (vRR - vLL)
    pM += 0.5 * Phys.Cpd * RhoThM * (ExpRR - ExpLL)  
    if vM > FT(0)
      F[RhoPos] = vM * VLL[RhoPos]
      F[uPos] = vM * VLL[uPos] + n1 * pM
      F[vPos] = vM * VLL[vPos] + n2 * pM
      F[wPos] = vM * VLL[wPos] + n3 * pM
      F[ThPos] = vM * VLL[ThPos]
    else
      F[RhoPos] = vM * VRR[RhoPos]
      F[uPos] = vM * VRR[uPos] + n1 * pM
      F[vPos] = vM * VRR[vPos] + n2 * pM
      F[wPos] = vM * VRR[wPos] + n3 * pM
      F[ThPos] = vM * VRR[ThPos]
    end
  end
  return RiemannByLMARSNonLin!
end

function (::RiemannExPLMARS)(Param,Phys,RhoPos,uPos,vPos,wPos,ThPos,pPos)
  @inline function RiemannByLMARSNonLin!(F,VLL,VRR,AuxL,AuxR,n1,n2,n3)

    FT = eltype(F)
    cS = Param.cS
    pLL = AuxL[pPos]
    pRR = AuxR[pPos]
    ExpLL = AuxL[3]
    ExpRR = AuxR[3]
    RhoM = FT(0.5) * (VLL[RhoPos] + VRR[RhoPos])
    RhoThM = FT(0.5) * (VLL[ThPos] + VRR[ThPos])
    vLL = (VLL[uPos] * n1 + VLL[vPos] * n2 + VLL[wPos] * n3) / VLL[RhoPos]
    vRR = (VRR[uPos] * n1 + VRR[vPos] * n2 + VRR[wPos] * n3) / VRR[RhoPos]
    pM = FT(0.5) * Phys.Cpd * (VLL[ThPos] * ExpLL + VRR[ThPos]* ExpRR) - FT(0.5) * cS * RhoM * (vRR - vLL)
    vM = FT(0.5) * (vRR + vLL) - FT(1.0) /(FT(2.0) * cS) * (pRR - pLL) / RhoM
    if vM > FT(0)
      F[RhoPos] = vM * VLL[RhoPos]
      F[uPos] = vM * VLL[uPos] + n1 * pM
      F[vPos] = vM * VLL[vPos] + n2 * pM
      F[wPos] = vM * VLL[wPos] + n3 * pM
      F[ThPos] = vM * VLL[ThPos]
    else
      F[RhoPos] = vM * VRR[RhoPos]
      F[uPos] = vM * VRR[uPos] + n1 * pM
      F[vPos] = vM * VRR[vPos] + n2 * pM
      F[wPos] = vM * VRR[wPos] + n3 * pM
      F[ThPos] = vM * VRR[ThPos]
    end
  end
  return RiemannByLMARSNonLin!
end


function (::RiemannLMARSFast)(Param,Phys,RhoPos,uPos,vPos,wPos,RhoThPos,dpdRhoThPos,ThPos)
  @inline function RiemannByLMARSSemi!(F,VLL,VRR,AuxL,AuxR,n1,n2,n3)

    FT = eltype(F)
    cS = Param.cS
    pLL = AuxL[dpdRhoThPos] * VLL[RhoThPos]
    pRR = AuxR[dpdRhoThPos] * VRR[RhoThPos]
    ThM = FT(0.5) * (AuxL[ThPos] + AuxR[ThPos])
    RhovLL = (VLL[uPos] * n1 + VLL[vPos] * n2 + VLL[wPos] * n3) 
    RhovRR = (VRR[uPos] * n1 + VRR[vPos] * n2 + VRR[wPos] * n3) 
    pM = FT(0.5) * (pLL + pRR) - FT(0.5) * cS * (RhovRR - RhovLL)
    RhovM = FT(0.5) * (RhovRR + RhovLL) - FT(1.0) /(FT(2.0) * cS) * (pRR - pLL)
    F[uPos] = n1 * pM
    F[vPos] = n2 * pM
    F[wPos] = n3 * pM
    F[RhoPos] = RhovM 
    F[RhoThPos] = RhovM * ThM
  end
  return RiemannByLMARSSemi!
end

function (::RiemannLMARSLinFast)(Param,Phys,RhoPos,uPos,vPos,wPos,RhoThPos,dpdRhoThPos,ThPos)
  @inline function RiemannByLMARSSemi!(F,VLL,VRR,AuxL,AuxR,n1,n2,n3)

    FT = eltype(F)
    cS = Param.cS
    pLL = AuxL[dpdRhoThPos] * VLL[RhoThPos]
    pRR = AuxR[dpdRhoThPos] * VRR[RhoThPos]
    vLL = (VLL[uPos] * n1 + VLL[vPos] * n2 + VLL[wPos] * n3) 
    vRR = (VRR[uPos] * n1 + VRR[vPos] * n2 + VRR[wPos] * n3) 
    pM = FT(0.5) * (pLL + pRR) - FT(0.5) * cS * (vRR - vLL)
    vM = FT(0.5) * (vRR + vLL) - FT(1.0) /(FT(2.0) * cS) * (pRR - pLL) 
    F[uPos] = n1 * pM
    F[vPos] = n2 * pM
    F[wPos] = n3 * pM
    F[RhoPos] = vM 
    F[RhoThPos] = FT(0.5) * F[RhoPos] * (AuxL[ThPos] + AuxR[ThPos])
  end
  return RiemannByLMARSSemi!
end


function (::RiemannLMARSFast)(Param,Phys,RhoPos,uPos,vPos,wPos,RhoThPos,pPos)
  @inline function RiemannByLMARSLin!(F,VLL,VRR,AuxL,AuxR,n1,n2,n3)

    FT = eltype(F)
    cS = Param.cS
    pLL = AuxL[pPos]
    pRR = AuxR[pPos]
    RhoM = FT(0.5) * (VLL[RhoPos] + VRR[RhoPos])
    vLL = (VLL[uPos] * n1 + VLL[vPos] * n2 + VLL[wPos] * n3) / VLL[RhoPos]
    vRR = (VRR[uPos] * n1 + VRR[vPos] * n2 + VRR[wPos] * n3) / VRR[RhoPos]
    pM = FT(0.5) * (pLL + pRR) - FT(0.5) * cS * RhoM * (vRR - vLL)
    vM = FT(0.5) * (vRR + vLL) - FT(1.0) /(FT(2.0) * cS) * (pRR - pLL) / RhoM
    if vM > FT(0)
      F[RhoPos] = vM * VLL[RhoPos]
      F[uPos] = n1 * pM
      F[vPos] = n2 * pM
      F[wPos] = n3 * pM
      F[RhoThPos] = vM * VLL[RhoThPos]
    else
      F[RhoPos] = vM * VRR[RhoPos]
      F[uPos] = n1 * pM
      F[vPos] = n2 * pM
      F[wPos] = n3 * pM
      F[RhoThPos] = vM * VRR[RhoThPos]
    end
  end
  return RiemannByLMARSLin!
end

Base.@kwdef struct RiemannBoussinesqLMARS <: RiemannSolver end

function (::RiemannBoussinesqLMARS)(Param,pPos,uPos,vPos,wPos,bPos)
  @inline function RiemannByLMARS(F,VLL,VRR,AuxL,AuxR,n1,n2,n3)

    FT = eltype(F)
    cS = Param.cS
    pLL = VLL[pPos]
    pRR = VRR[pPos]
    U = n1 * Param.U
    vLL = VLL[uPos] * n1 + VLL[vPos] * n2 + VLL[wPos] * n3 
    vRR = VRR[uPos] * n1 + VRR[vPos] * n2 + VRR[wPos] * n3
    pM = FT(0.5) * (pLL + pRR) - FT(0.5) * cS * (vRR - vLL)
    vM = FT(0.5) * (vRR + vLL) - FT(1.0) /(FT(2.0) * cS) * (pRR - pLL) 
    if Param.U > 0
      F[pPos] = U * VLL[pPos] + Param.cS^2*vM 
      F[uPos] = U * VLL[uPos] + n1 * pM
      F[vPos] = U * VLL[vPos] + n2 * pM
      F[wPos] = U * VLL[wPos] + n3 * pM
      F[bPos] = U * VLL[bPos]
    else
      F[pPos] = U * VRR[pPos] + Param.cS^2*vM 
      F[uPos] = U * VRR[uPos] + n1 * pM 
      F[vPos] = U * VRR[vPos] + n2 * pM
      F[wPos] = U * VRR[wPos] + n3 * pM
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






