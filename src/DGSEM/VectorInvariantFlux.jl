################################################################################
# Vector invariant formulation with potential temperature,
# Artiano, Ricardo, Knoth, Spichtinger, Ranocha, arXiv:2610.05820
#
#   rho_t   + div(rho v)                                    = 0
#   v_t     + (curl v) x v + grad(|v|^2/2) + cp theta grad(Pi) + grad(phi) = 0
#   (rho theta)_t + div(rho theta v)                        = 0
#
# Prognostic variables (rho, u, v, w, rho theta) at the positions 1:5.
# In the reference direction m with the metric row n = J a^m the
# nonconservative terms of the momentum row i are (q = n . v)
#   q d_m v_i - n_i (v . d_m v) + n_i (cp theta d_m Pi + d_m phi),
# which reduces to H^1, H^3 of the paper in two dimensions. The two-point
# volume function returns f + 1/2 H^num [[g]] (antisymmetric jumps vanish on
# the diagonal, so the split-form volume kernels are used unchanged).
################################################################################

# Numerical theta: LogMean 1/{{1/theta}}_log (entropy conservative,
# well-balanced for T = const and theta = const) or ArithMean {{theta}}
# (total energy conservative, well-balanced for theta = const)
@inline ThetaNum(::LogMean, ThL, ThR) = inv_ln_mean(one(ThL) / ThL, one(ThR) / ThR)
@inline ThetaNum(::ArithMean, ThL, ThR) = typeof(ThL)(0.5) * (ThL + ThR)

# Well-balanced, entropy and total energy conservative volume flux (Eq. 58),
# primitive states (rho, u, v, w, theta)
Base.@kwdef struct VectorInvariantFlux <: AverageFlux end

function (::VectorInvariantFlux)(::Val{ExPos}, ::Val{GPPos}, Cpd, ThMean = LogMean()) where {ExPos, GPPos}
  @inline function FluxNonLinAver(VL, VR, AuxL, AuxR, mL, mR)
    FT = eltype(VL)
    n1 = FT(0.5) * (mL[1] + mR[1])
    n2 = FT(0.5) * (mL[2] + mR[2])
    n3 = FT(0.5) * (mL[3] + mR[3])
    RhoL = VL[1]
    RhoR = VR[1]
    RhoAv = FT(0.5) * (RhoL + RhoR)
    # {{rho v}}
    m1 = FT(0.5) * (RhoL * VL[2] + RhoR * VR[2])
    m2 = FT(0.5) * (RhoL * VL[3] + RhoR * VR[3])
    m3 = FT(0.5) * (RhoL * VL[4] + RhoR * VR[4])
    fRho = n1 * m1 + n2 * m2 + n3 * m3
    ThNum = ThetaNum(ThMean, VL[5], VR[5])
    fRhoTh = fRho * ThNum
    KAv = FT(0.25) * (VL[2] * VL[2] + VL[3] * VL[3] + VL[4] * VL[4] +
      VR[2] * VR[2] + VR[3] * VR[3] + VR[4] * VR[4])
    # nonconservative part with v^num = {{rho v}} / {{rho}}
    v1 = m1 / RhoAv
    v2 = m2 / RhoAv
    v3 = m3 / RhoAv
    du = VR[2] - VL[2]
    dv = VR[3] - VL[3]
    dw = VR[4] - VL[4]
    q = n1 * v1 + n2 * v2 + n3 * v3
    vdv = v1 * du + v2 * dv + v3 * dw
    PG = Cpd * ThNum * (AuxR[ExPos] - AuxL[ExPos]) + (AuxR[GPPos] - AuxL[GPPos])
    s = vdv - PG
    return SVector(fRho,
      n1 * KAv + FT(0.5) * (q * du - n1 * s),
      n2 * KAv + FT(0.5) * (q * dv - n2 * s),
      n3 * KAv + FT(0.5) * (q * dw - n3 * s),
      fRhoTh)
  end
  return FluxNonLinAver
end

# Entropy stable and well-balanced interface flux (Eq. 72) with the
# nonconservative surface term: FL = f + 1/2 H [[g]], FR = f - 1/2 H [[g]].
# Stored states (rho, u, v, w, rho theta), unit normal n from L to R.
# cFix > 0: fixed sound speed cFix instead of {{c}} = {{sqrt(gamma p / rho)}}.
Base.@kwdef struct RiemannVectorInvariant <: RiemannSolver end

function (::RiemannVectorInvariant)(::Val{pPos}, ::Val{ExPos}, ::Val{GPPos}, Cpd, Gamma,
  ThMean = LogMean(), cFix = zero(Cpd)) where {pPos, ExPos, GPPos}
  @inline function RiemannVI!(FL, FR, VLL, VRR, AuxL, AuxR, n1, n2, n3)
    FT = eltype(FL)
    RhoL = VLL[1]
    RhoR = VRR[1]
    uL = VLL[2]; vL = VLL[3]; wL = VLL[4]
    uR = VRR[2]; vR = VRR[3]; wR = VRR[4]
    ThL = VLL[5] / RhoL
    ThR = VRR[5] / RhoR
    RhoAv = FT(0.5) * (RhoL + RhoR)
    if cFix > FT(0)
      cAv = cFix
    else
      cAv = FT(0.5) * (sqrt(Gamma * AuxL[pPos] / RhoL) + sqrt(Gamma * AuxR[pPos] / RhoR))
    end
    ThNum = ThetaNum(ThMean, ThL, ThR)
    PG = Cpd * ThNum * (AuxR[ExPos] - AuxL[ExPos]) + (AuxR[GPPos] - AuxL[GPPos])

    Vn = FT(0.5) * ((uL + uR) * n1 + (vL + vR) * n2 + (wL + wR) * n3)
    VNum = Vn - PG / (FT(2) * cAv)
    sg = sign(VNum)
    fRho = (RhoAv - FT(0.5) * (RhoR - RhoL) * sg) * VNum
    fRhoTh = fRho * (ThNum - min(ThL, ThR) / (ThL + ThR) * (ThR - ThL) * sg)

    KAv = FT(0.25) * (uL * uL + vL * vL + wL * wL + uR * uR + vR * vR + wR * wR)
    dmu = RhoR * uR - RhoL * uL
    dmv = RhoR * vR - RhoL * vL
    dmw = RhoR * wR - RhoL * wL
    a = KAv - cAv / (FT(2) * RhoAv) * (dmu * n1 + dmv * n2 + dmw * n3)
    b = abs(Vn) / (FT(2) * RhoAv)
    fu = a * n1 - b * dmu
    fv = a * n2 - b * dmv
    fw = a * n3 - b * dmw

    # 1/2 H^num [[g]] with v^num = {{rho v}} / {{rho}}
    v1 = FT(0.5) * (RhoL * uL + RhoR * uR) / RhoAv
    v2 = FT(0.5) * (RhoL * vL + RhoR * vR) / RhoAv
    v3 = FT(0.5) * (RhoL * wL + RhoR * wR) / RhoAv
    du = uR - uL
    dv = vR - vL
    dw = wR - wL
    q = n1 * v1 + n2 * v2 + n3 * v3
    s = v1 * du + v2 * dv + v3 * dw - PG
    hu = FT(0.5) * (q * du - n1 * s)
    hv = FT(0.5) * (q * dv - n2 * s)
    hw = FT(0.5) * (q * dw - n3 * s)

    FL[1] = fRho
    FL[2] = fu + hu
    FL[3] = fv + hv
    FL[4] = fw + hw
    FL[5] = fRhoTh
    FR[1] = fRho
    FR[2] = fu - hu
    FR[3] = fv - hv
    FR[4] = fw - hw
    FR[5] = fRhoTh
  end
  return RiemannVI!
end
