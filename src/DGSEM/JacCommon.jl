################################################################################
# Common parts of the vertically implicit HEVI solvers (JacHDG, JacTwoSidedSplit,
# JacOneSided)
#
# The vertical Jacobian is formulated in the third contravariant direction
# (arXiv:2610.15775, Sec. 3.3): at every node the vertical momentum is
# rhow = nu . rhov with the unit vector nu = Ja^3 / n, n = |Ja^3|
# (Metric.NVV and Metric.VolSurfVV). The horizontal momentum is not coupled.
#
# Linearization state per node: theta = rhotheta / rho and
# alpha = dp / drhotheta, computed from the current state U (Jac!).
################################################################################

@kernel inbounds = true function ThetaAlphaKernel!(Th,Al,@Const(U),kfac,kexp,Rdp0)
  i,iz,ID = @index(Global, NTuple)
  ND = @uniform @ndrange()[3]
  if ID <= ND
    RhoTh = U[i,iz,ID,5]
    Th[i,iz,ID] = RhoTh / U[i,iz,ID,1]
    Al[i,iz,ID] = kfac * (Rdp0 * RhoTh)^kexp
  end
end

function FillThetaAlpha!(Th,Al,U,Phys)
  backend = get_backend(Th)
  FT = eltype(Th)
  M, nz, ND = size(Th)
  kappa = FT(Phys.kappa)
  kexp  = kappa / (one(FT) - kappa)
  kfac  = FT(Phys.Rd) / (one(FT) - kappa)
  Rdp0  = FT(Phys.Rd) / FT(Phys.p0)
  NDG = max(1,min(div(256,M*nz),ND))
  KThetaAlphaKernel! = ThetaAlphaKernel!(backend,(M,nz,NDG))
  KThetaAlphaKernel!(Th,Al,U,kfac,kexp,Rdp0;ndrange=(M,nz,ND))
end

# Vector invariant form: Th = rhotheta, Al = dPi / drhotheta, Rm = rho,
# Ps = cp theta, the Exner pressure Pi = (Rd rhotheta / p0)^(Rd/Cvd) and the
# sound speed Cs = sqrt(gamma p / rho)
@kernel inbounds = true function CoeffVIKernel!(Th,Al,Rm,Ps,Pi,Cs,@Const(U),Cpd,ex,Rdp0,p0,Gamma)
  i,iz,ID = @index(Global, NTuple)
  ND = @uniform @ndrange()[3]
  if ID <= ND
    Rho = U[i,iz,ID,1]
    RhoTh = U[i,iz,ID,5]
    PiLoc = (Rdp0 * RhoTh)^ex
    Th[i,iz,ID] = RhoTh
    Al[i,iz,ID] = ex * PiLoc / RhoTh
    Rm[i,iz,ID] = Rho
    Ps[i,iz,ID] = Cpd * RhoTh / Rho
    Pi[i,iz,ID] = PiLoc
    Cs[i,iz,ID] = sqrt(Gamma * p0 * (Rdp0 * RhoTh)^Gamma / Rho)
  end
end

function FillCoeffVI!(Th,Al,Rm,Ps,Pi,Cs,U,Phys)
  backend = get_backend(Th)
  FT = eltype(Th)
  M, nz, ND = size(Th)
  Cpd  = FT(Phys.Cpd)
  ex   = FT(Phys.Rd) / FT(Phys.Cvd)
  Rdp0 = FT(Phys.Rd) / FT(Phys.p0)
  p0   = FT(Phys.p0)
  Gamma = FT(Phys.Cpd) / FT(Phys.Cvd)
  NDG = max(1,min(div(256,M*nz),ND))
  KCoeffVIKernel! = CoeffVIKernel!(backend,(M,nz,NDG))
  KCoeffVIKernel!(Th,Al,Rm,Ps,Pi,Cs,U,Cpd,ex,Rdp0,p0,Gamma;ndrange=(M,nz,ND))
end

################################################################################
# Interface fluxes, variable order (rho, rhotheta, rhow)
################################################################################

# Derivatives of the LMARS flux (rhow, theta rhow, p) with the constant sound
# speed cS with respect to the lower (L) and the upper (R) state
@inline function LMARSdL(ThL,ThR,AlL,cS)
  FT = typeof(ThL)
  ThB = FT(0.5) * (ThL + ThR)
  return SMatrix{3,3,FT}(0, 0, 0,
    AlL / (2 * cS), ThB * AlL / (2 * cS), FT(0.5) * AlL,
    FT(0.5), FT(0.5) * ThL, FT(0.5) * cS)
end

@inline function LMARSdR(ThL,ThR,AlR,cS)
  FT = typeof(ThL)
  ThB = FT(0.5) * (ThL + ThR)
  return SMatrix{3,3,FT}(0, 0, 0,
    -AlR / (2 * cS), -ThB * AlR / (2 * cS), FT(0.5) * AlR,
    FT(0.5), FT(0.5) * ThR, -FT(0.5) * cS)
end

# Derivatives of the vector invariant interface flux with respect to the lower
# (L) and the upper (R) state; rows (rho, rhotheta, w), columns (rho, rhotheta, w).
# Side = 1: flux of the lower element (FL), Side = 2: of the upper element (FR).
@inline function VIdL(ThL,ThR,RmL,RmR,AlL,PsL,PsR,c,Side)
  FT = typeof(ThL)
  RmB = FT(0.5) * (RmL + RmR)
  ThB = FT(0.5) * (ThL + ThR)
  PsB = FT(0.5) * (PsL + PsR)
  dw = Side == 1 ? (PsL - FT(0.5) * PsB) * AlL : FT(0.5) * PsB * AlL
  return SMatrix{3,3,FT}(0, 0, 0,
    RmB * PsB * AlL / (2 * c), ThB * PsB * AlL / (2 * c), dw,
    FT(0.5) * RmB, FT(0.5) * ThB, FT(0.5) * c)
end

@inline function VIdR(ThL,ThR,RmL,RmR,AlR,PsL,PsR,c,Side)
  FT = typeof(ThL)
  RmB = FT(0.5) * (RmL + RmR)
  ThB = FT(0.5) * (ThL + ThR)
  PsB = FT(0.5) * (PsL + PsR)
  dw = Side == 1 ? FT(0.5) * PsB * AlR : (PsR - FT(0.5) * PsB) * AlR
  return SMatrix{3,3,FT}(0, 0, 0,
    -RmB * PsB * AlR / (2 * c), -ThB * PsB * AlR / (2 * c), dw,
    FT(0.5) * RmB, FT(0.5) * ThB, -FT(0.5) * c)
end

# Reflecting wall, flux (0, 0, p -+ cS rhow) of the interior state;
# sgn = -1 at the bottom, +1 at the top
@inline function WallFlux(Al,cS,sgn)
  FT = typeof(Al)
  return SMatrix{3,3,FT}(0, 0, 0, 0, 0, Al, 0, 0, sgn * cS)
end

################################################################################
# Element data for the contravariant formulation
#   rJ  : 1 / J
#   nS  : n = |Ja^3|
#   Gr  : (n / J) d phi / d xi3 (gravity in the contravariant direction)
################################################################################

@inline function LoadElementC(Th,Al,J,Surf,Geo,DW,iz,ID,::Val{M}) where {M}
  FT = eltype(Th)
  ThL = MVector{M,FT}(undef)
  AlL = MVector{M,FT}(undef)
  rJ  = MVector{M,FT}(undef)
  nS  = MVector{M,FT}(undef)
  Gr  = MVector{M,FT}(undef)
  @unroll for i = 1 : M
    ThL[i] = Th[i,iz,ID]
    AlL[i] = Al[i,iz,ID]
    rJ[i]  = FT(1) / J[i,iz,ID]
    nS[i]  = Surf[i,iz,ID]
  end
  @unroll for i = 1 : M
    acc = FT(0)
    @unroll for k = 1 : M
      acc += DW[i,k] * (Geo[k,iz,ID] - Geo[i,iz,ID])
    end
    Gr[i] = rJ[i] * nS[i] * acc
  end
  return ThL, AlL, rJ, nS, Gr
end

# Split form gravity: (n_i / J_i) G_ik with
# G_ik = 1/2 DW_ik (phi_k - phi_i) + delta_ik 1/2 sum_m DW_im (phi_m - phi_i)
@inline function GravSplit(Geo,DW,rJ,nS,i,k,iz,ID,::Val{M}) where {M}
  FT = eltype(Geo)
  Phii = Geo[i,iz,ID]
  val = FT(0.5) * DW[i,k] * (Geo[k,iz,ID] - Phii)
  if i == k
    acc = FT(0)
    @unroll for m = 1 : M
      acc += DW[i,m] * (Geo[m,iz,ID] - Phii)
    end
    val += FT(0.5) * acc
  end
  return rJ[i] * nS[i] * val
end

# Vertical momentum component of the right hand side
@inline NuDot(NV,v,i,iz,ID) =
  NV[1,i,iz,ID] * v[i,iz,ID,2] + NV[2,i,iz,ID] * v[i,iz,ID,3] + NV[3,i,iz,ID] * v[i,iz,ID,4]

# Write the solution: rho, rhotheta and the vertical momentum xw of the
# vertical system, the horizontal momentum is facH * (I - nu nu^T) v
@inline function StoreSolution!(k,v,NV,i,iz,ID,xr,xt,xw,facH)
  vw = NuDot(NV,v,i,iz,ID)
  bb = xw - facH * vw
  k[i,iz,ID,2] = facH * v[i,iz,ID,2] + NV[1,i,iz,ID] * bb
  k[i,iz,ID,3] = facH * v[i,iz,ID,3] + NV[2,i,iz,ID] * bb
  k[i,iz,ID,4] = facH * v[i,iz,ID,4] + NV[3,i,iz,ID] * bb
  k[i,iz,ID,1] = xr
  k[i,iz,ID,5] = xt
end

# The geopotential is static, it is copied once from the auxiliary variables
function InitGeo!(Jac,Cache,DG,Global)
  if Jac.grav_do
    @views Jac.Geo .= Cache.Aux[:,:,1:DG.NumI,Global.Model.GPAuxPos]
    Jac.grav_do = false
  end
end
