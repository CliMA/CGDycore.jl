################################################################################
# Vertical HEVI solver with one-sided interface condensation
# (Artiano, Rueda-Ramirez, Knoth, Ranocha, Gassner, arXiv:2610.15775, Sec. 4.2.2)
#
# Linear system per column: (fac I - J) x = b with the HEVI Jacobian J of the
# paper: strong-form DGSEM with the weak derivative DW in the contravariant
# vertical direction (see JacCommon.jl), frozen theta, linearized pressure
# dp = alpha drhotheta, LMARS interface flux with the constant sound speed cS,
# pointwise gravity, reflecting walls at bottom and top.
#
# Per element e the unknowns are split into
#   y_e : rho, rhotheta, rhow at node M          (reduced unknowns)
#   x_e : rho, rhotheta, rhow at nodes 1 ... M-1 (eliminated element-wise)
# Only node 1 of element e couples to y_{e-1} and only y_e couples to node 1
# of element e+1. Eliminating x_e gives a block tridiagonal system for y with
# 3x3 blocks, which is solved with the block Thomas algorithm.
#
# Element solve A_e x = r: rho and rhotheta are eliminated analytically, the
# Schur complement for rhow is factored with LUFull! (stored in SLU).
#
# Vector invariant form (VI = true, variables rho, rhotheta, w): node
# coefficients of JacHDG, interface blocks VIdL/VIdR, and in the spirit of the
# pointwise gravity the pointwise buoyancy Bt (delta rhotheta - theta delta rho)
# with Bt = (n / J) cp (D Pi) / rho.
#
# Memory layouts
#   Geo, Th, Al, Rm, Ps, Pi, Cs     : (M, nz, NumI)
#   SLU                             : (M-1, M-1, nz, NumI)
#   Dself, Dprev, Lo, Up, Binv, Cp  : (NumI, 3, 3, nz)
#   rsA, rsB, rs                    : (NumI, 3, nz)
# Each entry has exactly one writer, contributions to the neighbour element
# are stored separately (Dprev, rsB) and gathered in the column kernels.
################################################################################

mutable struct JacOneSided{FT<:AbstractFloat,
                           AT3<:AbstractArray,
                           AT4<:AbstractArray}
  CompTri::Bool
  grav_do::Bool
  M::Int
  nz::Int
  fac::FT
  cS::FT
  Geo::AT3
  Th::AT3
  Al::AT3
  Rm::AT3
  Ps::AT3
  Pi::AT3
  Cs::AT3
  VI::Bool
  SLU::AT4
  Dself::AT4
  Dprev::AT4
  Lo::AT4
  Up::AT4
  Binv::AT4
  Cp::AT4
  rsA::AT3
  rsB::AT3
  rs::AT3
end

function JacOneSided(backend,FT,M,nz,DG)
  CompTri = false
  grav_do = true
  fac = zero(FT)
  cS = zero(FT)
  NumI = DG.NumI
  Geo = KernelAbstractions.zeros(backend,FT,M,nz,NumI)
  Th = KernelAbstractions.zeros(backend,FT,M,nz,NumI)
  Al = KernelAbstractions.zeros(backend,FT,M,nz,NumI)
  Rm = KernelAbstractions.ones(backend,FT,M,nz,NumI)
  Ps = KernelAbstractions.ones(backend,FT,M,nz,NumI)
  Pi = KernelAbstractions.zeros(backend,FT,M,nz,NumI)
  Cs = KernelAbstractions.zeros(backend,FT,M,nz,NumI)
  SLU = KernelAbstractions.zeros(backend,FT,M-1,M-1,nz,NumI)
  Dself = KernelAbstractions.zeros(backend,FT,NumI,3,3,nz)
  Dprev = KernelAbstractions.zeros(backend,FT,NumI,3,3,nz)
  Lo = KernelAbstractions.zeros(backend,FT,NumI,3,3,nz)
  Up = KernelAbstractions.zeros(backend,FT,NumI,3,3,nz)
  Binv = KernelAbstractions.zeros(backend,FT,NumI,3,3,nz)
  Cp = KernelAbstractions.zeros(backend,FT,NumI,3,3,nz)
  rsA = KernelAbstractions.zeros(backend,FT,NumI,3,nz)
  rsB = KernelAbstractions.zeros(backend,FT,NumI,3,nz)
  rs = KernelAbstractions.zeros(backend,FT,NumI,3,nz)
  return JacOneSided{FT,
                     typeof(Th),
                     typeof(SLU)}(
    CompTri,
    grav_do,
    M,
    nz,
    fac,
    cS,
    Geo,
    Th,
    Al,
    Rm,
    Ps,
    Pi,
    Cs,
    false,
    SLU,
    Dself,
    Dprev,
    Lo,
    Up,
    Binv,
    Cp,
    rsA,
    rsB,
    rs,
  )
end

# Contribution of the lower surface to the rows of node 1 (3x3 block on
# rho_1, rhotheta_1, rhow_1 of the own element)
@inline function CornerBottom(Th1,Al1,ThMm,iz,sig1,cS)
  if iz > 1
    return -sig1 * LMARSdR(ThMm,Th1,Al1,cS)
  else
    return -sig1 * WallFlux(Al1,cS,-1)
  end
end

@inline CornerBottom(Th1,Al1,ThMm,iz,sig1,cS,Th,Rm,Ps,Cs,ID,::Val{M},::Val{false}) where {M} =
  CornerBottom(Th1,Al1,ThMm,iz,sig1,cS)

@inline function CornerBottom(Th1,Al1,ThMm,iz,sig1,cS,Th,Rm,Ps,Cs,ID,::Val{M},::Val{true}) where {M}
  FT = typeof(Th1)
  if iz > 1
    c = FT(0.5) * (Cs[M,iz-1,ID] + Cs[1,iz,ID])
    return -sig1 * VIdR(ThMm,Th1,Rm[M,iz-1,ID],Rm[1,iz,ID],Al1,Ps[M,iz-1,ID],Ps[1,iz,ID],c,2)
  else
    return -sig1 * WallFlux(Ps[1,iz,ID] * Al1,Cs[1,iz,ID],-1)
  end
end

# Interface block of the face below element iz (lower element iz-1, upper iz):
# derivative of the flux of the lower element (Side 1) or of the upper
# element (Side 2) with respect to the lower (L) or upper (R) state
@inline function FaceBelow(LR,Side,Th,Al,ThL1,AlL1,iz,ID,cS,Rm,Ps,Cs,::Val{M},::Val{false}) where {M}
  return LR == 1 ? LMARSdL(Th[M,iz-1,ID],ThL1,Al[M,iz-1,ID],cS) : LMARSdR(Th[M,iz-1,ID],ThL1,AlL1,cS)
end

@inline function FaceBelow(LR,Side,Th,Al,ThL1,AlL1,iz,ID,cS,Rm,Ps,Cs,::Val{M},::Val{true}) where {M}
  FT = eltype(Th)
  c = FT(0.5) * (Cs[M,iz-1,ID] + Cs[1,iz,ID])
  a = (Th[M,iz-1,ID],ThL1,Rm[M,iz-1,ID],Rm[1,iz,ID])
  return LR == 1 ? VIdL(a...,Al[M,iz-1,ID],Ps[M,iz-1,ID],Ps[1,iz,ID],c,Side) :
    VIdR(a...,AlL1,Ps[M,iz-1,ID],Ps[1,iz,ID],c,Side)
end

# Top block of element iz (flux of the lower element, Side 1, derivative w.r.t. own state)
@inline function FaceTop(ThLM,AlLM,Th,iz,ID,nz,cS,Rm,Ps,Cs,::Val{M},::Val{false}) where {M}
  return iz < nz ? LMARSdL(ThLM,Th[1,iz+1,ID],AlLM,cS) : WallFlux(AlLM,cS,1)
end

@inline function FaceTop(ThLM,AlLM,Th,iz,ID,nz,cS,Rm,Ps,Cs,::Val{M},::Val{true}) where {M}
  FT = eltype(Th)
  if iz < nz
    c = FT(0.5) * (Cs[M,iz,ID] + Cs[1,iz+1,ID])
    return VIdL(ThLM,Th[1,iz+1,ID],Rm[M,iz,ID],Rm[1,iz+1,ID],AlLM,Ps[M,iz,ID],Ps[1,iz+1,ID],c,1)
  else
    return WallFlux(Ps[M,iz,ID] * AlLM,Cs[M,iz,ID],1)
  end
end

# Pointwise gravity / buoyancy on rho (row i)
@inline GravPW(Gr,Pi,Rm,Ps,Th,DW,rJ,nS,i,iz,ID,::Val{M},::Val{false}) where {M} = Gr[i]
@inline GravPW(Gr,Pi,Rm,Ps,Th,DW,rJ,nS,i,iz,ID,::Val{M},::Val{true}) where {M} =
  -BuoyPW(Pi,Rm,Ps,Th,DW,rJ,nS,i,iz,ID,Val(M)) * Th[i,iz,ID] / Rm[i,iz,ID]

# Bt_i = (n / J) cp (D Pi)_i / rho_i
@inline function BuoyPW(Pi,Rm,Ps,Th,DW,rJ,nS,i,iz,ID,::Val{M}) where {M}
  FT = eltype(Pi)
  acc = FT(0)
  @unroll for m = 1 : M
    acc += DW[i,m] * (Pi[m,iz,ID] - Pi[i,iz,ID])
  end
  return rJ[i] * nS[i] * Ps[i,iz,ID] / Th[i,iz,ID] * acc
end

# Pressure coupling of row i to rhotheta_j
@inline PresW1(Al,Pi,Rm,Ps,Th,DW,rJ,nS,i,j,iz,ID,::Val{M},::Val{false}) where {M} =
  rJ[i] * DW[i,j] * nS[j] * Al[j]
@inline function PresW1(Al,Pi,Rm,Ps,Th,DW,rJ,nS,i,j,iz,ID,::Val{M},::Val{true}) where {M}
  val = rJ[i] * nS[i] * Ps[i,iz,ID] * DW[i,j] * Al[j]
  if i == j
    val += BuoyPW(Pi,Rm,Ps,Th,DW,rJ,nS,i,iz,ID,Val(M))
  end
  return val
end

# (mass flux operator)_ij
@inline MassOp1(Rm,DW,rJ,nS,i,j,iz,ID,::Val{false}) = rJ[i] * DW[i,j] * nS[j]
@inline MassOp1(Rm,DW,rJ,nS,i,j,iz,ID,::Val{true}) = rJ[i] * DW[i,j] * nS[j] * Rm[j,iz,ID]

################################################################################
# Element solve A_e x = r, x and r are (M-1) x 3 (columns rho, rhotheta, rhow)
################################################################################

@inline function SolveElement!(x, r, SLU, iz, ID, Th, Al, rJ, nS, Gr, K1, DW, fac, Co, ::Val{M}, ::Val{VI}) where {M, VI}
  FT = eltype(x)
  ThA, Rm, Ps, Pi = Co
  Mm = M - 1
  invfac = FT(1) / fac
  d21 = FT(1) / (fac + K1[2,2])
  # rhotheta without the rhow part: q2 = A22^{-1} r2
  @unroll for i = 1 : Mm
    x[i,2] = r[i,2] * invfac
  end
  x[1,2] = r[1,2] * d21
  # rho without the rhow part: q1 = A11^{-1} (r1 - A12 q2)
  @unroll for i = 1 : Mm
    x[i,1] = r[i,1] * invfac
  end
  x[1,1] -= K1[1,2] * x[1,2] * invfac
  # right hand side of the rhow Schur complement: r3 - A31 q1 - A32 q2
  wv = MVector{Mm,FT}(undef)
  @unroll for i = 1 : Mm
    acc = r[i,3] - GravPW(Gr,Pi,Rm,Ps,ThA,DW,rJ,nS,i,iz,ID,Val(M),Val(VI)) * x[i,1]
    @unroll for j = 1 : Mm
      acc -= PresW1(Al,Pi,Rm,Ps,ThA,DW,rJ,nS,i,j,iz,ID,Val(M),Val(VI)) * x[j,2]
    end
    wv[i] = acc
  end
  wv[1] -= K1[3,2] * x[1,2]
  ldivFull!(iz,ID,SLU,wv,Val(Mm))
  @unroll for i = 1 : Mm
    x[i,3] = wv[i]
  end
  # rhotheta = q2 - A22^{-1} A23 w, rho = q1 - A11^{-1} (A12 (rhotheta - q2) + A13 w)
  @unroll for i = 1 : Mm
    acc = FT(0)
    @unroll for j = 1 : Mm
      acc += rJ[i] * DW[i,j] * nS[j] * Th[j] * x[j,3]
    end
    if i == 1
      acc += K1[2,3] * x[1,3]
      x[i,2] -= acc * d21
    else
      x[i,2] -= acc * invfac
    end
  end
  x[1,1] = (r[1,1] - K1[1,2] * x[1,2]) * invfac
  @unroll for i = 1 : Mm
    acc = FT(0)
    @unroll for j = 1 : Mm
      acc += MassOp1(Rm,DW,rJ,nS,i,j,iz,ID,Val(VI)) * x[j,3]
    end
    if i == 1
      acc += K1[1,3] * x[1,3]
    end
    x[i,1] -= acc * invfac
  end
end

# Cself x: rows (rho_M, rhotheta_M, rhow_M) applied to the interior unknowns
@inline function ApplyCself(x, Th, Al, rJ, nS, DW, Co, iz, ID, ::Val{M}, ::Val{VI}) where {M, VI}
  FT = eltype(x)
  ThA, Rm, Ps, Pi = Co
  c1 = FT(0)
  c2 = FT(0)
  c3 = FT(0)
  @unroll for j = 1 : M - 1
    d = rJ[M] * DW[M,j] * nS[j]
    c1 += MassOp1(Rm,DW,rJ,nS,M,j,iz,ID,Val(VI)) * x[j,3]
    c2 += d * Th[j] * x[j,3]
    c3 += PresW1(Al,Pi,Rm,Ps,ThA,DW,rJ,nS,M,j,iz,ID,Val(M),Val(VI)) * x[j,2]
  end
  return SVector{3,FT}(c1, c2, c3)
end

################################################################################
# Assembly: element factorization and the blocks of the reduced system
################################################################################

@kernel inbounds = true function FillJacOneSidedKernel!(SLU,Dself,Dprev,Lo,Up,
  @Const(Th),@Const(Al),@Const(Rm),@Const(Ps),@Const(Pi),@Const(Cs),
  @Const(Geo),@Const(J),@Const(Surf),@Const(DW),@Const(w),
  fac,cS,::Val{M},::Val{VI}) where {M,VI}

  iz,ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  ND = @uniform @ndrange()[2]

  if ID <= ND
    FT = eltype(Th)
    Mm = M - 1
    ThL, AlL, rJ, nS, Gr = LoadElementC(Th,Al,J,Surf,Geo,DW,iz,ID,Val(M))
    Co = (Th, Rm, Ps, Pi)
    invfac = FT(1) / fac
    sig1 = nS[1] * rJ[1] / w[1]
    sigM = nS[M] * rJ[M] / w[M]
    if iz > 1
      ThMm = Th[M,iz-1,ID]
      AlMm = Al[M,iz-1,ID]
    else
      ThMm = ThL[1]
      AlMm = AlL[1]
    end
    K1 = CornerBottom(ThL[1],AlL[1],ThMm,iz,sig1,cS,Th,Rm,Ps,Cs,ID,Val(M),Val(VI))

    # Schur complement for rhow:
    # S = A33 - A31 A11^{-1} A13 + A31 A11^{-1} A12 H - A32 H,  H = A22^{-1} A23
    d21 = FT(1) / (fac + K1[2,2])
    SL = MMatrix{Mm,Mm,FT}(undef)
    @unroll for i = 1 : Mm
      @unroll for j = 1 : Mm
        A13ij = MassOp1(Rm,DW,rJ,nS,i,j,iz,ID,Val(VI))
        if i == 1 && j == 1
          A13ij += K1[1,3]
        end
        Gri = GravPW(Gr,Pi,Rm,Ps,Th,DW,rJ,nS,i,iz,ID,Val(M),Val(VI))
        acc = -Gri * invfac * A13ij
        @unroll for k = 1 : Mm
          Hkj = rJ[k] * DW[k,j] * nS[j] * ThL[j]
          if k == 1 && j == 1
            Hkj += K1[2,3]
          end
          Hkj *= (k == 1 ? d21 : invfac)
          A32ik = PresW1(AlL,Pi,Rm,Ps,Th,DW,rJ,nS,i,k,iz,ID,Val(M),Val(VI))
          if i == 1 && k == 1
            A32ik += K1[3,2]
          end
          acc -= A32ik * Hkj
          if i == 1 && k == 1
            acc += Gri * invfac * K1[1,2] * Hkj
          end
        end
        SL[i,j] = acc
      end
      SL[i,i] += fac
    end
    SL[1,1] += K1[3,3]
    LUFull!(SL,Val(Mm))
    @unroll for i = 1 : Mm
      @unroll for j = 1 : Mm
        SLU[i,j,iz,ID] = SL[i,j]
      end
    end

    # Diagonal block D_e of the reduced unknowns
    D = MMatrix{3,3,FT}(undef)
    @unroll for i = 1 : 3
      @unroll for j = 1 : 3
        D[i,j] = FT(0)
      end
      D[i,i] = fac
    end
    dMM = rJ[M] * DW[M,M] * nS[M]
    D[3,1] += GravPW(Gr,Pi,Rm,Ps,Th,DW,rJ,nS,M,iz,ID,Val(M),Val(VI))
    D[1,3] += MassOp1(Rm,DW,rJ,nS,M,M,iz,ID,Val(VI))
    D[2,3] += dMM * ThL[M]
    D[3,2] += PresW1(AlL,Pi,Rm,Ps,Th,DW,rJ,nS,M,M,iz,ID,Val(M),Val(VI))
    Dtop = sigM * FaceTop(ThL[M],AlL[M],Th,iz,ID,nz,cS,Rm,Ps,Cs,Val(M),Val(VI))
    @unroll for i = 1 : 3
      @unroll for j = 1 : 3
        D[i,j] += Dtop[i,j]
      end
    end

    # Z = A_e^{-1} B for the columns rhotheta_M and rhow_M of Bself
    # (and of Bprev for iz > 1), the rho_M columns are zero
    r = MMatrix{Mm,3,FT}(undef)
    x = MMatrix{Mm,3,FT}(undef)
    CZs = MMatrix{3,3,FT}(undef)
    CZp = MMatrix{3,3,FT}(undef)
    NZs = MMatrix{3,3,FT}(undef)
    NZp = MMatrix{3,3,FT}(undef)
    @unroll for i = 1 : 3
      CZs[i,1] = FT(0)
      CZp[i,1] = FT(0)
      NZs[i,1] = FT(0)
      NZp[i,1] = FT(0)
    end
    if iz > 1
      sigMm = Surf[M,iz-1,ID] / (J[M,iz-1,ID] * w[M])
      Nm = sigMm * FaceBelow(2,1,Th,Al,ThL[1],AlL[1],iz,ID,cS,Rm,Ps,Cs,Val(M),Val(VI))
      P = -sig1 * FaceBelow(1,2,Th,Al,ThL[1],AlL[1],iz,ID,cS,Rm,Ps,Cs,Val(M),Val(VI))
    end
    @unroll for col = 2 : 3
      # Bself column
      @unroll for i = 1 : Mm
        d = rJ[i] * DW[i,M] * nS[M]
        if col == 2
          r[i,1] = FT(0)
          r[i,2] = FT(0)
          r[i,3] = PresW1(AlL,Pi,Rm,Ps,Th,DW,rJ,nS,i,M,iz,ID,Val(M),Val(VI))
        else
          r[i,1] = MassOp1(Rm,DW,rJ,nS,i,M,iz,ID,Val(VI))
          r[i,2] = d * ThL[M]
          r[i,3] = FT(0)
        end
      end
      SolveElement!(x, r, SLU, iz, ID, ThL, AlL, rJ, nS, Gr, K1, DW, fac, Co, Val(M), Val(VI))
      c = ApplyCself(x, ThL, AlL, rJ, nS, DW, Co, iz, ID, Val(M), Val(VI))
      @unroll for i = 1 : 3
        CZs[i,col] = c[i]
      end
      if iz > 1
        @unroll for i = 1 : 3
          NZs[i,col] = Nm[i,1] * x[1,1] + Nm[i,2] * x[1,2] + Nm[i,3] * x[1,3]
        end
        # Bprev column
        @unroll for i = 1 : Mm
          r[i,1] = FT(0)
          r[i,2] = FT(0)
          r[i,3] = FT(0)
        end
        r[1,1] = P[1,col]
        r[1,2] = P[2,col]
        r[1,3] = P[3,col]
        SolveElement!(x, r, SLU, iz, ID, ThL, AlL, rJ, nS, Gr, K1, DW, fac, Co, Val(M), Val(VI))
        c = ApplyCself(x, ThL, AlL, rJ, nS, DW, Co, iz, ID, Val(M), Val(VI))
        @unroll for i = 1 : 3
          CZp[i,col] = c[i]
          NZp[i,col] = Nm[i,1] * x[1,1] + Nm[i,2] * x[1,2] + Nm[i,3] * x[1,3]
        end
      end
    end

    @unroll for i = 1 : 3
      @unroll for j = 1 : 3
        Dself[ID,i,j,iz] = D[i,j] - CZs[i,j]
        if iz > 1
          Lo[ID,i,j,iz] = -CZp[i,j]
          Dprev[ID,i,j,iz] = -NZp[i,j]
          Up[ID,i,j,iz-1] = -NZs[i,j]
        end
      end
    end
  end
end

################################################################################
# Solve: forward elimination, reduced solve, back substitution
################################################################################

@kernel inbounds = true function ldivOneSidedFKernel!(rsA,rsB,@Const(v),@Const(SLU),
  @Const(Th),@Const(Al),@Const(Rm),@Const(Ps),@Const(Pi),@Const(Cs),
  @Const(Geo),@Const(J),@Const(Surf),@Const(NV),@Const(DW),@Const(w),
  fac,cS,::Val{M},::Val{VI}) where {M,VI}

  iz,ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  ND = @uniform @ndrange()[2]

  if ID <= ND
    FT = eltype(v)
    Mm = M - 1
    ThL, AlL, rJ, nS, Gr = LoadElementC(Th,Al,J,Surf,Geo,DW,iz,ID,Val(M))
    Co = (Th, Rm, Ps, Pi)
    sig1 = nS[1] * rJ[1] / w[1]
    ThMm = iz > 1 ? Th[M,iz-1,ID] : ThL[1]
    K1 = CornerBottom(ThL[1],AlL[1],ThMm,iz,sig1,cS,Th,Rm,Ps,Cs,ID,Val(M),Val(VI))
    r = MMatrix{Mm,3,FT}(undef)
    x = MMatrix{Mm,3,FT}(undef)
    @unroll for i = 1 : Mm
      r[i,1] = v[i,iz,ID,1]
      r[i,2] = v[i,iz,ID,5]
      r[i,3] = NuDot(NV,v,i,iz,ID)
    end
    SolveElement!(x, r, SLU, iz, ID, ThL, AlL, rJ, nS, Gr, K1, DW, fac, Co, Val(M), Val(VI))
    c = ApplyCself(x, ThL, AlL, rJ, nS, DW, Co, iz, ID, Val(M), Val(VI))
    rsA[ID,1,iz] = v[M,iz,ID,1] - c[1]
    rsA[ID,2,iz] = v[M,iz,ID,5] - c[2]
    rsA[ID,3,iz] = NuDot(NV,v,M,iz,ID) - c[3]
    if iz > 1
      sigMm = Surf[M,iz-1,ID] / (J[M,iz-1,ID] * w[M])
      Nm = sigMm * FaceBelow(2,1,Th,Al,ThL[1],AlL[1],iz,ID,cS,Rm,Ps,Cs,Val(M),Val(VI))
      @unroll for i = 1 : 3
        rsB[ID,i,iz] = -(Nm[i,1] * x[1,1] + Nm[i,2] * x[1,2] + Nm[i,3] * x[1,3])
      end
    end
  end
end

@kernel inbounds = true function ldivOneSidedBKernel!(k,@Const(v),@Const(rs),@Const(SLU),
  @Const(Th),@Const(Al),@Const(Rm),@Const(Ps),@Const(Pi),@Const(Cs),
  @Const(Geo),@Const(J),@Const(Surf),@Const(NV),@Const(DW),@Const(w),
  fac,facH,cS,::Val{M},::Val{VI}) where {M,VI}

  iz,ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  ND = @uniform @ndrange()[2]

  if ID <= ND
    FT = eltype(v)
    Mm = M - 1
    ThL, AlL, rJ, nS, Gr = LoadElementC(Th,Al,J,Surf,Geo,DW,iz,ID,Val(M))
    Co = (Th, Rm, Ps, Pi)
    sig1 = nS[1] * rJ[1] / w[1]
    ThMm = iz > 1 ? Th[M,iz-1,ID] : ThL[1]
    K1 = CornerBottom(ThL[1],AlL[1],ThMm,iz,sig1,cS,Th,Rm,Ps,Cs,ID,Val(M),Val(VI))
    y1 = rs[ID,1,iz]
    y2 = rs[ID,2,iz]
    y3 = rs[ID,3,iz]
    r = MMatrix{Mm,3,FT}(undef)
    x = MMatrix{Mm,3,FT}(undef)
    # r = f - Bself y_e - Bprev y_{e-1}
    @unroll for i = 1 : Mm
      d = rJ[i] * DW[i,M] * nS[M]
      r[i,1] = v[i,iz,ID,1] - MassOp1(Rm,DW,rJ,nS,i,M,iz,ID,Val(VI)) * y3
      r[i,2] = v[i,iz,ID,5] - d * ThL[M] * y3
      r[i,3] = NuDot(NV,v,i,iz,ID) - PresW1(AlL,Pi,Rm,Ps,Th,DW,rJ,nS,i,M,iz,ID,Val(M),Val(VI)) * y2
    end
    if iz > 1
      P = -sig1 * FaceBelow(1,2,Th,Al,ThL[1],AlL[1],iz,ID,cS,Rm,Ps,Cs,Val(M),Val(VI))
      ym2 = rs[ID,2,iz-1]
      ym3 = rs[ID,3,iz-1]
      @unroll for i = 1 : 3
        r[1,i] -= P[i,2] * ym2 + P[i,3] * ym3
      end
    end
    SolveElement!(x, r, SLU, iz, ID, ThL, AlL, rJ, nS, Gr, K1, DW, fac, Co, Val(M), Val(VI))
    @unroll for i = 1 : Mm
      StoreSolution!(k,v,NV,i,iz,ID,x[i,1],x[i,2],x[i,3],facH)
    end
    StoreSolution!(k,v,NV,M,iz,ID,y1,y2,y3,facH)
  end
end

################################################################################
# Driver routines
################################################################################

function FillJacOneSided!(Jac::JacOneSided,DG,Metric,fac,Phys)
  backend = get_backend(Jac.Th)
  M = Jac.M
  nz = Jac.nz
  ND = DG.NumI
  Jac.fac = fac
  Jac.cS = Phys.cS
  NDG = max(1,min(div(256,nz),ND))
  KFillJacOneSidedKernel! = FillJacOneSidedKernel!(backend,(nz,NDG))
  KFillJacOneSidedKernel!(Jac.SLU,Jac.Dself,Jac.Dprev,Jac.Lo,Jac.Up,
    Jac.Th,Jac.Al,Jac.Rm,Jac.Ps,Jac.Pi,Jac.Cs,Jac.Geo,Metric.J,Metric.VolSurfVV,DG.V.DWZ,DG.V.wZ,
    Jac.fac,Jac.cS,Val(M),Val(Jac.VI);ndrange=(nz,ND))
  KFactorBlockTriKernel! = FactorBlockTriKernel!(backend,(32,))
  KFactorBlockTriKernel!(Jac.Binv,Jac.Cp,Jac.Dself,Jac.Dprev,Jac.Lo,Jac.Up,Val(3),Val(nz),Val(true);
    ndrange=(ND,))
end

# k = (fac I - J)^{-1} v, the horizontal momentum is k_h = facH v_h
function SolveOneSided!(Jac::JacOneSided,k,v,DG,Metric,facH)
  backend = get_backend(v)
  M = Jac.M
  nz = Jac.nz
  ND = DG.NumI
  NDG = max(1,min(div(256,nz),ND))
  KldivOneSidedFKernel! = ldivOneSidedFKernel!(backend,(nz,NDG))
  KldivOneSidedFKernel!(Jac.rsA,Jac.rsB,v,Jac.SLU,Jac.Th,Jac.Al,Jac.Rm,Jac.Ps,Jac.Pi,Jac.Cs,Jac.Geo,Metric.J,
    Metric.VolSurfVV,Metric.NVV,DG.V.DWZ,DG.V.wZ,Jac.fac,Jac.cS,Val(M),Val(Jac.VI);ndrange=(nz,ND))
  KldivBlockTriKernel! = ldivBlockTriKernel!(backend,(32,))
  KldivBlockTriKernel!(Jac.rs,Jac.Binv,Jac.Cp,Jac.Lo,Jac.rsA,Jac.rsB,Val(3),Val(nz),Val(true);
    ndrange=(ND,))
  KldivOneSidedBKernel! = ldivOneSidedBKernel!(backend,(nz,NDG))
  KldivOneSidedBKernel!(k,v,Jac.rs,Jac.SLU,Jac.Th,Jac.Al,Jac.Rm,Jac.Ps,Jac.Pi,Jac.Cs,Jac.Geo,Metric.J,
    Metric.VolSurfVV,Metric.NVV,DG.V.DWZ,DG.V.wZ,Jac.fac,facH,Jac.cS,Val(M),Val(Jac.VI);ndrange=(nz,ND))
end

function Jac!(U,fac,DG,Metric,Phys,Cache,JCache::JacOneSided,Global)
  InitGeo!(JCache,Cache,DG,Global)
  if Global.Model.ModelType == "VectorInvariant"
    JCache.VI = true
    FillCoeffVI!(JCache.Th,JCache.Al,JCache.Rm,JCache.Ps,JCache.Pi,JCache.Cs,U,Phys,Global.Model.SoundSpeedDG)
  else
    FillThetaAlpha!(JCache.Th,JCache.Al,U,Phys)
  end
  FillJacOneSided!(JCache,DG,Metric,eltype(U)(1) / fac,Phys)
end

function Solve!(k,v,Jac::JacOneSided,fac,DG::FiniteElements.DGElement,Metric,Global)
  SolveOneSided!(Jac,k,v,DG,Metric,fac)
end
