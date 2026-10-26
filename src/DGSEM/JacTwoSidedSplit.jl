################################################################################
# Vertical HEVI solver with two-sided interface condensation and split-form
# gravity (arXiv:2610.15775, Sec. 4.2.1, gravity as in JacHDG)
#
# Linear system per column: (fac I - J) x = b with the HEVI Jacobian of
# JacOneSided (contravariant vertical direction, see JacCommon.jl), but the
# gravity term is the split form (n_i / J_i) G_ik rho_k, G from the geopotential.
#
# Per element e the unknowns are split into
#   y_e : rhotheta_1, rhow_1, rhotheta_M, rhow_M          (reduced unknowns)
#   x_e : rho_1..M, rhotheta_2..M-1, rhow_2..M-1          (eliminated element-wise)
# The interface fluxes only depend on rhotheta and rhow, therefore the rows of
# y_e couple to x_e, y_{e-1}, y_e, y_{e+1} and the rows of x_e to y_{e-1}, y_e,
# y_{e+1} (rho_1 and rho_M). Eliminating x_e gives a block tridiagonal system
# for y with 4x4 blocks, each element computes its own block row.
#
# Element solve A_e x = r: rho and rhotheta are eliminated analytically, the
# Schur complement for rhow_2..M-1 is factored with LUFull! (stored in SLU).
# Requires M >= 3.
#
# Vector invariant form (VI = true, variables rho, rhotheta, w): the same
# elimination with the node coefficients of JacHDG (Rm, Th = rhotheta,
# Al = dPi/drhotheta, Ps = cp theta, Pi, Cs). The nonconservative pressure term
# is written as weak form plus the end node terms, which together with
# 1/2 H [[g]] give the face pressure Ps_own Pi_own -+ 1/2 Ps_face [[Pi]].
#
# Memory layouts
#   Geo, Th, Al, Rm, Ps, Pi, Cs     : (M, nz, NumI)
#   SLU                             : (M-2, M-2, nz, NumI)
#   Dself, Lo, Up, Binv, Cp         : (NumI, 4, 4, nz)
#   rsA, rs                         : (NumI, 4, nz)
################################################################################

mutable struct JacTwoSidedSplit{FT<:AbstractFloat,
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
  Lo::AT4
  Up::AT4
  Binv::AT4
  Cp::AT4
  rsA::AT3
  rs::AT3
end

function JacTwoSidedSplit(backend,FT,M,nz,DG)
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
  SLU = KernelAbstractions.zeros(backend,FT,M-2,M-2,nz,NumI)
  Dself = KernelAbstractions.zeros(backend,FT,NumI,4,4,nz)
  Lo = KernelAbstractions.zeros(backend,FT,NumI,4,4,nz)
  Up = KernelAbstractions.zeros(backend,FT,NumI,4,4,nz)
  Binv = KernelAbstractions.zeros(backend,FT,NumI,4,4,nz)
  Cp = KernelAbstractions.zeros(backend,FT,NumI,4,4,nz)
  rsA = KernelAbstractions.zeros(backend,FT,NumI,4,nz)
  rs = KernelAbstractions.zeros(backend,FT,NumI,4,nz)
  return JacTwoSidedSplit{FT,
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
    Lo,
    Up,
    Binv,
    Cp,
    rsA,
    rs,
  )
end

################################################################################
# Element data and element operators
################################################################################

# Split form gravity rows (n_i / J_i) G_ik for all nodes
@inline function GravSplitMatrix(Geo,DW,rJ,nS,iz,ID,::Val{M}) where {M}
  FT = eltype(Geo)
  Gm = MMatrix{M,M,FT}(undef)
  @unroll for i = 1 : M
    @unroll for k = 1 : M
      Gm[i,k] = GravSplit(Geo,DW,rJ,nS,i,k,iz,ID,Val(M))
    end
  end
  return Gm
end

# Vector invariant split buoyancy Bs_ik (see build_GP in JacHDG.jl)
@inline function BuoySplit(Pi,Rm,Ps,Th,DW,i,k,iz,ID,::Val{M}) where {M}
  FT = eltype(Pi)
  Pii = Pi[i,iz,ID]
  val = FT(0.5) * DW[i,k] * (Pi[k,iz,ID] - Pii)
  if i == k
    acc = FT(0)
    @unroll for m = 1 : M
      acc += DW[i,m] * (Pi[m,iz,ID] - Pii)
    end
    val += FT(0.5) * acc
  end
  return Ps[i,iz,ID] * Rm[i,iz,ID] / Th[i,iz,ID] * val
end

# Gravity matrix (acts on rho), conservative: split form gravity
@inline GravMatrix(Geo,Pi,Rm,Ps,Th,DW,rJ,nS,iz,ID,::Val{M},::Val{false}) where {M} =
  GravSplitMatrix(Geo,DW,rJ,nS,iz,ID,Val(M))

@inline function GravMatrix(Geo,Pi,Rm,Ps,Th,DW,rJ,nS,iz,ID,::Val{M},::Val{true}) where {M}
  FT = eltype(Geo)
  Gm = MMatrix{M,M,FT}(undef)
  @unroll for i = 1 : M
    @unroll for k = 1 : M
      rRho = FT(1) / Rm[k,iz,ID]
      Gm[i,k] = -rJ[i] * nS[i] * BuoySplit(Pi,Rm,Ps,Th,DW,i,k,iz,ID,Val(M)) *
        Th[k,iz,ID] * rRho * rRho
    end
  end
  return Gm
end

# Pressure coupling of the vertical momentum row i to rhotheta_j
@inline PresW(Al,Pi,Rm,Ps,Th,DW,rJ,nS,i,j,iz,ID,::Val{M},::Val{false}) where {M} =
  rJ[i] * DW[i,j] * nS[j] * Al[j]
@inline PresW(Al,Pi,Rm,Ps,Th,DW,rJ,nS,i,j,iz,ID,::Val{M},::Val{true}) where {M} =
  rJ[i] * nS[i] * (Ps[i,iz,ID] * DW[i,j] * Al[j] +
    BuoySplit(Pi,Rm,Ps,Th,DW,i,j,iz,ID,Val(M)) / Rm[j,iz,ID])

# x * (mass flux operator)_kj
@inline MulMass(x,Rm,DW,rJ,nS,k,j,iz,ID,::Val{false}) = x * rJ[k] * DW[k,j] * nS[j]
@inline MulMass(x,Rm,DW,rJ,nS,k,j,iz,ID,::Val{true}) = x * rJ[k] * DW[k,j] * nS[j] * Rm[j,iz,ID]
# (mass flux operator)_kj
@inline MassOp(Rm,DW,rJ,nS,k,j,iz,ID,::Val{false}) = rJ[k] * DW[k,j] * nS[j]
@inline MassOp(Rm,DW,rJ,nS,k,j,iz,ID,::Val{true}) = rJ[k] * DW[k,j] * nS[j] * Rm[j,iz,ID]

# Interface blocks of element e: K1 (lower surface, own node 1),
# KM (upper surface, own node M), P (node 1 rows, node M of e-1),
# Q (node M rows, node 1 of e+1)
@inline function InterfaceBlocks(Th,Al,ThL,AlL,sig1,sigM,iz,ID,nz,cS,::Val{M}) where {M}
  if iz > 1
    ThMm = Th[M,iz-1,ID]
    AlMm = Al[M,iz-1,ID]
    K1 = -sig1 * LMARSdR(ThMm,ThL[1],AlL[1],cS)
    P = -sig1 * LMARSdL(ThMm,ThL[1],AlMm,cS)
  else
    K1 = -sig1 * WallFlux(AlL[1],cS,-1)
    P = zero(K1)
  end
  if iz < nz
    ThP = Th[1,iz+1,ID]
    AlP = Al[1,iz+1,ID]
    KM = sigM * LMARSdL(ThL[M],ThP,AlL[M],cS)
    Q = sigM * LMARSdR(ThL[M],ThP,AlP,cS)
  else
    KM = sigM * WallFlux(AlL[M],cS,1)
    Q = zero(KM)
  end
  return K1, KM, P, Q
end

@inline InterfaceBlocks(Th,Al,Rm,Ps,Cs,ThL,AlL,sig1,sigM,iz,ID,nz,cS,::Val{M},::Val{false}) where {M} =
  InterfaceBlocks(Th,Al,ThL,AlL,sig1,sigM,iz,ID,nz,cS,Val(M))

@inline function InterfaceBlocks(Th,Al,Rm,Ps,Cs,ThL,AlL,sig1,sigM,iz,ID,nz,cS,::Val{M},::Val{true}) where {M}
  FT = eltype(Th)
  if iz > 1
    c = FT(0.5) * (Cs[M,iz-1,ID] + Cs[1,iz,ID])
    a = (Th[M,iz-1,ID],ThL[1],Rm[M,iz-1,ID],Rm[1,iz,ID])
    K1 = -sig1 * VIdR(a...,AlL[1],Ps[M,iz-1,ID],Ps[1,iz,ID],c,2)
    P = -sig1 * VIdL(a...,Al[M,iz-1,ID],Ps[M,iz-1,ID],Ps[1,iz,ID],c,2)
  else
    K1 = -sig1 * WallFlux(Ps[1,iz,ID] * AlL[1],Cs[1,iz,ID],-1)
    P = zero(K1)
  end
  if iz < nz
    c = FT(0.5) * (Cs[M,iz,ID] + Cs[1,iz+1,ID])
    a = (ThL[M],Th[1,iz+1,ID],Rm[M,iz,ID],Rm[1,iz+1,ID])
    KM = sigM * VIdL(a...,AlL[M],Ps[M,iz,ID],Ps[1,iz+1,ID],c,1)
    Q = sigM * VIdR(a...,Al[1,iz+1,ID],Ps[M,iz,ID],Ps[1,iz+1,ID],c,1)
  else
    KM = sigM * WallFlux(Ps[M,iz,ID] * AlL[M],Cs[M,iz,ID],1)
    Q = zero(KM)
  end
  return K1, KM, P, Q
end

# Element solve A_e x = r with x = (xr[1:M], xt[1:M-2], xw[1:M-2])
@inline function SolveElementTS!(xr, xt, xw, rr, rt, rw, SLU, iz, ID, Th, Al, rJ, nS, Gm,
  DW, fac, Co, ::Val{M}, ::Val{VI}) where {M, VI}
  ThA, Rm, Ps, Pi = Co
  FT = eltype(xr)
  Mi = M - 2
  invfac = FT(1) / fac
  @unroll for k = 1 : M
    xr[k] = rr[k] * invfac
  end
  @unroll for a = 1 : Mi
    xt[a] = rt[a] * invfac
  end
  @unroll for a = 1 : Mi
    i = a + 1
    acc = rw[a]
    @unroll for k = 1 : M
      acc -= Gm[i,k] * xr[k]
    end
    @unroll for b = 1 : Mi
      j = b + 1
      acc -= PresW(Al,Pi,Rm,Ps,ThA,DW,rJ,nS,i,j,iz,ID,Val(M),Val(VI)) * xt[b]
    end
    xw[a] = acc
  end
  ldivFull!(iz,ID,SLU,xw,Val(Mi))
  @unroll for a = 1 : Mi
    i = a + 1
    acc = FT(0)
    @unroll for b = 1 : Mi
      j = b + 1
      acc += rJ[i] * DW[i,j] * nS[j] * Th[j] * xw[b]
    end
    xt[a] -= acc * invfac
  end
  @unroll for k = 1 : M
    acc = FT(0)
    @unroll for b = 1 : Mi
      j = b + 1
      acc += MassOp(Rm,DW,rJ,nS,k,j,iz,ID,Val(VI)) * xw[b]
    end
    xr[k] -= acc * invfac
  end
end

# Cself x: rows (rhotheta_1, rhow_1, rhotheta_M, rhow_M)
@inline function ApplyCselfTS(xr, xt, xw, Th, Al, rJ, nS, Gm, DW, Co, iz, ID, ::Val{M}, ::Val{VI}) where {M, VI}
  FT = eltype(xr)
  ThA, Rm, Ps, Pi = Co
  c1 = FT(0)
  c2 = FT(0)
  c3 = FT(0)
  c4 = FT(0)
  @unroll for k = 1 : M
    c2 += Gm[1,k] * xr[k]
    c4 += Gm[M,k] * xr[k]
  end
  @unroll for b = 1 : M - 2
    j = b + 1
    d1 = rJ[1] * DW[1,j] * nS[j]
    dM = rJ[M] * DW[M,j] * nS[j]
    c1 += d1 * Th[j] * xw[b]
    c2 += PresW(Al,Pi,Rm,Ps,ThA,DW,rJ,nS,1,j,iz,ID,Val(M),Val(VI)) * xt[b]
    c3 += dM * Th[j] * xw[b]
    c4 += PresW(Al,Pi,Rm,Ps,ThA,DW,rJ,nS,M,j,iz,ID,Val(M),Val(VI)) * xt[b]
  end
  return SVector{4,FT}(c1, c2, c3, c4)
end

# (rr, rt, rw) = Bself y for y = (rhotheta_1, rhow_1, rhotheta_M, rhow_M)
@inline function ApplyBselfTS!(rr, rt, rw, y, Th, Al, rJ, nS, K1, KM, DW, Co, iz, ID, ::Val{M}, ::Val{VI}) where {M, VI}
  ThA, Rm, Ps, Pi = Co
  @unroll for k = 1 : M
    if VI
      rr[k] = MassOp(Rm,DW,rJ,nS,k,1,iz,ID,Val(VI)) * y[2] + MassOp(Rm,DW,rJ,nS,k,M,iz,ID,Val(VI)) * y[4]
    else
      rr[k] = rJ[k] * (DW[k,1] * nS[1] * y[2] + DW[k,M] * nS[M] * y[4])
    end
  end
  rr[1] += K1[1,2] * y[1] + K1[1,3] * y[2]
  rr[M] += KM[1,2] * y[3] + KM[1,3] * y[4]
  @unroll for a = 1 : M - 2
    i = a + 1
    rt[a] = rJ[i] * (DW[i,1] * nS[1] * Th[1] * y[2] + DW[i,M] * nS[M] * Th[M] * y[4])
    if VI
      rw[a] = PresW(Al,Pi,Rm,Ps,ThA,DW,rJ,nS,i,1,iz,ID,Val(M),Val(VI)) * y[1] +
        PresW(Al,Pi,Rm,Ps,ThA,DW,rJ,nS,i,M,iz,ID,Val(M),Val(VI)) * y[3]
    else
      rw[a] = rJ[i] * (DW[i,1] * nS[1] * Al[1] * y[1] + DW[i,M] * nS[M] * Al[M] * y[3])
    end
  end
end

@inline UnitVec4(FT,j) = SVector{4,FT}(j == 1, j == 2, j == 3, j == 4)

################################################################################
# Assembly
################################################################################

@kernel inbounds = true function FillJacTwoSidedSplitKernel!(SLU,Dself,Lo,Up,
  @Const(Th),@Const(Al),@Const(Rm),@Const(Ps),@Const(Pi),@Const(Cs),
  @Const(Geo),@Const(J),@Const(Surf),@Const(DW),@Const(w),
  fac,cS,::Val{M},::Val{VI}) where {M,VI}

  iz,ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  ND = @uniform @ndrange()[2]

  if ID <= ND
    FT = eltype(Th)
    Mi = M - 2
    ThL, AlL, rJ, nS, _ = LoadElementC(Th,Al,J,Surf,Geo,DW,iz,ID,Val(M))
    Co = (Th, Rm, Ps, Pi)
    Gm = GravMatrix(Geo,Pi,Rm,Ps,Th,DW,rJ,nS,iz,ID,Val(M),Val(VI))
    invfac = FT(1) / fac
    sig1 = nS[1] * rJ[1] / w[1]
    sigM = nS[M] * rJ[M] / w[M]
    K1, KM, P, Q = InterfaceBlocks(Th,Al,Rm,Ps,Cs,ThL,AlL,sig1,sigM,iz,ID,nz,cS,Val(M),Val(VI))

    # Schur complement for rhow_2..M-1: S = fac I - (A31 A13 + A32 A23) / fac
    SL = MMatrix{Mi,Mi,FT}(undef)
    @unroll for a = 1 : Mi
      i = a + 1
      @unroll for b = 1 : Mi
        j = b + 1
        acc = FT(0)
        @unroll for k = 1 : M
          acc += MulMass(Gm[i,k],Rm,DW,rJ,nS,k,j,iz,ID,Val(VI))
        end
        @unroll for c = 1 : Mi
          l = c + 1
          acc += PresW(AlL,Pi,Rm,Ps,Th,DW,rJ,nS,i,l,iz,ID,Val(M),Val(VI)) * rJ[l] * DW[l,j] * nS[j] * ThL[j]
        end
        SL[a,b] = -acc * invfac
      end
      SL[a,a] += fac
    end
    LUFull!(SL,Val(Mi))
    @unroll for a = 1 : Mi
      @unroll for b = 1 : Mi
        SLU[a,b,iz,ID] = SL[a,b]
      end
    end

    # Coupling of the reduced unknowns among themselves
    D = MMatrix{4,4,FT}(undef)
    L = MMatrix{4,4,FT}(undef)
    U = MMatrix{4,4,FT}(undef)
    @unroll for i = 1 : 4
      @unroll for j = 1 : 4
        D[i,j] = FT(0)
        L[i,j] = FT(0)
        U[i,j] = FT(0)
      end
      D[i,i] = fac
    end
    d11 = rJ[1] * DW[1,1] * nS[1]
    d1M = rJ[1] * DW[1,M] * nS[M]
    dM1 = rJ[M] * DW[M,1] * nS[1]
    dMM = rJ[M] * DW[M,M] * nS[M]
    D[1,2] += d11 * ThL[1]
    D[1,4] += d1M * ThL[M]
    D[2,1] += PresW(AlL,Pi,Rm,Ps,Th,DW,rJ,nS,1,1,iz,ID,Val(M),Val(VI))
    D[2,3] += PresW(AlL,Pi,Rm,Ps,Th,DW,rJ,nS,1,M,iz,ID,Val(M),Val(VI))
    D[3,2] += dM1 * ThL[1]
    D[3,4] += dMM * ThL[M]
    D[4,1] += PresW(AlL,Pi,Rm,Ps,Th,DW,rJ,nS,M,1,iz,ID,Val(M),Val(VI))
    D[4,3] += PresW(AlL,Pi,Rm,Ps,Th,DW,rJ,nS,M,M,iz,ID,Val(M),Val(VI))
    @unroll for i = 1 : 2
      @unroll for j = 1 : 2
        D[i,j] += K1[i+1,j+1]
        D[i+2,j+2] += KM[i+1,j+1]
        L[i,j+2] = P[i+1,j+1]
        U[i+2,j] = Q[i+1,j+1]
      end
    end

    # Elimination of x_e: S_ef = D_ef - Cself A_e^{-1} B_ef
    rr = MVector{M,FT}(undef)
    rt = MVector{Mi,FT}(undef)
    rw = MVector{Mi,FT}(undef)
    xr = MVector{M,FT}(undef)
    xt = MVector{Mi,FT}(undef)
    xw = MVector{Mi,FT}(undef)
    @unroll for col = 1 : 4
      ApplyBselfTS!(rr, rt, rw, UnitVec4(FT,col), ThL, AlL, rJ, nS, K1, KM, DW, Co, iz, ID, Val(M), Val(VI))
      SolveElementTS!(xr, xt, xw, rr, rt, rw, SLU, iz, ID, ThL, AlL, rJ, nS, Gm, DW, fac, Co, Val(M), Val(VI))
      c = ApplyCselfTS(xr, xt, xw, ThL, AlL, rJ, nS, Gm, DW, Co, iz, ID, Val(M), Val(VI))
      @unroll for i = 1 : 4
        D[i,col] -= c[i]
      end
    end
    # Bprev: rho_1 row, columns rhotheta_M and rhow_M of element e-1
    # Bnext: rho_M row, columns rhotheta_1 and rhow_1 of element e+1
    @unroll for col = 1 : 4
      if (col >= 3 && iz > 1) || (col <= 2 && iz < nz)
        @unroll for k = 1 : M
          rr[k] = FT(0)
        end
        @unroll for a = 1 : Mi
          rt[a] = FT(0)
          rw[a] = FT(0)
        end
        if col >= 3
          rr[1] = P[1,col-1]
        else
          rr[M] = Q[1,col+1]
        end
        SolveElementTS!(xr, xt, xw, rr, rt, rw, SLU, iz, ID, ThL, AlL, rJ, nS, Gm, DW, fac, Co, Val(M), Val(VI))
        c = ApplyCselfTS(xr, xt, xw, ThL, AlL, rJ, nS, Gm, DW, Co, iz, ID, Val(M), Val(VI))
        @unroll for i = 1 : 4
          if col >= 3
            L[i,col] -= c[i]
          else
            U[i,col] -= c[i]
          end
        end
      end
    end

    @unroll for i = 1 : 4
      @unroll for j = 1 : 4
        Dself[ID,i,j,iz] = D[i,j]
        Lo[ID,i,j,iz] = L[i,j]
        Up[ID,i,j,iz] = U[i,j]
      end
    end
  end
end

################################################################################
# Solve: forward elimination, reduced solve, back substitution
################################################################################

@kernel inbounds = true function ldivTwoSidedSplitFKernel!(rsA,@Const(v),@Const(SLU),
  @Const(Th),@Const(Al),@Const(Rm),@Const(Ps),@Const(Pi),
  @Const(Geo),@Const(J),@Const(Surf),@Const(NV),@Const(DW),
  fac,::Val{M},::Val{VI}) where {M,VI}

  iz,ID = @index(Global, NTuple)

  ND = @uniform @ndrange()[2]

  if ID <= ND
    FT = eltype(v)
    Mi = M - 2
    ThL, AlL, rJ, nS, _ = LoadElementC(Th,Al,J,Surf,Geo,DW,iz,ID,Val(M))
    Co = (Th, Rm, Ps, Pi)
    Gm = GravMatrix(Geo,Pi,Rm,Ps,Th,DW,rJ,nS,iz,ID,Val(M),Val(VI))
    rr = MVector{M,FT}(undef)
    rt = MVector{Mi,FT}(undef)
    rw = MVector{Mi,FT}(undef)
    xr = MVector{M,FT}(undef)
    xt = MVector{Mi,FT}(undef)
    xw = MVector{Mi,FT}(undef)
    @unroll for k = 1 : M
      rr[k] = v[k,iz,ID,1]
    end
    @unroll for a = 1 : Mi
      rt[a] = v[a+1,iz,ID,5]
      rw[a] = NuDot(NV,v,a+1,iz,ID)
    end
    SolveElementTS!(xr, xt, xw, rr, rt, rw, SLU, iz, ID, ThL, AlL, rJ, nS, Gm, DW, fac, Co, Val(M), Val(VI))
    c = ApplyCselfTS(xr, xt, xw, ThL, AlL, rJ, nS, Gm, DW, Co, iz, ID, Val(M), Val(VI))
    rsA[ID,1,iz] = v[1,iz,ID,5] - c[1]
    rsA[ID,2,iz] = NuDot(NV,v,1,iz,ID) - c[2]
    rsA[ID,3,iz] = v[M,iz,ID,5] - c[3]
    rsA[ID,4,iz] = NuDot(NV,v,M,iz,ID) - c[4]
  end
end

@kernel inbounds = true function ldivTwoSidedSplitBKernel!(k,@Const(v),@Const(rs),@Const(SLU),
  @Const(Th),@Const(Al),@Const(Rm),@Const(Ps),@Const(Pi),@Const(Cs),
  @Const(Geo),@Const(J),@Const(Surf),@Const(NV),@Const(DW),@Const(w),
  fac,facH,cS,::Val{M},::Val{VI}) where {M,VI}

  iz,ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  ND = @uniform @ndrange()[2]

  if ID <= ND
    FT = eltype(v)
    Mi = M - 2
    ThL, AlL, rJ, nS, _ = LoadElementC(Th,Al,J,Surf,Geo,DW,iz,ID,Val(M))
    Co = (Th, Rm, Ps, Pi)
    Gm = GravMatrix(Geo,Pi,Rm,Ps,Th,DW,rJ,nS,iz,ID,Val(M),Val(VI))
    sig1 = nS[1] * rJ[1] / w[1]
    sigM = nS[M] * rJ[M] / w[M]
    K1, KM, P, Q = InterfaceBlocks(Th,Al,Rm,Ps,Cs,ThL,AlL,sig1,sigM,iz,ID,nz,cS,Val(M),Val(VI))
    y = SVector{4,FT}(rs[ID,1,iz], rs[ID,2,iz], rs[ID,3,iz], rs[ID,4,iz])
    rr = MVector{M,FT}(undef)
    rt = MVector{Mi,FT}(undef)
    rw = MVector{Mi,FT}(undef)
    xr = MVector{M,FT}(undef)
    xt = MVector{Mi,FT}(undef)
    xw = MVector{Mi,FT}(undef)
    # r = f - Bself y_e - Bprev y_{e-1} - Bnext y_{e+1}
    ApplyBselfTS!(rr, rt, rw, y, ThL, AlL, rJ, nS, K1, KM, DW, Co, iz, ID, Val(M), Val(VI))
    @unroll for k = 1 : M
      rr[k] = v[k,iz,ID,1] - rr[k]
    end
    @unroll for a = 1 : Mi
      rt[a] = v[a+1,iz,ID,5] - rt[a]
      rw[a] = NuDot(NV,v,a+1,iz,ID) - rw[a]
    end
    if iz > 1
      rr[1] -= P[1,2] * rs[ID,3,iz-1] + P[1,3] * rs[ID,4,iz-1]
    end
    if iz < nz
      rr[M] -= Q[1,2] * rs[ID,1,iz+1] + Q[1,3] * rs[ID,2,iz+1]
    end
    SolveElementTS!(xr, xt, xw, rr, rt, rw, SLU, iz, ID, ThL, AlL, rJ, nS, Gm, DW, fac, Co, Val(M), Val(VI))
    StoreSolution!(k,v,NV,1,iz,ID,xr[1],y[1],y[2],facH)
    @unroll for a = 1 : Mi
      StoreSolution!(k,v,NV,a+1,iz,ID,xr[a+1],xt[a],xw[a],facH)
    end
    StoreSolution!(k,v,NV,M,iz,ID,xr[M],y[3],y[4],facH)
  end
end

################################################################################
# Driver routines
################################################################################

function FillJacTwoSidedSplit!(Jac::JacTwoSidedSplit,DG,Metric,fac,Phys)
  backend = get_backend(Jac.Th)
  M = Jac.M
  nz = Jac.nz
  ND = DG.NumI
  Jac.fac = fac
  Jac.cS = Phys.cS
  NDG = max(1,min(div(256,nz),ND))
  KFillJacTwoSidedSplitKernel! = FillJacTwoSidedSplitKernel!(backend,(nz,NDG))
  KFillJacTwoSidedSplitKernel!(Jac.SLU,Jac.Dself,Jac.Lo,Jac.Up,
    Jac.Th,Jac.Al,Jac.Rm,Jac.Ps,Jac.Pi,Jac.Cs,Jac.Geo,Metric.J,Metric.VolSurfVV,DG.V.DWZ,DG.V.wZ,
    Jac.fac,Jac.cS,Val(M),Val(Jac.VI);ndrange=(nz,ND))
  KFactorBlockTriKernel! = FactorBlockTriKernel!(backend,(32,))
  KFactorBlockTriKernel!(Jac.Binv,Jac.Cp,Jac.Dself,Jac.Dself,Jac.Lo,Jac.Up,Val(4),Val(nz),Val(false);
    ndrange=(ND,))
end

# k = (fac I - J)^{-1} v, the horizontal momentum is k_h = facH v_h
function SolveTwoSidedSplit!(Jac::JacTwoSidedSplit,k,v,DG,Metric,facH)
  backend = get_backend(v)
  M = Jac.M
  nz = Jac.nz
  ND = DG.NumI
  NDG = max(1,min(div(256,nz),ND))
  KldivTwoSidedSplitFKernel! = ldivTwoSidedSplitFKernel!(backend,(nz,NDG))
  KldivTwoSidedSplitFKernel!(Jac.rsA,v,Jac.SLU,Jac.Th,Jac.Al,Jac.Rm,Jac.Ps,Jac.Pi,Jac.Geo,Metric.J,
    Metric.VolSurfVV,Metric.NVV,DG.V.DWZ,Jac.fac,Val(M),Val(Jac.VI);ndrange=(nz,ND))
  KldivBlockTriKernel! = ldivBlockTriKernel!(backend,(32,))
  KldivBlockTriKernel!(Jac.rs,Jac.Binv,Jac.Cp,Jac.Lo,Jac.rsA,Jac.rsA,Val(4),Val(nz),Val(false);
    ndrange=(ND,))
  KldivTwoSidedSplitBKernel! = ldivTwoSidedSplitBKernel!(backend,(nz,NDG))
  KldivTwoSidedSplitBKernel!(k,v,Jac.rs,Jac.SLU,Jac.Th,Jac.Al,Jac.Rm,Jac.Ps,Jac.Pi,Jac.Cs,Jac.Geo,Metric.J,
    Metric.VolSurfVV,Metric.NVV,DG.V.DWZ,DG.V.wZ,Jac.fac,facH,Jac.cS,Val(M),Val(Jac.VI);ndrange=(nz,ND))
end

function Jac!(U,fac,DG,Metric,Phys,Cache,JCache::JacTwoSidedSplit,Global)
  InitGeo!(JCache,Cache,DG,Global)
  if Global.Model.ModelType == "VectorInvariant"
    JCache.VI = true
    FillCoeffVI!(JCache.Th,JCache.Al,JCache.Rm,JCache.Ps,JCache.Pi,JCache.Cs,U,Phys)
  else
    FillThetaAlpha!(JCache.Th,JCache.Al,U,Phys)
  end
  FillJacTwoSidedSplit!(JCache,DG,Metric,eltype(U)(1) / fac,Phys)
end

function Solve!(k,v,Jac::JacTwoSidedSplit,fac,DG::FiniteElements.DGElement,Metric,Global)
  SolveTwoSidedSplit!(Jac,k,v,DG,Metric,fac)
end
