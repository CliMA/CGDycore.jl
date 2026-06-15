# ============================================================
# Refactored AverageFlux functions — no @views, no slices
#
# KEY CHANGE: instead of passing pre-sliced vectors
#   FluxNonLinAver!(flux, VL, VR, AuxL, AuxR, m_L, m_R)
# we now pass the parent shared-memory arrays + indices:
#   FluxNonLinAver!(flux, VLoc, AuxLoc, dXdxILoc,
#                   K, Iz, iD,   <- left state
#                   l, Iz, iD,   <- right state (vertical kernel example)
#                   Val(NV))
#
# All element accesses become direct array indexing,
# which the GPU compiler can fully inline and unroll.
# ============================================================

abstract type AverageFlux end

# ------------------------------------------------------------
# KennedyGruberGrav  (full refactor, drop-in replacement)
# ------------------------------------------------------------
Base.@kwdef struct KennedyGruberGrav <: AverageFlux end

function (::KennedyGruberGrav)(RhoPos, uPos, vPos, wPos, ThPos, pPos, GPPos)

  @inline function FluxNonLinAver!(flux,
      VLoc, AuxLoc, dXdxILoc,
      K1, Iz1, iD1,   # left state indices  (localidx into @localmem)
      K2, Iz2, iD2,   # right state indices
      ::Val{NV}) where {NV}

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
    m_L1  = dXdxILoc[1, K1, Iz1, iD1]
    m_L2  = dXdxILoc[2, K1, Iz1, iD1]
    m_L3  = dXdxILoc[3, K1, Iz1, iD1]
    m_R1  = dXdxILoc[1, K2, Iz2, iD2]
    m_R2  = dXdxILoc[2, K2, Iz2, iD2]
    m_R3  = dXdxILoc[3, K2, Iz2, iD2]

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


# ============================================================
# Updated kernel — FluxSplitVolumeNonLinV3Kernel!
# Only the inner loop changes; everything else is identical.
# ============================================================

@kernel inbounds = true function FluxSplitVolumeNonLinV3Kernel!(FluxAver!,
    F, @Const(V), @Const(Aux), @Const(dXdxI),
    @Const(DVT), @Const(Glob),
    ::Val{NV}, ::Val{NAUX}) where {NV, NAUX}

  K, Iz, iD     = @index(Local, NTuple)
  _, _, ID, IF  = @index(Global, NTuple)

  TilesDim = @uniform @groupsize()[3]
  M        = @uniform @groupsize()[1]
  Nz       = @uniform @groupsize()[2]
  ND       = @uniform @ndrange()[3]

  VLoc     = @localmem eltype(F)      (M, Nz, TilesDim, NV)
  AuxLoc   = @localmem eltype(F)      (M, Nz, TilesDim, NAUX)
  FLoc     = @localmem eltype(F)      (M, Nz, TilesDim, NV)
  # Store only the 3 metric coefficients needed for the vertical direction.
  # Layout: (3, M, Nz, TilesDim)  — direction index first for coalescing.
  dXdxILoc = @localmem eltype(dXdxI)  (3, M, Nz, TilesDim)

  hTilde = @private eltype(F) (NV,)

  # ---- load phase ----
  if ID <= ND
    ind = Glob[ID, IF]
    @unroll for iaux = 1:NAUX
      AuxLoc[K, Iz, iD, iaux] = Aux[K, Iz, ind, iaux]
    end
    @unroll for iv = 1:NV
      VLoc[K, Iz, iD, iv] = V[K, Iz, ind, iv]
      FLoc[K, Iz, iD, iv] = 0
    end
    @unroll for j = 1:3
      dXdxILoc[j, K, Iz, iD] = dXdxI[3, j, K, ID, Iz, IF]
    end
  end

  @synchronize

  # ---- compute phase — NO @views, NO slices ----
  if ID <= ND
    ind = Glob[ID, IF]
    @unroll for l = 1:M
      # Pass parent arrays + left/right indices directly.
      # FluxAver! reads VLoc[K,Iz,iD,iv] and VLoc[l,Iz,iD,iv] internally.
      FluxAver!(hTilde,
        VLoc, AuxLoc, dXdxILoc,
        K, Iz, iD,    # left  state: vary K
        l, Iz, iD,    # right state: vary l
        Val(NV))
      @unroll for iv = 1:NV
        FLoc[K, Iz, iD, iv] += -DVT[l, K] * hTilde[iv]
      end
    end
    @unroll for iv = 1:NV
      F[K, Iz, ind, iv] += FLoc[K, Iz, iD, iv]
    end
  end
end


# ============================================================
# Quad kernel — same pattern, two directions (I and J)
# ============================================================

@kernel inbounds = true function FluxSplitVolumeNonLinHQuadKernel!(FluxAver!,
    F, @Const(V), @Const(Aux), @Const(dXdxI),
    @Const(DVT), @Const(Glob),
    ::Val{NV}, ::Val{NAUX}) where {NV, NAUX}

  I, J, iz      = @index(Local, NTuple)
  _, _, IZ, IF  = @index(Global, NTuple)

  TilesDim = @uniform @groupsize()[3]
  N        = @uniform @groupsize()[1]
  NZ       = @uniform @ndrange()[3]
  M        = @uniform size(dXdxI, 3)

  VLoc     = @localmem eltype(F)      (N, N, TilesDim, NV)
  AuxLoc   = @localmem eltype(F)      (N, N, TilesDim, NAUX)
  FLoc     = @localmem eltype(F)      (N, N, TilesDim, NV)
  # 2 directions × 3 components
  dXdxILoc = @localmem eltype(dXdxI)  (2, 3, N, N, TilesDim)

  fTilde = @private eltype(F) (NV,)
  gTilde = @private eltype(F) (NV,)

  # ---- load phase ----
  if IZ <= NZ
    K  = mod(IZ - 1, M) + 1
    Iz = div(IZ - 1, M) + 1
    ID = I + (J - 1) * N
    ind = Glob[ID, IF]
    @unroll for iaux = 1:NAUX
      AuxLoc[I, J, iz, iaux] = Aux[K, Iz, ind, iaux]
    end
    @unroll for iv = 1:NV
      VLoc[I, J, iz, iv] = V[K, Iz, ind, iv]
      FLoc[I, J, iz, iv] = 0
    end
    @unroll for j = 1:3
      @unroll for i = 1:2
        dXdxILoc[i, j, I, J, iz] = dXdxI[i, j, K, ID, Iz, IF]
      end
    end
  end

  @synchronize

  # ---- compute phase ----
  if IZ <= NZ
    @unroll for l = 1:N
      # x-direction: left=(I,J,iz), right=(l,J,iz)
      FluxAver!(fTilde,
        VLoc, AuxLoc, dXdxILoc,
        I, J, iz,
        l, J, iz,
        Val(NV))
      # y-direction: left=(I,J,iz), right=(I,l,iz)
      FluxAver!(gTilde,
        VLoc, AuxLoc, dXdxILoc,
        I, J, iz,
        I, l, iz,
        Val(NV))
      @unroll for iv = 1:NV
        FLoc[I, J, iz, iv] += -DVT[l, I] * fTilde[iv] - DVT[l, J] * gTilde[iv]
      end
    end
    K   = mod(IZ - 1, M) + 1
    Iz  = div(IZ - 1, M) + 1
    ID  = I + (J - 1) * N
    ind = Glob[ID, IF]
    @unroll for iv = 1:NV
      F[K, Iz, ind, iv] += FLoc[I, J, iz, iv]
    end
  end
end


# ============================================================
# Migration guide for the other AverageFlux structs
# ============================================================
#
#  Every FluxNonLinAver! closure currently has this signature:
#
#    function FluxNonLinAver!(flux, VL, VR, AuxL, AuxR, m_L, m_R)
#
#  Replace it with:
#
#    function FluxNonLinAver!(flux,
#        VLoc, AuxLoc, dXdxILoc,
#        K1, Iz1, iD1,
#        K2, Iz2, iD2,
#        ::Val{NV}) where {NV}
#
#  Then replace every  VL[xPos]   →  VLoc[K1, Iz1, iD1, xPos]
#                      VR[xPos]   →  VLoc[K2, Iz2, iD2, xPos]
#                      AuxL[pos]  →  AuxLoc[K1, Iz1, iD1, pos]
#                      AuxR[pos]  →  AuxLoc[K2, Iz2, iD2, pos]
#                      m_L[i]     →  dXdxILoc[i, K1, Iz1, iD1]
#                      m_R[i]     →  dXdxILoc[i, K2, Iz2, iD2]
#
#  Structs to migrate (all follow the exact same pattern):
#    KennedyGruber         ✓ (shown above as KennedyGruberGrav)
#    ArtianoGrav
#    ArtianoExGrav
#    KennedyGruberExPGrav
#    ArtianoExner
#    KennedyGruberGravFast
#    KennedyGruberGravLinFast
#    KennedyGruberGravSlow
#    KennedyGruberIEGrav
#
#  NonConservativeFlux structs (BuoyancyFlux etc.) have a different
#  call site (face integrals, not volume kernels) — check before migrating.
# ============================================================
