function Vorticity(Vort,U,DG,dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX,::Grids.Quad)
  backend = get_backend(F)
  DoF = DG.DoF
  DoFE = DG.DoFE
  N = DG.OrdPoly + 1 
  M = DG.OrdPolyZ + 1 
  group = (N,N,M,1,1)
  ndrange = (N,N,M,Nz,NF)
  KVorticityQuadKernel! = VorticityQuadKernel!(backend,group)
  KVorticityQuadKernel!(Vort,U,dXdxI,DG.DW,DG.DWZ,DG.Glob,
    Val(N), Val(M), Val(NV);ndrange=ndrange)
end


@kernel inbounds = true function FluxSplitVolumeNonLinHVQuadKernel!(
    F, @Const(V), @Const(dXdxI), @Const(DW),
    @Const(DWZ), @Const(Glob),
    ::Val{N}, ::Val{M}, ::Val{NV}) where {N, M, NV}

  I, J, K      = @index(Local, NTuple)
  _, _, _, IZ, IF  = @index(Global, NTuple)

  NZ       = @uniform @ndrange()[4]

  VLoc     = @localmem eltype(F)      (N, N, M, 3)
  # 2 directions × 3 components
  dXdxILoc = @localmem eltype(dXdxI)  (3, 3, N, N, M)


  # ---- load phase ----
  ID = I + (J - 1) * N
  ind = Glob[ID, IF]
  VLoc[I, J, K, 1] = V[K, IZ, ind, 2]
  VLoc[I, J, K, 2] = V[K, IZ, ind, 3]
  VLoc[I, J, K, 3] = V[K, IZ, ind, 4]
  @unroll for j = 1:3
    @unroll for i = 1:3
      dXdxILoc[i, j, I, J, K] = dXdxI[i, j, K, ID, IZ, IF]
    end
  end

  @synchronize

  # ---- compute phase ----
  @unroll for l = 1:N
    # x-direction: left=(I,J,K), right=(l,J,K)
    FluxAver!(fTilde,
      VLoc, AuxLoc, dXdxILoc,
      I, J, K,
      l, J, K,
      Val(1))
    # y-direction: left=(I,J,K), right=(I,l,K)
    FluxAver!(gTilde,
      VLoc, AuxLoc, dXdxILoc,
      I, J, K,
      I, l, K,
      Val(2))
    @unroll for iv = 1:NV
      FLoc[iv] += -DVT[l, I] * fTilde[iv] - DVT[l, J] * gTilde[iv]
    end
  end  
  @unroll for l = 1:M
    # z-direction: left=(I,J,K), right=(I,J,l)
    FluxAver!(hTilde,
      VLoc, AuxLoc, dXdxILoc,
      I, J, K,
      I, J, l,
      Val(3))
    @unroll for iv = 1:NV
      FLoc[iv] += -DVZT[l, K] * hTilde[iv]
    end
  end  
  ID = I + (J - 1) * N
  ind = Glob[ID, IF]
  @unroll for iv = 1:NV
    F[K, IZ, ind, iv] += FLoc[iv]
  end
end


@inline function Vorticity!(vort,
      VLoc, dXdxILoc,
      K, Iz, iD,   # state indices  (localidx into @localmem)
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
