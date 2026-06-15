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
        Val(1))
      # y-direction: left=(I,J,iz), right=(I,l,iz)
      FluxAver!(gTilde,
        VLoc, AuxLoc, dXdxILoc,
        I, J, iz,
        I, l, iz,
        Val(2))
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
