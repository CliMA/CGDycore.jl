function FluxSplitVolumeNonLinH(FluxAverage,F,U,Aux,DG,dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX,::Grids.Quad)
  backend = get_backend(F)
  DoF = DG.DoF 
  DoFE = DG.DoFE 
  N = DG.OrdPoly + 1
  M = DG.OrdPolyZ + 1
  NzG = min(div(NumberThreadGPU,DoF),M*Nz)
  group = (N,N,NzG,1)
  ndrange = (N,N,M*Nz,NF)
  KFluxSplitVolumeNonLinHKernel! = FluxSplitVolumeNonLinHQuadKernel!(backend,group)
  KFluxSplitVolumeNonLinHKernel!(FluxAverage,F,U,Aux,dXdxI,DG.DVT,DG.Glob,
    Val(NV),Val(NAUX);ndrange=ndrange)
end  

function FluxVolumeNonLinH(Flux,F,U,Aux,DG,dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX,::Grids.Quad)
  backend = get_backend(F)
  DoF = DG.DoF
  DoFE = DG.DoFE
  N = DG.OrdPoly + 1
  M = DG.OrdPolyZ + 1
  NzG = min(div(NumberThreadGPU,DoF),M*Nz)
  group = (N,N,NzG,1)
  ndrange = (N,N,M*Nz,NF)
  KFluxVolumeNonLinHKernel! = FluxVolumeNonLinHQuadKernel!(backend,group)
  KFluxVolumeNonLinHKernel!(Flux,F,U,Aux,dXdxI,DG.DW,DG.Glob,
    Val(NV),Val(NAUX);ndrange=ndrange)
end

function FluxSplitVolumeNonLinH(FluxAverage,F,U,Aux,DG,dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX,::Grids.Tri)
  backend = get_backend(F)
  DoF = DG.DoF
  M = DG.OrdPolyZ + 1
  NzG = min(div(NumberThreadGPU,DoF),M*Nz)
  group = (DoF,NzG,1)
  ndrange = (DoF,M*Nz,NF)
  KFluxSplitVolumeNonLinHKernel! = FluxSplitVolumeNonLinHTriKernel!(backend,group)
  KFluxSplitVolumeNonLinHKernel!(FluxAverage,F,U,Aux,dXdxI,DG.DSx1,DG.DSx2,DG.Glob,
    Val(NV),Val(NAUX);ndrange=ndrange)
end

function FluxSplitVolumeNonLinV(FluxAverage,F,U,Aux,DG,dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX)
  backend = get_backend(F)
  DoF = DG.DoF 
  DoFE = DG.DoFE 
  M = DG.OrdPolyZ + 1
  NDG = min(div(NumberThreadGPU,M*Nz),DoF)
  group = (M,Nz,NDG,1)
  ndrange = (M,Nz,DoF,NF)
  KFluxSplitVolumeNonLinV3Kernel! = FluxSplitVolumeNonLinV3Kernel!(backend,group)
  KFluxSplitVolumeNonLinV3Kernel!(FluxAverage,F,U,Aux,dXdxI,DG.DVZT,DG.Glob,
    Val(NV),Val(NAUX);ndrange=ndrange)
end  

function FluxVolumeNonLinV(Flux,F,U,Aux,DG,dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX)
  backend = get_backend(F)
  DoF = DG.DoF
  DoFE = DG.DoFE
  M = DG.OrdPolyZ + 1
  NDG = min(div(NumberThreadGPU,M*Nz),DoF)
  group = (M,Nz,NDG,1)
  ndrange = (M,Nz,DoF,NF)
  KFluxVolumeNonLinV3Kernel! = FluxVolumeNonLinV3Kernel!(backend,group)
  KFluxVolumeNonLinV3Kernel!(Flux,F,U,Aux,dXdxI,DG.DWZ,DG.Glob,
    Val(NV),Val(NAUX);ndrange=ndrange)
end

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
  dXdxILoc = @localmem eltype(dXdxI)  (1, 3, M, Nz, TilesDim)

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
      dXdxILoc[1, j, K, Iz, iD] = dXdxI[3, j, K, ID, Iz, IF]
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
        Val(1))
      @unroll for iv = 1:NV
        FLoc[K, Iz, iD, iv] += -DVT[l, K] * hTilde[iv]
      end
    end
    @unroll for iv = 1:NV
      F[K, Iz, ind, iv] += FLoc[K, Iz, iD, iv]
    end
  end
end

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

@kernel inbounds = true function FluxSplitVolumeNonLinHTriKernel!(FluxAver!,
    F, @Const(V), @Const(Aux), @Const(dXdxI),
    @Const(Dx1), @Const(Dx2), @Const(Glob),
    ::Val{NV}, ::Val{NAUX}) where {NV, NAUX}

  ID, iz        = @index(Local, NTuple)
  _, IZ, IF     = @index(Global, NTuple)

  TilesDim = @uniform @groupsize()[2]
  N        = @uniform @groupsize()[1]
  NZ       = @uniform @ndrange()[2]
  M        = @uniform size(dXdxI, 3)

  VLoc     = @localmem eltype(F)      (N, TilesDim, NV)
  AuxLoc   = @localmem eltype(F)      (N, TilesDim, NAUX)
  FLoc     = @localmem eltype(F)      (N, TilesDim, NV)
  # 2 horizontal directions × 3 metric components
  dXdxILoc = @localmem eltype(dXdxI)  (2, 3, N, TilesDim)

  fTilde = @private eltype(F) (NV,)
  gTilde = @private eltype(F) (NV,)

  # ---- load phase ----
  if IZ <= NZ
    K   = mod(IZ - 1, M) + 1
    Iz  = div(IZ - 1, M) + 1
    ind = Glob[ID, IF]
    @unroll for iaux = 1:NAUX
      AuxLoc[ID, iz, iaux] = Aux[K, Iz, ind, iaux]
    end
    @unroll for iv = 1:NV
      VLoc[ID, iz, iv] = V[K, Iz, ind, iv]
      FLoc[ID, iz, iv] = 0
    end
    @unroll for j = 1:3
      @unroll for i = 1:2
        dXdxILoc[i, j, ID, iz] = dXdxI[i, j, K, ID, Iz, IF]
      end
    end
  end

  @synchronize

  # ---- compute phase — no @views, no slices, no Val{NV} in FluxAver! ----
  if IZ <= NZ
    @unroll for l = 1:N
      # x1-direction: left=(ID,iz), right=(l,iz)
      FluxAver!(fTilde,
        VLoc, AuxLoc, dXdxILoc,
        ID, iz,
        l,  iz, Val(1))
      # x2-direction: same pair, different metric row
      FluxAver!(gTilde,
        VLoc, AuxLoc, dXdxILoc,
        ID, iz,
        l,  iz, Val(2))
      @unroll for iv = 1:NV
        FLoc[ID, iz, iv] += -Dx1[ID, l] * fTilde[iv] - Dx2[ID, l] * gTilde[iv]
      end
    end
    K   = mod(IZ - 1, M) + 1
    Iz  = div(IZ - 1, M) + 1
    ind = Glob[ID, IF]
    @unroll for iv = 1:NV
      F[K, Iz, ind, iv] += FLoc[ID, iz, iv]
    end
  end
end

@kernel inbounds = true function FluxVolumeNonLinV3Kernel!(Flux,F,@Const(V),@Const(Aux),@Const(dXdxI),
  @Const(DW),@Const(Glob), ::Val{NV}, ::Val{NAUX}) where {NV,NAUX}

  K, Iz, iD,    = @index(Local, NTuple)
  _,_,ID,IF = @index(Global, NTuple)

  TilesDim = @uniform @groupsize()[3]
  M = @uniform @groupsize()[1]
  Nz = @uniform @groupsize()[2]
  ND = @uniform @ndrange()[3]

  Con = @localmem eltype(F) (M,Nz,TilesDim,NV)
  FLoc = @private eltype(F) (3,NV)

  if ID <= ND
    ind = Glob[ID,IF]
    @views Flux(FLoc,V[K,Iz,ind,:],Aux[K,Iz,ind,:])
    @unroll for iv = 1 : NV
      Con[K,Iz,iD,iv] = dXdxI[3,1,K,ID,Iz,IF] * FLoc[1,iv] +
        dXdxI[3,2,K,ID,Iz,IF] * FLoc[2,iv] +
        dXdxI[3,3,K,ID,Iz,IF] * FLoc[3,iv]
    end  
  end
  @synchronize
  if ID <= ND
    ind = Glob[ID,IF]
    @unroll for iv = 1 : NV
      FF = DW[K,1] * Con[1,Iz,iD,iv]
      @unroll for k = 2 : M
        FF += DW[K,k] * Con[k,Iz,iD,iv]
      end
      F[K,Iz,ind,iv] -= FF
    end
  end  
end  


@kernel inbounds = true function FluxVolumeNonLinHQuadKernel!(Flux,F,@Const(V),@Const(Aux),@Const(dXdxI),
  @Const(DW),@Const(Glob), ::Val{NV}, ::Val{NAUX}) where {NV,NAUX}

  I, J, iz,   = @index(Local, NTuple)
  _,_,IZ,IF = @index(Global, NTuple)

  TilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  NZ = @uniform @ndrange()[3]
  M = @uniform size(dXdxI,3)
  nz = @uniform size(dXdxI,5)

  ConX = @localmem eltype(F) (N,N,TilesDim,NV)
  ConY = @localmem eltype(F) (N,N,TilesDim,NV)
  FLoc = @private eltype(F) (3,NV)

  if IZ <= NZ
    K = mod(IZ-1,M)+1
    Iz = div(IZ-1,M) + 1
    ID = I + (J - 1) * N
    ind = Glob[ID,IF]
    @views Flux(FLoc,V[K,Iz,ind,:],Aux[K,Iz,ind,:])
    for iv = 1 : NV
      ConX[I,J,iz,iv] = dXdxI[1,1,K,ID,Iz,IF] * FLoc[1,iv] +
        dXdxI[1,2,K,ID,Iz,IF] * FLoc[2,iv] +
        dXdxI[1,3,K,ID,Iz,IF] * FLoc[3,iv]
      ConY[I,J,iz,iv] = dXdxI[2,1,K,ID,Iz,IF] * FLoc[1,iv] +
        dXdxI[2,2,K,ID,Iz,IF] * FLoc[2,iv] +
        dXdxI[2,3,K,ID,Iz,IF] * FLoc[3,iv]
    end  
  end
  @synchronize

  if IZ <= NZ
    K = mod(IZ-1,M)+1
    Iz = div(IZ-1,M) + 1
    ID = I + (J - 1) * N
    ind = Glob[ID,IF]
    for iv = 1 : NV
      FF = DW[I,1] * ConX[1,J,iz,iv] + DW[J,1] * ConY[I,1,iz,iv]
      for k = 2 : N
        FF = FF + DW[I,k] * ConX[k,J,iz,iv] + DW[J,k] * ConY[I,k,iz,iv]
      end
      F[K,Iz,ind,iv] -= FF
    end
  end
end  


