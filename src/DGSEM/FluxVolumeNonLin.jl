function FluxVolumeNonLinH(Flux,F,U,Aux,DG::FiniteElements.DGQuad,dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX)
  backend = get_backend(F)
  DoF = DG.DoF
  DoFE = DG.DoFE
  N = DG.OrdPoly + 1
  M = DG.V.OrdPolyZ + 1
  NzG = min(div(NumberThreadGPU,DoF),M*Nz)
  group = (N,N,NzG,1)
  ndrange = (N,N,M*Nz,NF)
  KFluxVolumeNonLinHKernel! = FluxVolumeNonLinHQuadKernel!(backend,group)
  KFluxVolumeNonLinHKernel!(Flux,F,U,Aux,dXdxI,DG.DW,DG.Glob,
    Val(NV),Val(NAUX);ndrange=ndrange)
end
function FluxVolumeNonLinV(Flux,F,U,Aux,DG,dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX)
  backend = get_backend(F)
  DoF = DG.DoF
  DoFE = DG.DoFE
  M = DG.V.OrdPolyZ + 1
  # one column has only M*Nz threads (65 for nz=13, order 4), take several columns
  # per workgroup, about 2*NumberThreadGPU threads, so that the warps are filled
  NDG = min(max(div(2*NumberThreadGPU,M*Nz),1),DoF)
  group = (M,Nz,NDG,1)
  ndrange = (M,Nz,DoF,NF)
  KFluxVolumeNonLinV3Kernel! = FluxVolumeNonLinV3Kernel!(backend,group)
  KFluxVolumeNonLinV3Kernel!(Flux,F,U,Aux,dXdxI,DG.V.DWZ,DG.Glob,
    Val(NV),Val(NAUX);ndrange=ndrange)
end

# Load the state of one node from shared memory as SVector; I are the node
# indices, the last index of A runs over the components.
@inline LoadLoc(A, ::Val{N}, I...) where {N} = SVector{N}(ntuple(n -> A[I..., n], Val(N)))
# Metric row J a^dir of one node
@inline LoadMetric(A, ::Val{dir}, I...) where {dir} = SVector{3}(ntuple(j -> A[dir, j, I...], Val(3)))

# Horizontal and vertical volume integral in one kernel per element type
# (hexahedra and prisms), both with half sweep.
function FluxSplitVolumeNonLin(FluxAverageH,FluxAverageV,F,U,Aux,DG::FiniteElements.DGQuad,dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX,Form=ConservativeForm())
  backend = get_backend(F)
  N = DG.OrdPoly + 1
  M = DG.V.OrdPolyZ + 1
  group = (N,N,M,1)
  ndrange = (N,N,M*Nz,NF)
  KFluxSplitVolumeNonLinQuadKernel! = FluxSplitVolumeNonLinQuadKernel!(backend,group)
  KFluxSplitVolumeNonLinQuadKernel!(FluxAverageH,FluxAverageV,F,U,Aux,dXdxI,DG.DVT,DG.V.DVZT,DG.Glob,
    Val(N),Val(M),Val(NV),Val(NAUX),Form;ndrange=ndrange)
end

function FluxSplitVolumeNonLin(FluxAverageH,FluxAverageV,F,U,Aux,DG::FiniteElements.DGTri,dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX,Form=ConservativeForm())
  backend = get_backend(F)
  DoF = DG.DoF
  M = DG.V.OrdPolyZ + 1
  # TZ prisms (DoF x M nodes) on top of each other per workgroup
  TZ = min(max(div(NumberThreadGPU,DoF*M),1),Nz)
  group = (DoF,M,TZ,1)
  ndrange = (DoF,M,Nz,NF)
  KFluxSplitVolumeNonLinTriKernel! = FluxSplitVolumeNonLinTriKernel!(backend,group)
  KFluxSplitVolumeNonLinTriKernel!(FluxAverageH,FluxAverageV,F,U,Aux,dXdxI,DG.DSx1,DG.DSx2,DG.V.DVZT,DG.Glob,
    Val(DoF),Val(M),Val(NV),Val(NAUX),Form;ndrange=ndrange)
end

# One workgroup per hexahedral element (N x N x M nodes), half sweep: the two-point
# flux of a pair of nodes is evaluated once, FluxAver(..,SymAnti()) -> S, A, the
# node gets S + A, its partner S - A. In round k the node I of a line pairs with
# I+k (cyclic, for even N the distance N/2 only from the first half), so all nodes
# evaluate the same number of pairs and the partners of a round are distinct;
# their contributions go through shared memory, one buffer per direction.
# The split operators DVT, DVZT have a zero diagonal (2 D - M^-1 B on Gauss-Lobatto
# nodes), the node with itself is skipped.
# Flux evaluations per node: 2 (N/2) + (M/2) instead of 2 N + M.
@kernel inbounds = true function FluxSplitVolumeNonLinQuadKernel!(FluxAverH, FluxAverV,
    F, @Const(V), @Const(Aux), @Const(dXdxI),
    @Const(DVT), @Const(DVZT), @Const(Glob),
    ::Val{N}, ::Val{M}, ::Val{NV}, ::Val{NAUX}, Form) where {N, M, NV, NAUX}

  I, J, K, _    = @index(Local, NTuple)
  _, _, IZ, IF  = @index(Global, NTuple)

  VLoc     = @localmem eltype(F)      (N, N, M, NV)
  AuxLoc   = @localmem eltype(F)      (N, N, M, NAUX)
  # 3 directions × 3 components
  dXdxILoc = @localmem eltype(dXdxI)  (3, 3, N, N, M)
  PartX    = @localmem eltype(F)      (N, N, M, NV)
  PartY    = @localmem eltype(F)      (N, N, M, NV)
  PartZ    = @localmem eltype(F)      (N, N, M, NV)
  # DVT, DVZT are indexed with the partner index, from the kernel argument (SMatrix)
  # every thread would keep its own copy in local memory
  DVTLoc   = @localmem eltype(F)      (N, N)
  DVZTLoc  = @localmem eltype(F)      (M, M)

  FLoc = @private eltype(F) (NV,)
  # ---- load phase ----
  Iz = div(IZ - 1, M) + 1
  ID = I + (J - 1) * N
  lid = ID + (K - 1) * N * N
  @unroll for n = 1:N*N
    if n == lid
      DVTLoc[n] = DVT[n]
    end
  end
  @unroll for n = 1:M*M
    if n == lid
      DVZTLoc[n] = DVZT[n]
    end
  end
  ind = Glob[ID, IF]
  @unroll for iaux = 1:NAUX
    AuxLoc[I, J, K, iaux] = Aux[K, Iz, ind, iaux]
  end
  VLoc[I, J, K, 1] = V[K, Iz, ind, 1]
  invRho = eltype(F)(1) / VLoc[I, J, K, 1]
  @unroll for iv = 2:NV
    VLoc[I, J, K, iv] = V[K, Iz, ind, iv] * PrimScale(Form, iv, invRho)
  end
  @unroll for j = 1:3
    @unroll for i = 1:3
      dXdxILoc[i, j, I, J, K] = dXdxI[i, j, K, ID, Iz, IF]
    end
  end

  @synchronize

  @unroll for iv = 1:NV
    FLoc[iv] = 0
  end

  # ---- rounds of pairs ----
  for k = 1 : max(N ÷ 2, M ÷ 2)
    VL = LoadLoc(VLoc, Val(NV), I, J, K)
    AuxL = LoadLoc(AuxLoc, Val(NAUX), I, J, K)
    if k <= N ÷ 2
      # x-direction: partner (Ip,J,K)
      Ip = mod(I - 1 + k, N) + 1
      if 2 * k < N || I <= k
        S, A = FluxAverH(VL, LoadLoc(VLoc, Val(NV), Ip, J, K),
          AuxL, LoadLoc(AuxLoc, Val(NAUX), Ip, J, K),
          LoadMetric(dXdxILoc, Val(1), I, J, K), LoadMetric(dXdxILoc, Val(1), Ip, J, K), SymAnti())
        @unroll for iv = 1:NV
          FLoc[iv] -= DVTLoc[Ip, I] * (S[iv] + A[iv])
          PartX[Ip, J, K, iv] = -DVTLoc[I, Ip] * (S[iv] - A[iv])
        end
      else
        @unroll for iv = 1:NV
          PartX[Ip, J, K, iv] = 0
        end
      end
      # y-direction: partner (I,Jp,K)
      Jp = mod(J - 1 + k, N) + 1
      if 2 * k < N || J <= k
        S, A = FluxAverH(VL, LoadLoc(VLoc, Val(NV), I, Jp, K),
          AuxL, LoadLoc(AuxLoc, Val(NAUX), I, Jp, K),
          LoadMetric(dXdxILoc, Val(2), I, J, K), LoadMetric(dXdxILoc, Val(2), I, Jp, K), SymAnti())
        @unroll for iv = 1:NV
          FLoc[iv] -= DVTLoc[Jp, J] * (S[iv] + A[iv])
          PartY[I, Jp, K, iv] = -DVTLoc[J, Jp] * (S[iv] - A[iv])
        end
      else
        @unroll for iv = 1:NV
          PartY[I, Jp, K, iv] = 0
        end
      end
    end
    if k <= M ÷ 2
      # z-direction: partner (I,J,Kp)
      Kp = mod(K - 1 + k, M) + 1
      if 2 * k < M || K <= k
        S, A = FluxAverV(VL, LoadLoc(VLoc, Val(NV), I, J, Kp),
          AuxL, LoadLoc(AuxLoc, Val(NAUX), I, J, Kp),
          LoadMetric(dXdxILoc, Val(3), I, J, K), LoadMetric(dXdxILoc, Val(3), I, J, Kp), SymAnti())
        @unroll for iv = 1:NV
          FLoc[iv] -= DVZTLoc[Kp, K] * (S[iv] + A[iv])
          PartZ[I, J, Kp, iv] = -DVZTLoc[K, Kp] * (S[iv] - A[iv])
        end
      else
        @unroll for iv = 1:NV
          PartZ[I, J, Kp, iv] = 0
        end
      end
    end

    @synchronize

    if k <= N ÷ 2
      @unroll for iv = 1:NV
        FLoc[iv] += PartX[I, J, K, iv] + PartY[I, J, K, iv]
      end
    end
    if k <= M ÷ 2
      @unroll for iv = 1:NV
        FLoc[iv] += PartZ[I, J, K, iv]
      end
    end

    @synchronize
  end

  Iz = div(IZ - 1, M) + 1
  ID = I + (J - 1) * N
  ind = Glob[ID, IF]
  @unroll for iv = 1:NV
    F[K, Iz, ind, iv] += FLoc[iv]
  end
end

# One workgroup per TZ prisms (DoF x M nodes), half sweep as in the hexahedral kernel:
# the horizontal pairs run cyclically over the DoF nodes of the triangle, each pair
# with both metric rows, f = Dx1 f(m1) + Dx2 f(m2); the vertical pairs as there.
# Dx1, Dx2 (DSx1, DSx2) have a zero diagonal, the node with itself is skipped.
# Flux evaluations per node: 2 (DoF/2) + (M/2) instead of 2 DoF + M.
@kernel inbounds = true function FluxSplitVolumeNonLinTriKernel!(FluxAverH, FluxAverV,
    F, @Const(V), @Const(Aux), @Const(dXdxI),
    @Const(Dx1), @Const(Dx2), @Const(DVZT), @Const(Glob),
    ::Val{DoF}, ::Val{M}, ::Val{NV}, ::Val{NAUX}, Form) where {DoF, M, NV, NAUX}

  ID, K, iz, _  = @index(Local, NTuple)
  _, _, Iz, IF  = @index(Global, NTuple)

  TZ = @uniform @groupsize()[3]
  Nz = @uniform @ndrange()[3]

  VLoc     = @localmem eltype(F)      (DoF, M, TZ, NV)
  AuxLoc   = @localmem eltype(F)      (DoF, M, TZ, NAUX)
  # 3 directions × 3 components
  dXdxILoc = @localmem eltype(dXdxI)  (3, 3, DoF, M, TZ)
  PartH    = @localmem eltype(F)      (DoF, M, TZ, NV)
  PartZ    = @localmem eltype(F)      (DoF, M, TZ, NV)
  # DVZT (SMatrix) in shared memory, see the hexahedral kernel
  DVZTLoc  = @localmem eltype(F)      (M, M)

  FLoc = @private eltype(F) (NV,)
  # ---- load phase ----
  lid = ID + (K - 1) * DoF + (iz - 1) * DoF * M
  @unroll for n = 1:M*M
    if n == lid
      DVZTLoc[n] = DVZT[n]
    end
  end
  if Iz <= Nz
    ind = Glob[ID, IF]
    @unroll for iaux = 1:NAUX
      AuxLoc[ID, K, iz, iaux] = Aux[K, Iz, ind, iaux]
    end
    VLoc[ID, K, iz, 1] = V[K, Iz, ind, 1]
    invRho = eltype(F)(1) / VLoc[ID, K, iz, 1]
    @unroll for iv = 2:NV
      VLoc[ID, K, iz, iv] = V[K, Iz, ind, iv] * PrimScale(Form, iv, invRho)
    end
    @unroll for j = 1:3
      @unroll for i = 1:3
        dXdxILoc[i, j, ID, K, iz] = dXdxI[i, j, K, ID, Iz, IF]
      end
    end
  end

  @synchronize

  @unroll for iv = 1:NV
    FLoc[iv] = 0
  end

  # ---- rounds of pairs ----
  for k = 1 : max(DoF ÷ 2, M ÷ 2)
    if Iz <= Nz
      VL = LoadLoc(VLoc, Val(NV), ID, K, iz)
      AuxL = LoadLoc(AuxLoc, Val(NAUX), ID, K, iz)
      if k <= DoF ÷ 2
        # horizontal: partner (Lp,K)
        Lp = mod(ID - 1 + k, DoF) + 1
        if 2 * k < DoF || ID <= k
          VR = LoadLoc(VLoc, Val(NV), Lp, K, iz)
          AuxR = LoadLoc(AuxLoc, Val(NAUX), Lp, K, iz)
          S1, A1 = FluxAverH(VL, VR, AuxL, AuxR,
            LoadMetric(dXdxILoc, Val(1), ID, K, iz), LoadMetric(dXdxILoc, Val(1), Lp, K, iz), SymAnti())
          S2, A2 = FluxAverH(VL, VR, AuxL, AuxR,
            LoadMetric(dXdxILoc, Val(2), ID, K, iz), LoadMetric(dXdxILoc, Val(2), Lp, K, iz), SymAnti())
          @unroll for iv = 1:NV
            FLoc[iv] -= Dx1[ID, Lp] * (S1[iv] + A1[iv]) + Dx2[ID, Lp] * (S2[iv] + A2[iv])
            PartH[Lp, K, iz, iv] = -(Dx1[Lp, ID] * (S1[iv] - A1[iv]) + Dx2[Lp, ID] * (S2[iv] - A2[iv]))
          end
        else
          @unroll for iv = 1:NV
            PartH[Lp, K, iz, iv] = 0
          end
        end
      end
      if k <= M ÷ 2
        # vertical: partner (ID,Kp)
        Kp = mod(K - 1 + k, M) + 1
        if 2 * k < M || K <= k
          S, A = FluxAverV(VL, LoadLoc(VLoc, Val(NV), ID, Kp, iz),
            AuxL, LoadLoc(AuxLoc, Val(NAUX), ID, Kp, iz),
            LoadMetric(dXdxILoc, Val(3), ID, K, iz), LoadMetric(dXdxILoc, Val(3), ID, Kp, iz), SymAnti())
          @unroll for iv = 1:NV
            FLoc[iv] -= DVZTLoc[Kp, K] * (S[iv] + A[iv])
            PartZ[ID, Kp, iz, iv] = -DVZTLoc[K, Kp] * (S[iv] - A[iv])
          end
        else
          @unroll for iv = 1:NV
            PartZ[ID, Kp, iz, iv] = 0
          end
        end
      end
    end

    @synchronize

    if Iz <= Nz
      if k <= DoF ÷ 2
        @unroll for iv = 1:NV
          FLoc[iv] += PartH[ID, K, iz, iv]
        end
      end
      if k <= M ÷ 2
        @unroll for iv = 1:NV
          FLoc[iv] += PartZ[ID, K, iz, iv]
        end
      end
    end

    @synchronize
  end

  if Iz <= Nz
    ind = Glob[ID, IF]
    @unroll for iv = 1:NV
      F[K, Iz, ind, iv] += FLoc[iv]
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
    @unroll for iv = 1 : NV
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
    @unroll for iv = 1 : NV
      FF = DW[I,1] * ConX[1,J,iz,iv] + DW[J,1] * ConY[I,1,iz,iv]
      for k = 2 : N
        FF = FF + DW[I,k] * ConX[k,J,iz,iv] + DW[J,k] * ConY[I,k,iz,iv]
      end
      F[K,Iz,ind,iv] -= FF
    end
  end
end  
