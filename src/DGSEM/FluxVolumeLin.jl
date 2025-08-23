@kernel inbounds = true function FluxVolumeLinH3Kernel!(F,@Const(V),@Const(Aux),@Const(dXdxI),
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

  if IZ <= NZ
    K = mod(IZ-1,M)+1
    Iz = div(IZ-1,M) + 1
    ID = I + (J - 1) * N
    ind = Glob[ID,IF]
    p = Aux[K,Iz,ind,1] * V[K,Iz,ind,5]
    Th = Aux[K,Iz,ind,2]
    ConX[I,J,iz,1] = dXdxI[1,1,K,ID,Iz,IF] * V[K,Iz,ind,2] +
      dXdxI[1,2,K,ID,Iz,IF] * V[K,Iz,ind,3] +
      dXdxI[1,3,K,ID,Iz,IF] * V[K,Iz,ind,4]
    ConX[I,J,iz,5] = ConX[I,J,iz,1] * Th  
    ConY[I,J,iz,1] = dXdxI[2,1,K,ID,Iz,IF] * V[K,Iz,ind,2] +
      dXdxI[2,2,K,ID,Iz,IF] * V[K,Iz,ind,3] +
      dXdxI[2,3,K,ID,Iz,IF] * V[K,Iz,ind,4]
    ConY[I,J,iz,5] = ConY[I,J,iz,1] * Th  
    ConX[I,J,iz,2] = dXdxI[1,1,K,ID,Iz,IF] * p
    ConY[I,J,iz,2] = dXdxI[2,1,K,ID,Iz,IF] * p
    ConX[I,J,iz,3] = dXdxI[1,2,K,ID,Iz,IF] * p
    ConY[I,J,iz,3] = dXdxI[2,2,K,ID,Iz,IF] * p
    ConX[I,J,iz,4] = dXdxI[1,3,K,ID,Iz,IF] * p
    ConY[I,J,iz,4] = dXdxI[2,3,K,ID,Iz,IF] * p
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

@kernel inbounds = true function FluxVolumeLinV3Kernel!(F,@Const(V),@Const(Aux),@Const(dXdxI),
  @Const(DW),@Const(Glob), ::Val{NV}, ::Val{NAUX}) where {NV,NAUX}

  K, Iz, iD,    = @index(Local, NTuple)
  _,_,ID,IF = @index(Global, NTuple)

  TilesDim = @uniform @groupsize()[3]
  M = @uniform @groupsize()[1]
  Nz = @uniform @groupsize()[2]
  ND = @uniform @ndrange()[3]

  Con = @localmem eltype(F) (M,Nz,TilesDim,NV)

  if ID <= ND
    ind = Glob[ID,IF]
    p = Aux[K,Iz,ind,1] * V[K,Iz,ind,5]
    Th = Aux[K,Iz,ind,2]
    ind = Glob[ID,IF]
    Con[K,Iz,iD,1] = dXdxI[3,1,K,ID,Iz,IF] * V[K,Iz,ind,2] +
      dXdxI[3,2,K,ID,Iz,IF] * V[K,Iz,ind,3] +
      dXdxI[3,3,K,ID,Iz,IF] * V[K,Iz,ind,4]
    Con[K,Iz,iD,5] = Con[K,Iz,iD,1] * Th
    Con[K,Iz,iD,2] = dXdxI[3,1,K,ID,Iz,IF] * p
    Con[K,Iz,iD,3] = dXdxI[3,2,K,ID,Iz,IF] * p
    Con[K,Iz,iD,4] = dXdxI[3,3,K,ID,Iz,IF] * p
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
