@kernel inbounds = true function FluxSplitVolumeNonLinVertKernel!(FluxAver!,F,@Const(V),@Const(Aux),@Const(dXdxI),
  @Const(DVT), ::Val{NV}, ::Val{NAUX}) where {NV,NAUX}

  K, iz    = @index(Local, NTuple)
  _, Iz    = @index(Global, NTuple)

  TilesDim = @uniform @groupsize()[2]
  M = @uniform @groupsize()[1]
  Nz = @uniform @groupsize()[2]

  VLoc = @localmem eltype(F) (M,TilesDim,NV)
  AuxLoc = @localmem eltype(F) (M,TilesDim,NAUX)
  dXdxILoc = @localmem eltype(F) (M,TilesDim)
  FLoc = @localmem eltype(F) (M,TilesDim,NV)
  hTilde = @private eltype(F) (NV,)

  if Iz <= Nz
    @unroll for iaux = 1 : NAUX
      AuxLoc[K,iz,iaux] = Aux[Iz,K,iaux]  
    end
    @unroll for iv = 1 : NV
      VLoc[K,iz,iv] = V[Iz,K,iv]  
      FLoc[K,iz,iv] = 0.0
    end

    dXdxILoc[K,iz] = dXdxI[K,Iz]
  end

  @synchronize

  if Iz <= Nz
    @unroll for l = 1 : M
      @views FluxAver!(hTilde,VLoc[K,iz,:],VLoc[l,iz,:],
        AuxLoc[K,iz,:],AuxLoc[l,iz,:],
        dXdxILoc[K,iz],dXdxILoc[l,iz])    
      @unroll for iv = 1 : NV
        FLoc[K,iz,iv] += -DVT[l,K] * hTilde[iv] 
      end  
    end  
    @unroll for iv = 1 : NV
      F[Iz,K,iv] += FLoc[K,iz,iv] 
    end
  end
end  

