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
      AuxLoc[K,iz,iaux] = Aux[K,Iz,iaux]  
    end
    @unroll for iv = 1 : NV
      VLoc[K,iz,iv] = V[K,Iz,iv]  
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
      F[K,Iz,iv] += FLoc[K,iz,iv] 
    end
  end
end  

@kernel inbounds = true function FluxVolumeNonLinVertKernel!(Flux,F,@Const(V),@Const(Aux),
  @Const(DW), ::Val{NV}, ::Val{NAUX}) where {NV,NAUX}

  K, iz    = @index(Local, NTuple)
  _,Iz = @index(Global, NTuple)

  TilesDim = @uniform @groupsize()[2]
  M = @uniform @groupsize()[1]
  NZ = @uniform @ndrange()[2]

  Con = @localmem eltype(F) (M,TilesDim,NV)
  FLoc = @private eltype(F) (NV)

  if Iz <= NZ
    @views Flux(FLoc,V[K,Iz,:],Aux[K,Iz,:])
    @. Con[K,Iz,:] = FLoc
  end
  @synchronize
  if Iz <= NZ
    @unroll for iv = 1 : NV
      FF = DW[K,1] * Con[1,Iz,iv]
      @unroll for k = 2 : M
        FF += DW[K,k] * Con[k,Iz,iv]
      end
      F[K,Iz,iv] = -FF
    end
  end
end  
