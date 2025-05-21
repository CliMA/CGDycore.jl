@kernel inbounds = true function FluxSplitVolumeNonLinVertKernel!(FluxAver!,F,@Const(V),@Const(Aux),@Const(dXdxI),
  @Const(DVT),@Const(Glob), ::Val{NV}, ::Val{NAUX}) where {NV,NAUX}

  K, Iz    = @index(Local, NTuple)

  TilesDim = @uniform @groupsize()[3]
  M = @uniform @groupsize()[1]
  Nz = @uniform @groupsize()[2]

  VLoc = @localmem eltype(F) (M,Nz,TilesDim,NV)
  AuxLoc = @localmem eltype(F) (M,Nz,TilesDim,NAUX)
  FLoc = @localmem eltype(F) (M,Nz,TilesDim,NV)
  hTilde = @private eltype(F) (NV,)

  if Iz <= Nz
    @unroll for iaux = 1 : NAUX
      AuxLoc[K,Iz,iaux] = Aux[Iz,K,iaux]  
    end
    @unroll for iv = 1 : NV
      VLoc[K,Iz,iv] = V[Iz,K,iv]  
      FLoc[K,Iz,iv] = 0.0
    end

    @unroll for j = 1 : 3
      dXdxILoc[j,K,Iz] = dXdxI[j,K,Iz]
    end  
  end

  @synchronize

  if Iz <= Nz
    @unroll for l = 1 : M
      @views FluxAver!(hTilde,VLoc[K,Iz,:],VLoc[l,Iz,:],
        AuxLoc[K,Iz,:],AuxLoc[l,Iz,:],
        dXdxILoc[K,Iz],dXdxILoc[l,Iz])    
      @unroll for iv = 1 : NV
        FLoc[K,Iz,iv] += -DVT[l,K] * hTilde[iv] 
      end  
    end  
    @unroll for iv = 1 : NV
      F[Iz,K,iv] += FLoc[K,Iz,iv] 
    end
  end
end  

