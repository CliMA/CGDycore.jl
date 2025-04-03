@kernel inbounds = true function FluxSplitVolumeNonLinV3Kernel!(FluxAver!,F,@Const(V),@Const(Aux),@Const(dXdxI),
  @Const(DVT),@Const(Glob), ::Val{NV}, ::Val{NAUX}) where {NV,NAUX}

  K, Iz, iD,    = @index(Local, NTuple)
  _,_,ID,IF = @index(Global, NTuple)

  TilesDim = @uniform @groupsize()[3]
  M = @uniform @groupsize()[1]
  Nz = @uniform @groupsize()[2]
  ND = @uniform @ndrange()[3]

  VLoc = @localmem eltype(F) (M,Nz,TilesDim,NV)
  AuxLoc = @localmem eltype(F) (M,Nz,TilesDim,NAUX)
  FLoc = @localmem eltype(F) (M,Nz,TilesDim,NV)
  dXdxILoc = @localmem eltype(F) (3,M,Nz,TilesDim)
  hTilde = @private eltype(F) (NV,)

  if ID <= ND
    ind = Glob[ID,IF]  
    @unroll for iaux = 1 : NAUX
      AuxLoc[K,Iz,iD,iaux] = Aux[Iz,K,ind,iaux]  
    end
    @unroll for iv = 1 : NV
      VLoc[K,Iz,iD,iv] = V[Iz,K,ind,iv]  
      FLoc[K,Iz,iD,iv] = 0.0
    end

    @unroll for j = 1 : 3
      dXdxILoc[j,K,Iz,iD] = dXdxI[3,j,K,ID,Iz,IF]
    end  
  end

  @synchronize

  if ID <= ND
    @unroll for l = 1 : M
      @views FluxAver!(hTilde,VLoc[K,Iz,iD,:],VLoc[l,Iz,iD,:],
        AuxLoc[K,Iz,iD,:],AuxLoc[l,Iz,iD,:],
        dXdxILoc[:,K,Iz,iD],dXdxILoc[:,l,Iz,iD])    
      @unroll for iv = 1 : NV
        FLoc[K,Iz,iD,iv] += -DVT[l,K] * hTilde[iv] 
      end  
    end  
    ind = Glob[ID,IF]  
    @unroll for iv = 1 : NV
      F[Iz,K,ind,iv] += FLoc[K,Iz,iD,iv] 
    end
  end
end  

@kernel inbounds = true function FluxSplitVolumeNonLinH3Kernel!(FluxAver!,F,@Const(V),@Const(Aux),@Const(dXdxI),
  @Const(DVT),@Const(Glob), ::Val{NV}, ::Val{NAUX}) where {NV,NAUX}

  I, J, iz,   = @index(Local, NTuple)
  _,_,IZ,IF = @index(Global, NTuple)

  TilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  NZ = @uniform @ndrange()[3]
  M = @uniform size(dXdxI,3)
  nz = @uniform size(dXdxI,5)

  VLoc = @localmem eltype(F) (N,N,TilesDim,NV)
  AuxLoc = @localmem eltype(F) (N,N,TilesDim,NAUX)
  FLoc = @localmem eltype(F) (N,N,TilesDim,NV)
  dXdxILoc = @localmem eltype(F) (2,3,N,N,TilesDim)
  fTilde = @private eltype(F) (NV,)
  gTilde = @private eltype(F) (NV,)

  if IZ <= NZ
    K = mod(IZ-1,M)+1  
    Iz = div(IZ-1,M) + 1
    ID = I + (J - 1) * N  
    ind = Glob[ID,IF]  
    @unroll for iaux = 1 : NAUX
      AuxLoc[I,J,iz,iaux] = Aux[Iz,K,ind,iaux]  
    end
    @unroll for iv = 1 : NV
      VLoc[I,J,iz,iv] = V[Iz,K,ind,iv]  
      FLoc[I,J,iz,iv] = 0.0
    end
    @unroll for j = 1 : 3
      @unroll for i = 1 : 2
        dXdxILoc[i,j,I,J,iz] = dXdxI[i,j,K,ID,Iz,IF]
      end   
    end  
  end

  @synchronize

  if IZ <= NZ
    @unroll for l = 1 : N
      @views FluxAver!(fTilde,VLoc[I,J,iz,:],VLoc[l,J,iz,:],
        AuxLoc[I,J,iz,:],AuxLoc[l,J,iz,:],
        dXdxILoc[1,:,I,J,iz],dXdxILoc[1,:,l,J,iz])    
      @views FluxAver!(gTilde,VLoc[I,J,iz,:],VLoc[I,l,iz,:],
        AuxLoc[I,J,iz,:],AuxLoc[I,l,iz,:],
        dXdxILoc[2,:,I,J,iz],dXdxILoc[2,:,I,l,iz])    
      @unroll for iv = 1 : NV
        FLoc[I,J,iz,iv] += -DVT[l,I] * fTilde[iv] - DVT[l,J] * gTilde[iv]
      end  
    end  
    K = mod(IZ-1,M)+1  
    Iz = div(IZ-1,M) + 1
    ID = I + (J - 1) * N  
    ind = Glob[ID,IF]  
    @unroll for iv = 1 : NV
      F[Iz,K,ind,iv] += FLoc[I,J,iz,iv] 
    end
  end
end  
