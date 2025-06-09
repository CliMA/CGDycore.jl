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
      AuxLoc[K,Iz,iD,iaux] = Aux[K,Iz,ind,iaux]  
    end
    @unroll for iv = 1 : NV
      VLoc[K,Iz,iD,iv] = V[K,Iz,ind,iv]  
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
      F[K,Iz,ind,iv] += FLoc[K,Iz,iD,iv] 
    end
  end
end  

@kernel inbounds = true function FluxSplitVolumeNonLinHQuadKernel!(FluxAver!,F,@Const(V),@Const(Aux),@Const(dXdxI),
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
      AuxLoc[I,J,iz,iaux] = Aux[K,Iz,ind,iaux]  
    end
    @unroll for iv = 1 : NV
      VLoc[I,J,iz,iv] = V[K,Iz,ind,iv]  
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
      F[K,Iz,ind,iv] += FLoc[I,J,iz,iv] 
    end
  end
end  

@kernel inbounds = true function FluxSplitVolumeNonLinHTriKernel!(FluxAver!,F,@Const(V),@Const(Aux),@Const(dXdxI),
  @Const(Dx1),@Const(Dx2),@Const(Glob), ::Val{NV}, ::Val{NAUX}) where {NV,NAUX}

  ID, iz,   = @index(Local, NTuple)
  _,IZ,IF = @index(Global, NTuple)

  TilesDim = @uniform @groupsize()[2]
  N = @uniform @groupsize()[1]
  NZ = @uniform @ndrange()[2]
  M = @uniform size(dXdxI,3)
  nz = @uniform size(dXdxI,5)

  VLoc = @localmem eltype(F) (N,TilesDim,NV)
  AuxLoc = @localmem eltype(F) (N,TilesDim,NAUX)
  FLoc = @localmem eltype(F) (N,TilesDim,NV)
  dXdxILoc = @localmem eltype(F) (2,3,N,TilesDim)
  fTilde = @private eltype(F) (NV,)
  gTilde = @private eltype(F) (NV,)

  if IZ <= NZ
    K = mod(IZ-1,M)+1  
    Iz = div(IZ-1,M) + 1
    ind = Glob[ID,IF]  
    @unroll for iaux = 1 : NAUX
      AuxLoc[ID,iz,iaux] = Aux[K,Iz,ind,iaux]  
    end
    @unroll for iv = 1 : NV
      VLoc[ID,iz,iv] = V[K,Iz,ind,iv]  
      FLoc[ID,iz,iv] = 0.0
    end
    @unroll for j = 1 : 3
      @unroll for i = 1 : 2
        dXdxILoc[i,j,ID,iz] = dXdxI[i,j,K,ID,Iz,IF]
      end   
    end  
  end

  @synchronize

  if IZ <= NZ
    @unroll for l = 1 : N
      @views FluxAver!(fTilde,VLoc[ID,iz,:],VLoc[l,iz,:],
        AuxLoc[ID,iz,:],AuxLoc[l,iz,:],
        dXdxILoc[1,:,ID,iz],dXdxILoc[1,:,l,iz])    
      @views FluxAver!(gTilde,VLoc[ID,iz,:],VLoc[l,iz,:],
        AuxLoc[ID,iz,:],AuxLoc[l,iz,:],
        dXdxILoc[2,:,ID,iz],dXdxILoc[2,:,l,iz])    
      @unroll for iv = 1 : NV
        FLoc[ID,iz,iv] += -Dx1[ID,l] * fTilde[iv] - Dx2[ID,l] * gTilde[iv]
      end  
    end  
    K = mod(IZ-1,M)+1  
    Iz = div(IZ-1,M) + 1
    ind = Glob[ID,IF]  
    @unroll for iv = 1 : NV
      F[K,Iz,ind,iv] += FLoc[ID,iz,iv] 
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
    @. Con[K,Iz,iD,:] = dXdxI[3,1,K,ID,Iz,IF] * FLoc[1,:] +
      dXdxI[3,2,K,ID,Iz,IF] * FLoc[2,:] +
      dXdxI[3,3,K,ID,Iz,IF] * FLoc[3,:]
  end
  @synchronize
  if ID <= ND
    ind = Glob[ID,IF]
    @unroll for iv = 1 : NV
      FF = DW[K,1] * Con[1,Iz,iD,iv]
      @unroll for k = 2 : M
        FF += DW[K,k] * Con[k,Iz,iD,iv]
      end
      F[K,Iz,ind,iv] += FF
    end
  end  
end  


@kernel inbounds = true function FluxVolumeNonLinH3Kernel!(Flux,F,@Const(V),@Const(Aux),@Const(dXdxI),
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
    @. ConX[I,J,iz,:] = dXdxI[1,1,K,ID,Iz,IF] * FLoc[1,:] +
      dXdxI[1,2,K,ID,Iz,IF] * FLoc[2,:] +
      dXdxI[1,3,K,ID,Iz,IF] * FLoc[3,:]
    @. ConY[I,J,iz,:] = dXdxI[2,1,K,ID,Iz,IF] * FLoc[1,:] +
      dXdxI[2,2,K,ID,Iz,IF] * FLoc[2,:] +
      dXdxI[2,3,K,ID,Iz,IF] * FLoc[3,:]
  end
  @synchronize

  if IZ <= NZ
    K = mod(IZ-1,M)+1
    Iz = div(IZ-1,M) + 1
    ID = I + (J - 1) * N
    ind = Glob[ID,IF]
    @unroll for iv = 1 : NV
      FF = DW[I,1] * ConX[1,J,iz,iv] + DW[1,J] * ConY[I,1,iz,iv]
      @unroll for k = 2 : M
        FF = DW[I,k] * ConX[k,J,iz,iv] + DW[k,J] * ConY[I,k,iz,iv]
      end
      F[K,Iz,ind,iv] += FF
    end
  end
end  


