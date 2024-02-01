@kernel function MomentumVectorInvariantCoriolisKernel!(F,@Const(U),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(X),@Const(MRho),@Const(M),@Const(Glob),CoriolisFun)

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N  
  @inbounds ind = Glob[ID,IF]

  uCol = @localmem eltype(F) (N,N,ColumnTilesDim)
  vCol = @localmem eltype(F) (N,N,ColumnTilesDim)
  wCol = @localmem eltype(F) (N,N,ColumnTilesDim+1)
  tempuZ = @localmem eltype(F) (N,N,2,ColumnTilesDim+2)
  tempvZ = @localmem eltype(F) (N,N,2,ColumnTilesDim+2)
  tempwZ = @localmem eltype(F) (N,N,2,ColumnTilesDim+2)
  tempuZ1 = @localmem eltype(F) (N,N,2,ColumnTilesDim)
  tempuZ2 = @localmem eltype(F) (N,N,2,ColumnTilesDim)
  tempvZ1 = @localmem eltype(F) (N,N,2,ColumnTilesDim)
  tempvZ2 = @localmem eltype(F) (N,N,2,ColumnTilesDim)
  tempwZ1 = @localmem eltype(F) (N,N,2,ColumnTilesDim)
  tempwZ2 = @localmem eltype(F) (N,N,2,ColumnTilesDim)

  if Iz <= Nz
    @inbounds uCol[I,J,iz] = U[Iz,ind,2]
    @inbounds vCol[I,J,iz] = U[Iz,ind,3]
    @inbounds wCol[I,J,iz+1] = U[Iz,ind,4]
    if iz == 1 && Iz == 1
      @inbounds wCol[I,J,1] = -(dXdxI[3,1,1,ID,1,IF] * U[Iz,ind,2] +
        dXdxI[3,2,1,ID,1,IF] * U[Iz,ind,3]) / dXdxI[3,3,1,ID,1,IF]
    elseif iz == 1
      @inbounds wCol[I,J,1] = U[Iz-1,ind,4]
    end   
  end 
  @synchronize 

  ID = I + (J - 1) * N  
  @inbounds ind = Glob[ID,IF]
  if Iz <= Nz
    # Dz*(dx33*v - dx32*w)
    @inbounds tempuZ[I,J,1,iz+1] = dXdxI[3,3,1,ID,Iz,IF] * vCol[I,J,iz] - dXdxI[3,2,1,ID,Iz,IF] * wCol[I,J,iz]
    @inbounds tempuZ[I,J,2,iz+1] = dXdxI[3,3,2,ID,Iz,IF] * vCol[I,J,iz] - dXdxI[3,2,2,ID,Iz,IF] * wCol[I,J,iz+1]
    # Dz*(dx33*u - dx31*w)
    @inbounds tempvZ[I,J,1,iz+1] = dXdxI[3,3,1,ID,Iz,IF] * uCol[I,J,iz] - dXdxI[3,1,1,ID,Iz,IF] * wCol[I,J,iz]
    @inbounds tempvZ[I,J,2,iz+1] = dXdxI[3,3,2,ID,Iz,IF] * uCol[I,J,iz] - dXdxI[3,1,2,ID,Iz,IF] * wCol[I,J,iz+1]
    # Dz*(dx32*u - dx31*v)
    @inbounds tempwZ[I,J,1,iz+1] = dXdxI[3,2,1,ID,Iz,IF] * uCol[I,J,iz] - dXdxI[3,1,1,ID,Iz,IF] * vCol[I,J,iz]
    @inbounds tempwZ[I,J,2,iz+1] = dXdxI[3,2,2,ID,Iz,IF] * uCol[I,J,iz] - dXdxI[3,1,2,ID,Iz,IF] * vCol[I,J,iz] 
    if iz == 1 && Iz > 1
      um1 = U[Iz-1,ind,2]  
      vm1 = U[Iz-1,ind,3]  
      @inbounds tempuZ[I,J,2,1] = dXdxI[3,3,2,ID,Iz-1,IF] * vm1 - dXdxI[3,2,2,ID,Iz-1,IF] * wCol[I,J,1]
      @inbounds tempvZ[I,J,2,1] = dXdxI[3,3,2,ID,Iz-1,IF] * um1 - dXdxI[3,1,2,ID,Iz-1,IF] * wCol[I,J,1]
      @inbounds tempwZ[I,J,2,1] = dXdxI[3,2,2,ID,Iz-1,IF] * um1 - dXdxI[3,1,2,ID,Iz-1,IF] * vm1 
    end  
    if iz == ColumnTilesDim && Iz < Nz
      up1 = U[Iz+1,ind,2]  
      vp1 = U[Iz+1,ind,3]  
      @inbounds tempuZ[I,J,1,iz+2] = dXdxI[3,3,1,ID,Iz+1,IF] * vp1 - dXdxI[3,2,1,ID,Iz+1,IF] * wCol[I,J,iz+1]
      @inbounds tempvZ[I,J,1,iz+2] = dXdxI[3,3,1,ID,Iz+1,IF] * up1 - dXdxI[3,1,1,ID,Iz+1,IF] * wCol[I,J,iz+1]
      @inbounds tempwZ[I,J,1,iz+2] = dXdxI[3,2,1,ID,Iz+1,IF] * up1 - dXdxI[3,1,1,ID,Iz+1,IF] * vp1 
    end  
  end  
  @synchronize 

  ID = I + (J - 1) * N  
  @inbounds ind = Glob[ID,IF]
  if Iz <= Nz

# uDot = - v*(Dx*(dx12*u - dx11*v) + Dy*(dx22*u - dx21*v) + Dz*(dx32*u - dx31*v))  
#        - w*(Dx*(dx13*u - dx11*w) + Dy*(dx23*u - dx21*w) + Dz*(dx33*u - dx31*w))
# vDot =   u*(Dx*(dx12*u - dx11*v) + Dy*(dx22*u - dx21*v) + Dz*(dx32*u - dx31*v))
#        - w*(Dx*(dx13*v - dx12*w) + Dy*(dx23*v - dx22*w) + Dz*(dx33*v - dx32*w))
# wDot =   u*(Dx*(dx13*u - dx11*w) + Dy*(dx23*u - dx21*w) + Dz*(dx33*u - dx31*w)) 
#          v*(Dx*(dx13*v - dx12*w) + Dy*(dx23*v - dx22*w) + Dz*(dx33*v - dx32*w))
#   U = Dx*(dx13*v - dx12*w) + Dy*(dx23*v - dx22*w) + Dz*(dx33*v - dx32*w)    
#   V = Dx*(dx13*u - dx11*w) + Dy*(dx23*u - dx21*w) + Dz*(dx33*u - dx31*w)
#   W = Dx*(dx12*u - dx11*v) + Dy*(dx22*u - dx21*v) + Dz*(dx32*u - dx31*v)
    @inbounds tempuZ1[I,J,1,iz] = - dXdxI[1,2,1,ID,Iz,IF] * wCol[I,J,iz]
    @inbounds tempuZ1[I,J,2,iz] = - dXdxI[1,2,2,ID,Iz,IF] * wCol[I,J,iz+1]
#   DerivativeX!(U,temp,D)
    @inbounds tempuZ2[I,J,1,iz] = - dXdxI[2,2,1,ID,Iz,IF] * wCol[I,J,iz]
    @inbounds tempuZ2[I,J,2,iz] = - dXdxI[2,2,2,ID,Iz,IF] * wCol[I,J,iz+1]
#   DerivativeY!(U,temp,D)
    @inbounds tempvZ1[I,J,1,iz] = - dXdxI[1,1,1,ID,Iz,IF] * wCol[I,J,iz]
    @inbounds tempvZ1[I,J,2,iz] = - dXdxI[1,1,2,ID,Iz,IF] * wCol[I,J,iz+1]
#   DerivativeX!(V,temp,D)
    @inbounds tempvZ2[I,J,1,iz] = - dXdxI[2,1,1,ID,Iz,IF] * wCol[I,J,iz]
    @inbounds tempvZ2[I,J,2,iz] = - dXdxI[2,1,2,ID,Iz,IF] * wCol[I,J,iz+1]
#   DerivativeY!(V,temp,D)
    @inbounds tempwZ1[I,J,1,iz] = dXdxI[1,2,1,ID,Iz,IF] * uCol[I,J,iz] - dXdxI[1,1,1,ID,Iz,IF] * vCol[I,J,iz]
    @inbounds tempwZ1[I,J,2,iz] = dXdxI[1,2,2,ID,Iz,IF] * uCol[I,J,iz] - dXdxI[1,1,2,ID,Iz,IF] * vCol[I,J,iz]
#   DerivativeX!(W,temp,D)
    @inbounds tempwZ2[I,J,1,iz] = dXdxI[2,2,1,ID,Iz,IF] * uCol[I,J,iz] - dXdxI[2,1,1,ID,Iz,IF] * vCol[I,J,iz]
    @inbounds tempwZ2[I,J,2,iz] = dXdxI[2,2,2,ID,Iz,IF] * uCol[I,J,iz] - dXdxI[2,1,2,ID,Iz,IF] * vCol[I,J,iz]
#   DerivativeY!(W,temp,D)
    
  end
  @synchronize 

  ID = I + (J - 1) * N  
  @inbounds ind = Glob[ID,IF]
  if Iz <= Nz
    FU = eltype(F)(0)  
    FV = eltype(F)(0)  
    FW1 = eltype(F)(0)  
    FW2 = eltype(F)(0)  
    if Iz > 1
      @inbounds FluxUZ = eltype(F)(0.5) * (tempuZ[I,J,1,iz+1] - tempuZ[I,J,2,iz])
      @inbounds FluxVZ = eltype(F)(0.5) * (tempvZ[I,J,1,iz+1] - tempvZ[I,J,2,iz])
      @inbounds FluxWZ = eltype(F)(0.5) * (tempwZ[I,J,1,iz+1] - tempwZ[I,J,2,iz])
      @inbounds FU += (-vCol[I,J,iz] * FluxWZ - wCol[I,J,iz] * FluxVZ)
      @inbounds FV += (uCol[I,J,iz] * FluxWZ - wCol[I,J,iz] * FluxUZ)
      @inbounds FW1 += ( uCol[I,J,iz] * FluxVZ + vCol[I,J,iz] * FluxUZ)
    end
    if Iz < Nz
      @inbounds FluxUZ = eltype(F)(0.5) * (tempuZ[I,J,1,iz+2] - tempuZ[I,J,2,iz+1])
      @inbounds FluxVZ = eltype(F)(0.5) * (tempvZ[I,J,1,iz+2] - tempvZ[I,J,2,iz+1])
      @inbounds FluxWZ = eltype(F)(0.5) * (tempwZ[I,J,1,iz+2] - tempwZ[I,J,2,iz+1])
      @inbounds FU += (-vCol[I,J,iz] * FluxWZ - wCol[I,J,iz+1] * FluxVZ)
      @inbounds FV += ( uCol[I,J,iz] * FluxWZ - wCol[I,J,iz+1] * FluxUZ)
      @inbounds FW2 += ( uCol[I,J,iz] * FluxVZ + vCol[I,J,iz] * FluxUZ)
    end
    @inbounds U1 = eltype(F)(0.5) * (tempuZ[I,J,2,iz+1] - tempuZ[I,J,1,iz+1])
    @inbounds U2 = eltype(F)(0.5) * (tempuZ[I,J,2,iz+1] - tempuZ[I,J,1,iz+1])
    @inbounds V1 = eltype(F)(0.5) * (tempvZ[I,J,2,iz+1] - tempvZ[I,J,1,iz+1])
    @inbounds V2 = eltype(F)(0.5) * (tempvZ[I,J,2,iz+1] - tempvZ[I,J,1,iz+1])
    @inbounds W1 = eltype(F)(0.5) * (tempwZ[I,J,2,iz+1] - tempwZ[I,J,1,iz+1])
    @inbounds W2 = eltype(F)(0.5) * (tempwZ[I,J,2,iz+1] - tempwZ[I,J,1,iz+1])
    @inbounds U1 += D[I,1] * tempuZ1[1,J,1,iz] + D[J,1] * tempuZ2[I,1,1,iz]
    @inbounds U2 += D[I,1] * tempuZ1[1,J,2,iz] + D[J,1] * tempuZ2[I,1,2,iz]
    @inbounds V1 += D[I,1] * tempvZ1[1,J,1,iz] + D[J,1] * tempvZ2[I,1,1,iz]
    @inbounds V2 += D[I,1] * tempvZ1[1,J,2,iz] + D[J,1] * tempvZ2[I,1,2,iz]
    @inbounds W1 += D[I,1] * tempwZ1[1,J,1,iz] + D[J,1] * tempwZ2[I,1,1,iz]
    @inbounds W2 += D[I,1] * tempwZ1[1,J,2,iz] + D[J,1] * tempwZ2[I,1,2,iz]
    for k = 2 : N
      @inbounds U1 += D[I,k] * tempuZ1[k,J,1,iz] + D[J,k] * tempuZ2[I,k,1,iz]
      @inbounds U2 += D[I,k] * tempuZ1[k,J,2,iz] + D[J,k] * tempuZ2[I,k,2,iz]
      @inbounds V1 += D[I,k] * tempvZ1[k,J,1,iz] + D[J,k] * tempvZ2[I,k,1,iz]
      @inbounds V2 += D[I,k] * tempvZ1[k,J,2,iz] + D[J,k] * tempvZ2[I,k,2,iz]
      @inbounds W1 += D[I,k] * tempwZ1[k,J,1,iz] + D[J,k] * tempwZ2[I,k,1,iz]
      @inbounds W2 += D[I,k] * tempwZ1[k,J,2,iz] + D[J,k] * tempwZ2[I,k,2,iz]
    end  

#   Coriolis
    x = eltype(F)(0.5) * (X[ID,1,1,Iz,IF] + X[ID,2,1,Iz,IF])
    y = eltype(F)(0.5) * (X[ID,1,2,Iz,IF] + X[ID,2,2,Iz,IF])
    z = eltype(F)(0.5) * (X[ID,1,3,Iz,IF] + X[ID,2,3,Iz,IF])
    FuCor, FvCor, FwCor = CoriolisFun(x,y,z,uCol[I,J,iz],vCol[I,J,iz],wCol[I,J,iz],wCol[I,J,iz+1])
    FU += FuCor * (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF])
    FV += FvCor * (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF])
    FW1 += FwCor * JJ[ID,1,Iz,IF]
    FW2 += FwCor * JJ[ID,2,Iz,IF]

    @inbounds @atomic F[Iz,ind,2] += (-vCol[I,J,iz] * W1 - wCol[I,J,iz] * V1 -
      vCol[I,J,iz] * W2 - wCol[I,J,iz+1] * V2 + FU) / M[Iz,ind]
    @inbounds @atomic F[Iz,ind,3] += (uCol[I,J,iz] * W1 - wCol[I,J,iz] * U1 +
      uCol[I,J,iz] * W2 - wCol[I,J,iz+1] * U2 + FV) / M[Iz,ind]
    RhoCol = U[Iz,ind,1]  
    if Iz > 1  
      @inbounds @atomic F[Iz-1,ind,4] += RhoCol * (uCol[I,J,iz] * V1 + vCol[I,J,iz] * U1 + FW1) / MRho[Iz-1,ind]
    end
    if Iz < Nz
      @inbounds @atomic F[Iz,ind,4] += RhoCol * (uCol[I,J,iz] * V2 + vCol[I,J,iz] * U2 + FW2) / MRho[Iz,ind]
    end  
  end
end
