function JacobiDG3GPU!(X,dXdxI,J,FE,F,z,zs)

  backend = get_backend(X)
  FT = eltype(X)

  NF = size(X,5)
  N = size(FE.xw,1)
  Nz = size(X,4)

  NzG = min(div(512,N*N*2),Nz)
  group = (N, N, 2, NzG, 1)
  ndrange = (N, N, 2, Nz, NF)

  KJacobiDG3Kernel! = JacobiDG3Kernel!(backend,group)

  KJacobiDG3Kernel!(X,dXdxI,J,FE.xw,FE.xwZ,FE.DS,FE.DSZ,F,z,zs,ndrange=ndrange)
end

@kernel function JacobiDG3Kernel!(X,dXdxI,JJ,@Const(ksi),@Const(zeta),@Const(D),
  @Const(DZ),@Const(F),@Const(z),@Const(zs))

  gi, gj, gk, gz, gF = @index(Group, NTuple)
  I, J, K, iz   = @index(Local, NTuple)
  _,_,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[4]
  N = @uniform @groupsize()[1]
  L = @uniform @groupsize()[3]
  Nz = @uniform @ndrange()[4]
  NF = @uniform @ndrange()[5]
  H = @uniform z[Nz+1]

  dXdx = @localmem eltype(X) (N,N,L,3,3,ColumnTilesDim)
  XLoc = @localmem eltype(X) (N,N,L,3,ColumnTilesDim)

  eta = ksi
  if Iz <= Nz
    ID = I + (J - 1) * N
    @inbounds z1 = z[Iz]
    @inbounds z2 = z[Iz+1]
    @views @inbounds JacobiDG3Loc!(XLoc[I,J,K,:,iz],
      ksi[I],eta[J],zeta[K],F[:,:,IF],z1,z2,H,zs[I,J,IF])
    X[ID,K,1,Iz,IF] = XLoc[I,J,K,1,iz]
    X[ID,K,2,Iz,IF] = XLoc[I,J,K,2,iz]
    X[ID,K,3,Iz,IF] = XLoc[I,J,K,3,iz]
  end

  @synchronize
  if Iz <= Nz
    ID = I + (J - 1) * N
    @inbounds dXdx[I,J,K,1,1,iz] = D[I,1] * XLoc[1,J,K,1,iz]
    @inbounds dXdx[I,J,K,2,1,iz] = D[I,1] * XLoc[1,J,K,2,iz]
    @inbounds dXdx[I,J,K,3,1,iz] = D[I,1] * XLoc[1,J,K,3,iz]
    @inbounds dXdx[I,J,K,1,2,iz] = D[J,1] * XLoc[I,1,K,1,iz]
    @inbounds dXdx[I,J,K,2,2,iz] = D[J,1] * XLoc[I,1,K,2,iz]
    @inbounds dXdx[I,J,K,3,2,iz] = D[J,1] * XLoc[I,1,K,3,iz]
    @inbounds dXdx[I,J,K,1,3,iz] = DZ[K,1] * XLoc[I,J,1,1,iz]
    @inbounds dXdx[I,J,K,2,3,iz] = DZ[K,1] * XLoc[I,J,1,2,iz]
    @inbounds dXdx[I,J,K,3,3,iz] = DZ[K,1] * XLoc[I,J,1,3,iz]
    for k = 2 : N
      @inbounds dXdx[I,J,K,1,1,iz] += D[I,k] * XLoc[k,J,K,1,iz]
      @inbounds dXdx[I,J,K,2,1,iz] += D[I,k] * XLoc[k,J,K,2,iz]
      @inbounds dXdx[I,J,K,3,1,iz] += D[I,k] * XLoc[k,J,K,3,iz]
      @inbounds dXdx[I,J,K,1,2,iz] += D[J,k] * XLoc[I,k,K,1,iz]
      @inbounds dXdx[I,J,K,2,2,iz] += D[J,k] * XLoc[I,k,K,2,iz]
      @inbounds dXdx[I,J,K,3,2,iz] += D[J,k] * XLoc[I,k,K,3,iz]
    end
    for k = 2 : 2
      @inbounds dXdx[I,J,K,1,3,iz] += DZ[K,k] * XLoc[I,J,k,1,iz]
      @inbounds dXdx[I,J,K,2,3,iz] += DZ[K,k] * XLoc[I,J,k,2,iz]
      @inbounds dXdx[I,J,K,3,3,iz] += DZ[K,k] * XLoc[I,J,k,3,iz]
    end  
    @views @inbounds JJ[ID,K,Iz,IF] = Det3(dXdx[I,J,K,:,:,iz])
    @views @inbounds Adjunct3!(dXdxI[:,:,K,ID,Iz,IF],dXdx[I,J,K,:,:,iz])
  end
end  

@inline function JacobiDG3Loc!(X,ksi1,ksi2,ksi3,F,z1,z2,H,zs)
  zero = eltype(X)(0)
  one = eltype(X)(1)
  half = eltype(X)(1/2)
  quarter = eltype(X)(1/4)
  X1 = quarter * (F[1,1] * (one-ksi1)*(one-ksi2) +
   F[2,1] * (one+ksi1)*(one-ksi2) +
   F[3,1] * (one+ksi1)*(one+ksi2) +
   F[4,1] * (one-ksi1)*(one+ksi2))
  X2 = quarter * (F[1,2] * (one-ksi1)*(one-ksi2) +
   F[2,2] * (one+ksi1)*(one-ksi2) +
   F[3,2] * (one+ksi1)*(one+ksi2) +
   F[4,2] * (one-ksi1)*(one+ksi2))
  zLoc = half * ((one-ksi3) * z1 + (one+ksi3) * z2)
  hR = zLoc + (H - zLoc) * zs / H
  X[1] = X1 
  X[2] = X2 
  X[3] = hR 
end  

