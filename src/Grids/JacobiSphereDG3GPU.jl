function JacobiCartDG3GPU!(AdaptGrid,X,dXdxI,J,Rotate,FE,F,z,zs,Rad)

  backend = get_backend(X)
  FT = eltype(X)

  NF = size(X,5)
  N = size(FE.xw,1)
  M = size(FE.xwZ,1)
  Nz = size(X,4)

  NzG = min(div(512,N*N*M),Nz)
  group = (N, N, M, NzG, 1)
  ndrange = (N, N, M, Nz, NF)

  KJacobiCartDG3Kernel! = JacobiCartDG3Kernel!(backend,group)

  KJacobiCartDG3Kernel!(AdaptGrid,X,dXdxI,J,Rotate,FE.xw,FE.xwZ,FE.DS,FE.DSZ,F,z,zs,Rad,ndrange=ndrange)
end

@kernel inbounds = true function JacobiCartDG3Kernel!(AdaptGrid,X,dXdxI,JJ,Rotate,@Const(ksi),@Const(zeta),@Const(D),
  @Const(DZ),@Const(F),@Const(z),@Const(zs),Rad)

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
    z1 = z[Iz]
    z2 = z[Iz+1]
    z1 = z[Iz]
    z2 = z[Iz+1]
    ksi1 = ksi[I]
    ksi2 = ksi[J]
    ksi3 = zeta[K]
    X1 = eltype(X)(1/4) * (F[1,1,IF] * (eltype(X)(1)-ksi1)*(eltype(X)(1)-ksi2) +
     F[2,1,IF] * (eltype(X)(1)+ksi1)*(eltype(X)(1)-ksi2) +
     F[3,1,IF] * (eltype(X)(1)+ksi1)*(eltype(X)(1)+ksi2) +
     F[4,1,IF] * (eltype(X)(1)-ksi1)*(eltype(X)(1)+ksi2))
    X2 = eltype(X)(1/4) * (F[1,2,IF] * (eltype(X)(1)-ksi1)*(eltype(X)(1)-ksi2) +
     F[2,2,IF] * (eltype(X)(1)+ksi1)*(eltype(X)(1)-ksi2) +
     F[3,2,IF] * (eltype(X)(1)+ksi1)*(eltype(X)(1)+ksi2) +
     F[4,2,IF] * (eltype(X)(1)-ksi1)*(eltype(X)(1)+ksi2))
    X3 = eltype(X)(1/4) * (F[1,3,IF] * (eltype(X)(1)-ksi1)*(eltype(X)(1)-ksi2) +
     F[2,3,IF] * (eltype(X)(1)+ksi1)*(eltype(X)(1)-ksi2) +
     F[3,3,IF] * (eltype(X)(1)+ksi1)*(eltype(X)(1)+ksi2) +
     F[4,3,IF] * (eltype(X)(1)-ksi1)*(eltype(X)(1)+ksi2))
    zLoc = eltype(X)(1/2) * ((eltype(X)(1)-ksi3) * z1 + (eltype(X)(1)+ksi3) * z2)
    hR = AdaptGrid(zLoc,zs[I,J,IF])
    XLoc[I,J,K,1,iz] = X1
    XLoc[I,J,K,2,iz] = X2
    XLoc[I,J,K,3,iz] = X3 + hR
    X[ID,K,1,Iz,IF] = XLoc[I,J,K,1,iz]
    X[ID,K,2,Iz,IF] = XLoc[I,J,K,2,iz]
    X[ID,K,3,Iz,IF] = XLoc[I,J,K,3,iz]
  end

  @synchronize
  if Iz <= Nz
    ID = I + (J - 1) * N
    dXdx[I,J,K,1,1,iz] = D[I,1] * XLoc[1,J,K,1,iz]
    dXdx[I,J,K,2,1,iz] = D[I,1] * XLoc[1,J,K,2,iz]
    dXdx[I,J,K,3,1,iz] = D[I,1] * XLoc[1,J,K,3,iz]
    dXdx[I,J,K,1,2,iz] = D[J,1] * XLoc[I,1,K,1,iz]
    dXdx[I,J,K,2,2,iz] = D[J,1] * XLoc[I,1,K,2,iz]
    dXdx[I,J,K,3,2,iz] = D[J,1] * XLoc[I,1,K,3,iz]
    dXdx[I,J,K,1,3,iz] = DZ[K,1] * XLoc[I,J,1,1,iz]
    dXdx[I,J,K,2,3,iz] = DZ[K,1] * XLoc[I,J,1,2,iz]
    dXdx[I,J,K,3,3,iz] = DZ[K,1] * XLoc[I,J,1,3,iz]
    for k = 2 : N
      dXdx[I,J,K,1,1,iz] += D[I,k] * XLoc[k,J,K,1,iz]
      dXdx[I,J,K,2,1,iz] += D[I,k] * XLoc[k,J,K,2,iz]
      dXdx[I,J,K,3,1,iz] += D[I,k] * XLoc[k,J,K,3,iz]
      dXdx[I,J,K,1,2,iz] += D[J,k] * XLoc[I,k,K,1,iz]
      dXdx[I,J,K,2,2,iz] += D[J,k] * XLoc[I,k,K,2,iz]
      dXdx[I,J,K,3,2,iz] += D[J,k] * XLoc[I,k,K,3,iz]
    end
    for k = 2 : L
      dXdx[I,J,K,1,3,iz] += DZ[K,k] * XLoc[I,J,k,1,iz]
      dXdx[I,J,K,2,3,iz] += DZ[K,k] * XLoc[I,J,k,2,iz]
      dXdx[I,J,K,3,3,iz] += DZ[K,k] * XLoc[I,J,k,3,iz]
    end  
    @views JJ[ID,K,Iz,IF] = Det3(dXdx[I,J,K,:,:,iz])
    dXdxI[1,1,K,ID,Iz,IF] = dXdx[I,J,K,2,2,iz] * dXdx[I,J,K,3,3,iz] -
      dXdx[I,J,K,2,3,iz] * dXdx[I,J,K,3,2,iz]
    dXdxI[2,1,K,ID,Iz,IF] = -(dXdx[I,J,K,2,1,iz] * dXdx[I,J,K,3,3,iz] -
      dXdx[I,J,K,2,3,iz] * dXdx[I,J,K,3,1,iz])
    dXdxI[3,1,K,ID,Iz,IF] = dXdx[I,J,K,2,1,iz] * dXdx[I,J,K,3,2,iz] -
      dXdx[I,J,K,2,2,iz] * dXdx[I,J,K,3,1,iz]

    dXdxI[1,2,K,ID,Iz,IF] = -(dXdx[I,J,K,1,2,iz] * dXdx[I,J,K,3,3,iz] -
      dXdx[I,J,K,1,3,iz] * dXdx[I,J,K,3,2,iz])
    dXdxI[2,2,K,ID,Iz,IF] = dXdx[I,J,K,1,1,iz] * dXdx[I,J,K,3,3,iz] -
      dXdx[I,J,K,1,3,iz] * dXdx[I,J,K,3,1,iz]
    dXdxI[3,2,K,ID,Iz,IF] = -(dXdx[I,J,K,1,1,iz] * dXdx[I,J,K,3,2,iz] -
      dXdx[I,J,K,1,2,iz] * dXdx[I,J,K,3,1,iz])

    dXdxI[1,3,K,ID,Iz,IF] = dXdx[I,J,K,1,2,iz] * dXdx[I,J,K,2,3,iz] -
      dXdx[I,J,K,1,3,iz] * dXdx[I,J,K,2,2,iz]
    dXdxI[2,3,K,ID,Iz,IF] = -(dXdx[I,J,K,1,1,iz] * dXdx[I,J,K,2,3,iz] -
      dXdx[I,J,K,1,3,iz] * dXdx[I,J,K,2,1,iz])
    dXdxI[3,3,K,ID,Iz,IF] = dXdx[I,J,K,1,1,iz] * dXdx[I,J,K,2,2,iz] -
      dXdx[I,J,K,1,2,iz] * dXdx[I,J,K,2,1,iz]
    lon,lat,_ = cart2sphere(XLoc[I,J,K,1,iz],XLoc[I,J,K,2,iz],XLoc[I,J,K,3,iz])
    Rotate[1,1,K,ID,Iz,IF] = eltype(X)(1.0)
    Rotate[2,1,K,ID,Iz,IF] = eltype(X)(0.0)
    Rotate[3,1,K,ID,Iz,IF] = eltype(X)(0.0)

    Rotate[1,2,K,ID,Iz,IF] = eltype(X)(0.0)
    Rotate[2,2,K,ID,Iz,IF] = eltype(X)(1.0)
    Rotate[3,2,K,ID,Iz,IF] = eltype(X)(0.0)

    Rotate[1,3,K,ID,Iz,IF] = eltype(X)(0.0)
    Rotate[2,3,K,ID,Iz,IF] = eltype(X)(0.0)
    Rotate[3,3,K,ID,Iz,IF] = eltype(X)(1.0)

  end
end  


function JacobiSphereDG3GPU!(AdaptGrid,X,dXdxI,J,Rotate,FE,F,z,zs,Rad)

  backend = get_backend(X)
  FT = eltype(X)

  NF = size(X,5)
  N = size(FE.xw,1)
  M = size(FE.xwZ,1)
  Nz = size(X,4)

  NzG = min(div(512,N*N*2),Nz)
  group = (N, N, M, NzG, 1)
  ndrange = (N, N, M, Nz, NF)

  KJacobiSphereDG3Kernel! = JacobiSphereDG3Kernel!(backend,group)

  KJacobiSphereDG3Kernel!(AdaptGrid,X,dXdxI,J,Rotate,FE.xw,FE.xwZ,FE.DS,FE.DSZ,F,z,zs,Rad,ndrange=ndrange)
end

@kernel inbounds = true function JacobiSphereDG3Kernel!(AdaptGrid,X,dXdxI,JJ,Rotate,@Const(ksi),@Const(zeta),@Const(D),
  @Const(DZ),@Const(F),@Const(z),@Const(zs),Rad)

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
    z1 = z[Iz]
    z2 = z[Iz+1]
    z1 = z[Iz]
    z2 = z[Iz+1]
    ksi1 = ksi[I]
    ksi2 = ksi[J]
    ksi3 = zeta[K]
    X1 = eltype(X)(1/4) * (F[1,1,IF] * (eltype(X)(1)-ksi1)*(eltype(X)(1)-ksi2) +
     F[2,1,IF] * (eltype(X)(1)+ksi1)*(eltype(X)(1)-ksi2) +
     F[3,1,IF] * (eltype(X)(1)+ksi1)*(eltype(X)(1)+ksi2) +
     F[4,1,IF] * (eltype(X)(1)-ksi1)*(eltype(X)(1)+ksi2))
    X2 = eltype(X)(1/4) * (F[1,2,IF] * (eltype(X)(1)-ksi1)*(eltype(X)(1)-ksi2) +
     F[2,2,IF] * (eltype(X)(1)+ksi1)*(eltype(X)(1)-ksi2) +
     F[3,2,IF] * (eltype(X)(1)+ksi1)*(eltype(X)(1)+ksi2) +
     F[4,2,IF] * (eltype(X)(1)-ksi1)*(eltype(X)(1)+ksi2))
    X3 = eltype(X)(1/4) * (F[1,3,IF] * (eltype(X)(1)-ksi1)*(eltype(X)(1)-ksi2) +
     F[2,3,IF] * (eltype(X)(1)+ksi1)*(eltype(X)(1)-ksi2) +
     F[3,3,IF] * (eltype(X)(1)+ksi1)*(eltype(X)(1)+ksi2) +
     F[4,3,IF] * (eltype(X)(1)-ksi1)*(eltype(X)(1)+ksi2))
    zLoc = eltype(X)(1/2) * ((eltype(X)(1)-ksi3) * z1 + (eltype(X)(1)+ksi3) * z2)
    hR = AdaptGrid(zLoc,zs[I,J,IF])
    r = sqrt(X1 * X1 + X2 * X2 + X3 * X3)
    XLoc[I,J,K,1,iz] = X1 / r * (Rad + hR)
    XLoc[I,J,K,2,iz] = X2 / r * (Rad + hR)
    XLoc[I,J,K,3,iz] = X3 / r * (Rad + hR)
    X[ID,K,1,Iz,IF] = XLoc[I,J,K,1,iz]
    X[ID,K,2,Iz,IF] = XLoc[I,J,K,2,iz]
    X[ID,K,3,Iz,IF] = XLoc[I,J,K,3,iz]
  end

  @synchronize
  if Iz <= Nz
    ID = I + (J - 1) * N
    dXdx[I,J,K,1,1,iz] = D[I,1] * XLoc[1,J,K,1,iz]
    dXdx[I,J,K,2,1,iz] = D[I,1] * XLoc[1,J,K,2,iz]
    dXdx[I,J,K,3,1,iz] = D[I,1] * XLoc[1,J,K,3,iz]
    dXdx[I,J,K,1,2,iz] = D[J,1] * XLoc[I,1,K,1,iz]
    dXdx[I,J,K,2,2,iz] = D[J,1] * XLoc[I,1,K,2,iz]
    dXdx[I,J,K,3,2,iz] = D[J,1] * XLoc[I,1,K,3,iz]
    dXdx[I,J,K,1,3,iz] = DZ[K,1] * XLoc[I,J,1,1,iz]
    dXdx[I,J,K,2,3,iz] = DZ[K,1] * XLoc[I,J,1,2,iz]
    dXdx[I,J,K,3,3,iz] = DZ[K,1] * XLoc[I,J,1,3,iz]
    for k = 2 : N
      dXdx[I,J,K,1,1,iz] += D[I,k] * XLoc[k,J,K,1,iz]
      dXdx[I,J,K,2,1,iz] += D[I,k] * XLoc[k,J,K,2,iz]
      dXdx[I,J,K,3,1,iz] += D[I,k] * XLoc[k,J,K,3,iz]
      dXdx[I,J,K,1,2,iz] += D[J,k] * XLoc[I,k,K,1,iz]
      dXdx[I,J,K,2,2,iz] += D[J,k] * XLoc[I,k,K,2,iz]
      dXdx[I,J,K,3,2,iz] += D[J,k] * XLoc[I,k,K,3,iz]
    end
    for k = 2 : L
      dXdx[I,J,K,1,3,iz] += DZ[K,k] * XLoc[I,J,k,1,iz]
      dXdx[I,J,K,2,3,iz] += DZ[K,k] * XLoc[I,J,k,2,iz]
      dXdx[I,J,K,3,3,iz] += DZ[K,k] * XLoc[I,J,k,3,iz]
    end  
    @views JJ[ID,K,Iz,IF] = Det3(dXdx[I,J,K,:,:,iz])
    dXdxI[1,1,K,ID,Iz,IF] = dXdx[I,J,K,2,2,iz] * dXdx[I,J,K,3,3,iz] -
      dXdx[I,J,K,2,3,iz] * dXdx[I,J,K,3,2,iz]
    dXdxI[2,1,K,ID,Iz,IF] = -(dXdx[I,J,K,2,1,iz] * dXdx[I,J,K,3,3,iz] -
      dXdx[I,J,K,2,3,iz] * dXdx[I,J,K,3,1,iz])
    dXdxI[3,1,K,ID,Iz,IF] = dXdx[I,J,K,2,1,iz] * dXdx[I,J,K,3,2,iz] -
      dXdx[I,J,K,2,2,iz] * dXdx[I,J,K,3,1,iz]

    dXdxI[1,2,K,ID,Iz,IF] = -(dXdx[I,J,K,1,2,iz] * dXdx[I,J,K,3,3,iz] -
      dXdx[I,J,K,1,3,iz] * dXdx[I,J,K,3,2,iz])
    dXdxI[2,2,K,ID,Iz,IF] = dXdx[I,J,K,1,1,iz] * dXdx[I,J,K,3,3,iz] -
      dXdx[I,J,K,1,3,iz] * dXdx[I,J,K,3,1,iz]
    dXdxI[3,2,K,ID,Iz,IF] = -(dXdx[I,J,K,1,1,iz] * dXdx[I,J,K,3,2,iz] -
      dXdx[I,J,K,1,2,iz] * dXdx[I,J,K,3,1,iz])

    dXdxI[1,3,K,ID,Iz,IF] = dXdx[I,J,K,1,2,iz] * dXdx[I,J,K,2,3,iz] -
      dXdx[I,J,K,1,3,iz] * dXdx[I,J,K,2,2,iz]
    dXdxI[2,3,K,ID,Iz,IF] = -(dXdx[I,J,K,1,1,iz] * dXdx[I,J,K,2,3,iz] -
      dXdx[I,J,K,1,3,iz] * dXdx[I,J,K,2,1,iz])
    dXdxI[3,3,K,ID,Iz,IF] = dXdx[I,J,K,1,1,iz] * dXdx[I,J,K,2,2,iz] -
      dXdx[I,J,K,1,2,iz] * dXdx[I,J,K,2,1,iz]
    lon,lat,_ = cart2sphere(XLoc[I,J,K,1,iz],XLoc[I,J,K,2,iz],XLoc[I,J,K,3,iz])
    Rotate[1,1,K,ID,Iz,IF] = -sin(lon)
    Rotate[2,1,K,ID,Iz,IF] = -sin(lat)*cos(lon)
    Rotate[3,1,K,ID,Iz,IF] = -sin(lat)*cos(lon)

    Rotate[1,2,K,ID,Iz,IF] = cos(lon)
    Rotate[2,2,K,ID,Iz,IF] = -sin(lat)*sin(lon)
    Rotate[3,2,K,ID,Iz,IF] = cos(lat) * sin(lon)

    Rotate[1,3,K,ID,Iz,IF] = 0.0
    Rotate[2,3,K,ID,Iz,IF] = cos(lat)  
    Rotate[3,3,K,ID,Iz,IF] = sin(lat)  

  end
end  


