function JacobiCartDG3GPU!(AdaptGrid,X,dXdxI,J,Rotate,FE,F,z,zs,Rad,::Grids.Quad)

  backend = get_backend(X)
  FT = eltype(X)

  NF = size(X,5)
  N = size(FE.xw,1)
  M = size(FE.xwZ,1)
  Nz = size(X,4)

  NzG = min(div(512,N*N*M),Nz)
  group = (N, N, M, NzG, 1)
  ndrange = (N, N, M, Nz, NF)

  KJacobiCartDG3Kernel! = JacobiCartDG3QuadKernel!(backend,group)

  KJacobiCartDG3Kernel!(AdaptGrid,X,dXdxI,J,Rotate,FE.xw,FE.xwZ,FE.DS,FE.DSZ,F,z,zs,Rad,ndrange=ndrange)
end

@kernel inbounds = true function JacobiCartDG3QuadKernel!(AdaptGrid,X,dXdxI,JJ,Rotate,@Const(ksi),@Const(zeta),@Const(D),
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
    hR,_,_ = AdaptGrid(zLoc,zs[ID,IF])
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


function JacobiDG1GPU!(X,dXdxI,J,FE,z)

  backend = get_backend(X)
  FT = eltype(X)

  M = size(FE.xwZ,1)
  Nz = size(X,2)

  NzG = min(div(512,M),Nz)
  group = (M, NzG)
  ndrange = (M, Nz)

  KJacobiDG1Kernel! = JacobiDG1Kernel!(backend,group)

  KJacobiDG1Kernel!(X,dXdxI,J,FE.xwZ,FE.DSZ,z,ndrange=ndrange)
end

@kernel inbounds = true function JacobiDG1Kernel!(X,dXdxI,JJ,@Const(zeta),@Const(DZ),@Const(z))

  gk, gz = @index(Group, NTuple)
  K, iz   = @index(Local, NTuple)
  _,Iz = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  L = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[2]
  H = @uniform z[Nz+1]

  dXdx = @localmem eltype(X) (L,ColumnTilesDim)
  XLoc = @localmem eltype(X) (L,ColumnTilesDim)

  if Iz <= Nz
    z1 = z[Iz]
    z2 = z[Iz+1]
    ksi3 = zeta[K]
    zLoc = eltype(X)(1/2) * ((eltype(X)(1)-ksi3) * z1 + (eltype(X)(1)+ksi3) * z2)
    XLoc[K,iz] = zLoc
    X[K,Iz] = XLoc[K,iz]
  end

  @synchronize
  if Iz <= Nz
    dXdx[K,iz] = DZ[K,1] * XLoc[1,iz]
    for k = 2 : L
      dXdx[K,iz] += DZ[K,k] * XLoc[k,iz]
    end  
    @views JJ[K,Iz] = abs(dXdx[K,iz])
    dXdxI[K,Iz] = eltype(X)(1)
  end
end  


function JacobiSphereDG3GPU!(AdaptGrid,X,dXdxI,J,Rotate,FE,F,z,zs,Rad,::Grids.Quad)

  backend = get_backend(X)
  FT = eltype(X)

  NF = size(X,5)
  N = size(FE.xw,1)
  M = size(FE.xwZ,1)
  Nz = size(X,4)

  NzG = min(div(512,N*N*M),Nz)
  group = (N, N, M, NzG, 1)
  ndrange = (N, N, M, Nz, NF)

  KJacobiSphereDG3Kernel! = JacobiSphereDG3QuadKernel!(backend,group)

  KJacobiSphereDG3Kernel!(AdaptGrid,X,dXdxI,J,Rotate,FE.xw,FE.xwZ,FE.DS,FE.DSZ,F,z,zs,Rad,ndrange=ndrange)
end

@kernel inbounds = true function JacobiSphereDG3QuadKernel!(AdaptGrid,X,dXdxI,JJ,Rotate,@Const(ksi),
  @Const(zeta),@Const(D),@Const(DZ),@Const(F),@Const(z),@Const(zs),Rad)

  gi, gj, gk, gz, gF = @index(Group, NTuple)
  I, J, K, iz   = @index(Local, NTuple)
  _,_,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[4]
  N = @uniform @groupsize()[1]
  L = @uniform @groupsize()[3]
  Nz = @uniform @ndrange()[4]
  NF = @uniform @ndrange()[5]
  H = @uniform z[Nz+1]
  zsLoc = @localmem eltype(zs) (N,N)

  dXdx = @private eltype(X) (3,3)

  if Iz <= Nz && K == 1
    ID = I + (J - 1) * N
    zsLoc[I,J] = zs[ID,IF]
  end

  @synchronize 

  eta = ksi
  if Iz <= Nz
    ID = I + (J - 1) * N
    z1 = z[Iz]
    z2 = z[Iz+1]
    ksi1 = ksi[I]
    ksi2 = ksi[J]
    ksi3 = zeta[K]
    XT1 = eltype(X)(1/4) * (F[1,1,IF] * (eltype(X)(1)-ksi1)*(eltype(X)(1)-ksi2) +
     F[2,1,IF] * (eltype(X)(1)+ksi1)*(eltype(X)(1)-ksi2) +
     F[3,1,IF] * (eltype(X)(1)+ksi1)*(eltype(X)(1)+ksi2) +
     F[4,1,IF] * (eltype(X)(1)-ksi1)*(eltype(X)(1)+ksi2))
    XT2 = eltype(X)(1/4) * (F[1,2,IF] * (eltype(X)(1)-ksi1)*(eltype(X)(1)-ksi2) +
     F[2,2,IF] * (eltype(X)(1)+ksi1)*(eltype(X)(1)-ksi2) +
     F[3,2,IF] * (eltype(X)(1)+ksi1)*(eltype(X)(1)+ksi2) +
     F[4,2,IF] * (eltype(X)(1)-ksi1)*(eltype(X)(1)+ksi2))
    XT3 = eltype(X)(1/4) * (F[1,3,IF] * (eltype(X)(1)-ksi1)*(eltype(X)(1)-ksi2) +
     F[2,3,IF] * (eltype(X)(1)+ksi1)*(eltype(X)(1)-ksi2) +
     F[3,3,IF] * (eltype(X)(1)+ksi1)*(eltype(X)(1)+ksi2) +
     F[4,3,IF] * (eltype(X)(1)-ksi1)*(eltype(X)(1)+ksi2))
    zLoc = eltype(X)(1/2) * ((eltype(X)(1)-ksi3) * z1 + (eltype(X)(1)+ksi3) * z2)
    hR,dhrdz,dhrdzs = AdaptGrid(zLoc,zs[ID,IF])
    r = sqrt(XT1 * XT1 + XT2 * XT2 + XT3 * XT3)
    f = eltype(X)(1) / r^3
    dX1dXT1 = f * (XT2^2 + XT3^2) 
    dX1dXT2 = -f * XT1 * XT2 
    dX1dXT3 = -f * XT1 * XT3
    dX2dXT1 = dX1dXT2
    dX2dXT2 = f * (XT1^2 + XT3^2)
    dX2dXT3 = -f * XT2 * XT3
    dX3dXT1 = dX1dXT3
    dX3dXT2 = dX2dXT3 
    dX3dXT3 = f * (XT1^2 + XT2^2)
    J1 = @SArray([dX1dXT1    dX1dXT2     dX1dXT3
                  dX2dXT1    dX2dXT2     dX2dXT3
                  dX3dXT1    dX3dXT2     dX3dXT3])
    dzsdksi1 = D[I,1] * zsLoc[1,J]
    dzsdksi2 = D[J,1] * zsLoc[I,1]
    for k = 2 : N
      dzsdksi1 += D[I,k] * zsLoc[k,J]
      dzsdksi2 += D[J,k] * zsLoc[I,k]
    end  
    B = @SArray([F[1,1,IF] F[2,1,IF] F[3,1,IF] F[4,1,IF]
                 F[1,2,IF] F[2,2,IF] F[3,2,IF] F[4,2,IF]
                 F[1,3,IF] F[2,3,IF] F[3,3,IF] F[4,3,IF]])
    C = @SArray([eltype(X)(-1)+ksi2  eltype(X)(-1)+ksi1
                 eltype(X)(1)-ksi2  eltype(X)(-1)-ksi1
                 eltype(X)(1)+ksi2   eltype(X)(1)+ksi1
                 eltype(X)(-1)-ksi2   eltype(X)(1)-ksi1])

    dXdx[:,1:2] .= eltype(X)(0.25) * (Rad + hR) * J1 * B * C

    dXdx[1,1] += dzsdksi1 * dhrdzs * XT1 / r 
    dXdx[2,1] += dzsdksi1 * dhrdzs * XT2 / r 
    dXdx[3,1] += dzsdksi1 * dhrdzs * XT3 / r 
    dXdx[1,2] += dzsdksi2 * dhrdzs * XT1 / r 
    dXdx[2,2] += dzsdksi2 * dhrdzs * XT2 / r 
    dXdx[3,2] += dzsdksi2 * dhrdzs * XT3 / r 
    dXdx[1,3] = XT1 * dhrdz * eltype(X)(0.5) / r * (z2 - z1)
    dXdx[2,3] = XT2 * dhrdz * eltype(X)(0.5) / r * (z2 - z1)
    dXdx[3,3] = XT3 * dhrdz * eltype(X)(0.5) / r * (z2 - z1)

    X[ID,K,1,Iz,IF] = XT1 / r * (Rad + hR)
    X[ID,K,2,Iz,IF] = XT2 / r * (Rad + hR)
    X[ID,K,3,Iz,IF] = XT3 / r * (Rad + hR)

    @views JJ[ID,K,Iz,IF] = Det3(dXdx)
    dXdxI[1,1,K,ID,Iz,IF] = dXdx[2,2] * dXdx[3,3] - dXdx[2,3] * dXdx[3,2]
    dXdxI[2,1,K,ID,Iz,IF] = -(dXdx[2,1] * dXdx[3,3] - dXdx[2,3] * dXdx[3,1])
    dXdxI[3,1,K,ID,Iz,IF] = dXdx[2,1] * dXdx[3,2] - dXdx[2,2] * dXdx[3,1]

    dXdxI[1,2,K,ID,Iz,IF] = -(dXdx[1,2] * dXdx[3,3] - dXdx[1,3] * dXdx[3,2])
    dXdxI[2,2,K,ID,Iz,IF] = dXdx[1,1] * dXdx[3,3] - dXdx[1,3] * dXdx[3,1]
    dXdxI[3,2,K,ID,Iz,IF] = -(dXdx[1,1] * dXdx[3,2] - dXdx[1,2] * dXdx[3,1])

    dXdxI[1,3,K,ID,Iz,IF] = dXdx[1,2] * dXdx[2,3] - dXdx[1,3] * dXdx[2,2]
    dXdxI[2,3,K,ID,Iz,IF] = -(dXdx[1,1] * dXdx[2,3] - dXdx[1,3] * dXdx[2,1])
    dXdxI[3,3,K,ID,Iz,IF] = dXdx[1,1] * dXdx[2,2] - dXdx[1,2] * dXdx[2,1]

    lon,lat,_ = cart2sphere(XT1,XT2,XT3)
    Rotate[1,1,K,ID,Iz,IF] = -sin(lon)
    Rotate[2,1,K,ID,Iz,IF] = -sin(lat)*cos(lon)
    Rotate[3,1,K,ID,Iz,IF] =  cos(lat)*cos(lon)

    Rotate[1,2,K,ID,Iz,IF] = cos(lon)
    Rotate[2,2,K,ID,Iz,IF] = -sin(lat)*sin(lon)
    Rotate[3,2,K,ID,Iz,IF] = cos(lat) * sin(lon)

    Rotate[1,3,K,ID,Iz,IF] = 0.0
    Rotate[2,3,K,ID,Iz,IF] = cos(lat)  
    Rotate[3,3,K,ID,Iz,IF] = sin(lat)  

  end
end  


function JacobiSphereDG3GPU!(AdaptGrid,X,dXdxI,J,Rotate,FE,F,z,zs,Rad,::Grids.Tri)

  backend = get_backend(X)
  FT = eltype(X)

  NF = size(X,5)
  N = FE.DoF
  M = size(FE.xwZ,1)
  Nz = size(X,4)

  NzG = min(div(512,N*M),Nz)
  group = (N, M, NzG, 1)
  ndrange = (N, M, Nz, NF)

  KJacobiSphereDGTriKernel! = JacobiSphereDGTriKernel!(backend,group)

  KJacobiSphereDGTriKernel!(AdaptGrid,X,dXdxI,J,Rotate,FE.Glob,FE.ksi,FE.xwZ,FE.Dx1,FE.Dx2,FE.DSZ,F,z,zs,Rad,ndrange=ndrange)
end

@kernel inbounds = true function JacobiSphereDGTriKernel!(AdaptGrid,X,dXdxI,JJ,Rotate,@Const(Glob),@Const(ksi),@Const(zeta),
  @Const(Dx1),@Const(Dx2),@Const(DZ),@Const(F),@Const(z),@Const(zs),Rad)

  ID, K, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[4]
  N = @uniform @groupsize()[1]
  L = @uniform @groupsize()[2]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]
  H = @uniform z[Nz+1]

  dXdx = @private eltype(X) (3,3)


  if Iz <= Nz
    z1 = z[Iz]
    z2 = z[Iz+1]
    ksi1 = ksi[1,ID]
    ksi2 = ksi[2,ID]
    ksi3 = zeta[K]
    XT1 = eltype(X)(1/2) * (F[1,1,IF] * (-ksi1 - ksi2) +
     F[2,1,IF] * (eltype(X)(1)+ksi1) +
     F[3,1,IF] * (eltype(X)(1)+ksi2)) 
    XT2 = eltype(X)(1/2) * (F[1,2,IF] * (-ksi1 - ksi2) +
     F[2,2,IF] * (eltype(X)(1)+ksi1) +
     F[3,2,IF] * (eltype(X)(1)+ksi2)) 
    XT3 = eltype(X)(1/2) * (F[1,3,IF] * (-ksi1 - ksi2) +
     F[2,3,IF] * (eltype(X)(1)+ksi1) +
     F[3,3,IF] * (eltype(X)(1)+ksi2)) 
    zLoc = eltype(X)(1/2) * ((eltype(X)(1)-ksi3) * z1 + (eltype(X)(1)+ksi3) * z2)
    hR,dhrdz,dhrdzs = AdaptGrid(zLoc,zs[ID,IF])
    r = sqrt(XT1 * XT1 + XT2 * XT2 + XT3 * XT3)
    f = eltype(X)(1) / r^3
    dX1dXT1 = f * (XT2^2 + XT3^2) 
    dX1dXT2 = -f * XT1 * XT2 
    dX1dXT3 = -f * XT1 * XT3
    dX2dXT1 = dX1dXT2
    dX2dXT2 = f * (XT1^2 + XT3^2)
    dX2dXT3 = -f * XT2 * XT3
    dX3dXT1 = dX1dXT3
    dX3dXT2 = dX2dXT3 
    dX3dXT3 = f * (XT1^2 + XT2^2)
    J1 = @SArray([dX1dXT1    dX1dXT2     dX1dXT3
                  dX2dXT1    dX2dXT2     dX2dXT3
                  dX3dXT1    dX3dXT2     dX3dXT3])

    dzsdksi1 = Dx1[ID,1] * zs[1,IF]
    dzsdksi2 = Dx2[ID,1] * zs[1,IF]
    for k = 2 : N
      dzsdksi1 += Dx1[ID,k] * zs[k,IF]
      dzsdksi2 += Dx2[ID,k] * zs[k,IF]
    end  

    B = @SArray([F[1,1,IF] F[2,1,IF] F[3,1,IF] 
                 F[1,2,IF] F[2,2,IF] F[3,2,IF]
                 F[1,3,IF] F[2,3,IF] F[3,3,IF]])
    C = @SArray([eltype(X)(-1)  eltype(X)(-1)
                 eltype(X)(1)   eltype(X)(0)
                 eltype(X)(0)   eltype(X)(1)])

    dXdx[:,1:2] .= eltype(X)(0.5) * (Rad + hR) * J1 * B * C

    dXdx[1,1] += dzsdksi1 * dhrdzs * XT1 / r 
    dXdx[2,1] += dzsdksi1 * dhrdzs * XT2 / r 
    dXdx[3,1] += dzsdksi1 * dhrdzs * XT3 / r 
    dXdx[1,2] += dzsdksi2 * dhrdzs * XT1 / r 
    dXdx[2,2] += dzsdksi2 * dhrdzs * XT2 / r 
    dXdx[3,2] += dzsdksi2 * dhrdzs * XT3 / r 
    dXdx[1,3] = XT1 * dhrdz * eltype(X)(0.5) / r * (z2 - z1)
    dXdx[2,3] = XT2 * dhrdz * eltype(X)(0.5) / r * (z2 - z1)
    dXdx[3,3] = XT3 * dhrdz * eltype(X)(0.5) / r * (z2 - z1)

    X[ID,K,1,Iz,IF] = XT1 / r * (Rad + hR)
    X[ID,K,2,Iz,IF] = XT2 / r * (Rad + hR)
    X[ID,K,3,Iz,IF] = XT3 / r * (Rad + hR)

    @views JJ[ID,K,Iz,IF] = Det3(dXdx)
    dXdxI[1,1,K,ID,Iz,IF] = dXdx[2,2] * dXdx[3,3] - dXdx[2,3] * dXdx[3,2]
    dXdxI[2,1,K,ID,Iz,IF] = -(dXdx[2,1] * dXdx[3,3] - dXdx[2,3] * dXdx[3,1])
    dXdxI[3,1,K,ID,Iz,IF] = dXdx[2,1] * dXdx[3,2] - dXdx[2,2] * dXdx[3,1]

    dXdxI[1,2,K,ID,Iz,IF] = -(dXdx[1,2] * dXdx[3,3] - dXdx[1,3] * dXdx[3,2])
    dXdxI[2,2,K,ID,Iz,IF] = dXdx[1,1] * dXdx[3,3] - dXdx[1,3] * dXdx[3,1]
    dXdxI[3,2,K,ID,Iz,IF] = -(dXdx[1,1] * dXdx[3,2] - dXdx[1,2] * dXdx[3,1])

    dXdxI[1,3,K,ID,Iz,IF] = dXdx[1,2] * dXdx[2,3] - dXdx[1,3] * dXdx[2,2]
    dXdxI[2,3,K,ID,Iz,IF] = -(dXdx[1,1] * dXdx[2,3] - dXdx[1,3] * dXdx[2,1])
    dXdxI[3,3,K,ID,Iz,IF] = dXdx[1,1] * dXdx[2,2] - dXdx[1,2] * dXdx[2,1]

    lon,lat,_ = cart2sphere(XT1,XT2,XT3)
    Rotate[1,1,K,ID,Iz,IF] = -sin(lon)
    Rotate[2,1,K,ID,Iz,IF] = -sin(lat)*cos(lon)
    Rotate[3,1,K,ID,Iz,IF] =  cos(lat)*cos(lon)

    Rotate[1,2,K,ID,Iz,IF] = cos(lon)
    Rotate[2,2,K,ID,Iz,IF] = -sin(lat)*sin(lon)
    Rotate[3,2,K,ID,Iz,IF] = cos(lat) * sin(lon)

    Rotate[1,3,K,ID,Iz,IF] = 0.0
    Rotate[2,3,K,ID,Iz,IF] = cos(lat)  
    Rotate[3,3,K,ID,Iz,IF] = sin(lat)  

  end
end  


