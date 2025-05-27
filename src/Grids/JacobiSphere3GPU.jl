#function JacobiSphere3GPU!(X,dXdxI,J,FE::CGTri,F,z,zs,Rad)
#
#  backend = get_backend(X)
#  FT = eltype(X)
#
#end
function JacobiSphere2GPU!(X,dXdxI,J,FE,F,Rad)

  backend = get_backend(X)
  FT = eltype(X)

  NF = size(X,3)
  N = size(FE.xw,1)

  NFG = min(div(256,N*N),NF)
  group = (N, N, NFG)
  ndrange = (N, N, NF)

  KJacobiSphere2Kernel! = JacobiSphere2Kernel!(backend,group)

  KJacobiSphere2Kernel!(X,dXdxI,J,FE.xw,FE.DS,F,Rad,ndrange=ndrange)
end

function JacobiSphere3GPU!(AdaptGrid,X,dXdxI,J,FE,F,z,zs,Rad,Equation::Models.CompressibleShallow)

  backend = get_backend(X)
  FT = eltype(X)

  NF = size(X,5)
  N = size(FE.xw,1)
  Nz = size(X,4)
 
  NzG = min(div(256,N*N*2),Nz)
  group = (N, N, 2, NzG, 1)
  ndrange = (N, N, 2, Nz, NF)

  KJacobiSphere3Kernel! = JacobiSphere3Kernel!(backend,group)

  KJacobiSphere3Kernel!(AdaptGrid,X,dXdxI,J,FE.xw,FE.xwZ,FE.DS,F,z,Rad,zs,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)
end

function JacobiSphere3GPU!(AdaptGrid,X,dXdxI,J,FE,F,z,zs,Rad,Equation::Models.Advection)

  backend = get_backend(X)
  FT = eltype(X)

  NF = size(X,5)
  N = size(FE.xw,1)
  Nz = size(X,4)
 
  NzG = min(div(256,N*N*2),Nz)
  group = (N, N, 2, NzG, 1)
  ndrange = (N, N, 2, Nz, NF)

  KJacobiSphere3Kernel! = JacobiSphere3Kernel!(backend,group)

  KJacobiSphere3Kernel!(AdaptGrid,X,dXdxI,J,FE.xw,FE.xwZ,FE.DS,F,z,Rad,zs,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)
end

function JacobiSphere3GPU!(AdaptGrid,X,dXdxI,J,FE,F,z,zs,Rad,Equation::Models.CompressibleDeep)

  backend = get_backend(X)
  FT = eltype(X)

  NF = size(X,5)
  N = size(FE.xw,1)
  Nz = size(X,4)

  NzG = min(div(256,N*N*2),Nz)
  group = (N, N, 2, NzG, 1)
  ndrange = (N, N, 2, Nz, NF)

  KJacobiSphereDeep3Kernel! = JacobiSphereDeep3Kernel!(backend,group)

  KJacobiSphereDeep3Kernel!(AdaptGrid,X,dXdxI,J,FE.xw,FE.xwZ,FE.DS,F,z,Rad,zs,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)
end

function JacobiSphere3DGGPU!(AdaptGrid,X,dXdxI,J,FE,F,z,zs,Rad)

  backend = get_backend(X)
  FT = eltype(X)

  NF = size(X,5)
  N = size(FE.xw,1)
  M = size(FE.xwZ,1)
  Nz = size(X,4)

  NzG = min(div(256,N*N*M),Nz)
  group = (N, N, M, NzG, 1)
  ndrange = (N, N, M, Nz, NF)

  KJacobiSphereDeep3Kernel! = JacobiSphere3DGKernel!(backend,group)

  KJacobiSphereDeep3Kernel!(AdaptGrid,X,dXdxI,J,FE.xw,FE.xwZ,FE.DS,FE.DSZ,F,z,Rad,zs,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)
end


@kernel inbounds = true function JacobiSphere2Kernel!(X,dXdxI,JJ,@Const(ksi),@Const(D),@Const(F),Rad)

  gi, gj, gF = @index(Group, NTuple)
  I, J,  iF   = @index(Local, NTuple)
  _,_,IF = @index(Global, NTuple)

  FaceTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  NF = @uniform @ndrange()[3]

  dXdxLoc = @localmem eltype(X) (N,N,2,2,FaceTilesDim)

  eta = ksi
  if IF <= NF
    ID = I + (J - 1) * N
    @views JacobiSphere2Loc!(X[ID,:,IF],dXdxLoc[I,J,:,:,iF],ksi[I],eta[J],F[:,:,IF],Rad)
    @views JJ[ID,IF] = Det2(dXdxLoc[I,J,:,:,iF])
    @views Adjunct2!(dXdxI[:,:,ID,IF],dXdxLoc[I,J,:,:,iF])
  end
end

@kernel inbounds = true function JacobiSphere3Kernel!(AdaptGrid,X,dXdxI,JJ,@Const(ksi),@Const(zeta),@Const(D),
  @Const(F),@Const(z),Rad,@Const(zs))

  gi, gj, gk, gz, gF = @index(Group, NTuple)
  I, J, K, iz   = @index(Local, NTuple)
  _,_,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[4]
  N = @uniform @groupsize()[1]
  L = @uniform @groupsize()[3]
  Nz = @uniform @ndrange()[4]
  NF = @uniform @ndrange()[5]
  H = @uniform z[Nz+1]

  hR = @localmem eltype(X) (N,N,L,ColumnTilesDim)
  dXdxLoc = @localmem eltype(X) (N,N,L,3,3,ColumnTilesDim)

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
    hR[I,J,K,iz],_,_ = AdaptGrid(zLoc,zs[ID,IF])
    r = sqrt(X1 * X1 + X2 * X2 + X3 * X3)
    f = Rad / r
    X[ID,K,1,Iz,IF] = X1 / r * (Rad + hR[I,J,K,iz])
    X[ID,K,2,Iz,IF] = X2 / r * (Rad + hR[I,J,K,iz])
    X[ID,K,3,Iz,IF] = X3 / r * (Rad + hR[I,J,K,iz])
    (lam,theta)=cart2sphere(X1,X2,X3)

    DD=@SArray([-sin(lam) cos(lam) eltype(X)(0);
      eltype(X)(0)       eltype(X)(0)     eltype(X)(1)])
    sinlam = sin(lam)
    coslam = cos(lam)
    sinth = sin(theta)
    costh = cos(theta)
    a11 = sinlam * sinlam * costh * costh + sinth * sinth
    a12 = -sinlam * coslam * costh * costh
    a13 = -coslam * sinth * costh
    a21 = a12
    a22 = coslam * coslam * costh * costh + sinth * sinth
    a23 = -sinlam * sinth * costh
    a31 = -coslam * sinth
    a32 = -sinlam * sinth
    a33 = costh
    A = @SArray([a11 a12 a13;
      a21 a22 a23;
      a31 a32 a33])
    B = @SArray([F[1,1,IF] F[2,1,IF] F[3,1,IF] F[4,1,IF];
       F[1,2,IF] F[2,2,IF] F[3,2,IF] F[4,2,IF];
       F[1,3,IF] F[2,3,IF] F[3,3,IF] F[4,3,IF]])

    C = @SArray([-eltype(X)(1)+ksi2  -eltype(X)(1)+ksi1;
              eltype(X)(1)-ksi2  -eltype(X)(1)-ksi1;
              eltype(X)(1)+ksi2   eltype(X)(1)+ksi1;
             -eltype(X)(1)-ksi2   eltype(X)(1)-ksi1])
    @views dXdxLoc[I,J,K,1:2,1:2,iz] .= eltype(X)(1/4) * f * DD * A * B * C
    dXdxLoc[I,J,K,1,3,iz] = eltype(X)(0)
    dXdxLoc[I,J,K,2,3,iz] = eltype(X)(0)
  end  

  @synchronize 

  if Iz <= Nz
    ID = I + (J - 1) * N
    DxhR = D[I,1] * hR[1,J,K,iz]
    DyhR = D[J,1] * hR[I,1,K,iz]
    for k = 2 : N
      DxhR += D[I,k] * hR[k,J,K,iz]
      DyhR += D[J,k] * hR[I,k,K,iz]
    end    
    DzhR = eltype(X)(0.5) * (hR[I,J,2,iz] - hR[I,J,1,iz])

    dXdxLoc[I,J,K,3,1,iz] = DxhR
    dXdxLoc[I,J,K,3,2,iz] = DyhR
    dXdxLoc[I,J,K,3,3,iz] = DzhR
    @views JJ[ID,K,Iz,IF] = Det3(dXdxLoc[I,J,K,:,:,iz])
    dXdxI[1,1,K,ID,Iz,IF] = dXdxLoc[I,J,K,2,2,iz] * dXdxLoc[I,J,K,3,3,iz] - 
      dXdxLoc[I,J,K,2,3,iz] * dXdxLoc[I,J,K,3,2,iz]
    dXdxI[2,1,K,ID,Iz,IF] = -(dXdxLoc[I,J,K,2,1,iz] * dXdxLoc[I,J,K,3,3,iz] - 
      dXdxLoc[I,J,K,2,3,iz] * dXdxLoc[I,J,K,3,1,iz])
    dXdxI[3,1,K,ID,Iz,IF] = dXdxLoc[I,J,K,2,1,iz] * dXdxLoc[I,J,K,3,2,iz] - 
      dXdxLoc[I,J,K,2,2,iz] * dXdxLoc[I,J,K,3,1,iz]

    dXdxI[1,2,K,ID,Iz,IF] = -(dXdxLoc[I,J,K,1,2,iz] * dXdxLoc[I,J,K,3,3,iz] - 
      dXdxLoc[I,J,K,1,3,iz] * dXdxLoc[I,J,K,3,2,iz])
    dXdxI[2,2,K,ID,Iz,IF] = dXdxLoc[I,J,K,1,1,iz] * dXdxLoc[I,J,K,3,3,iz] - 
      dXdxLoc[I,J,K,1,3,iz] * dXdxLoc[I,J,K,3,1,iz]
    dXdxI[3,2,K,ID,Iz,IF] = -(dXdxLoc[I,J,K,1,1,iz] * dXdxLoc[I,J,K,3,2,iz] - 
      dXdxLoc[I,J,K,1,2,iz] * dXdxLoc[I,J,K,3,1,iz])

    dXdxI[1,3,K,ID,Iz,IF] = dXdxLoc[I,J,K,1,2,iz] * dXdxLoc[I,J,K,2,3,iz] - 
      dXdxLoc[I,J,K,1,3,iz] * dXdxLoc[I,J,K,2,2,iz]
    dXdxI[2,3,K,ID,Iz,IF] = -(dXdxLoc[I,J,K,1,1,iz] * dXdxLoc[I,J,K,2,3,iz] - 
      dXdxLoc[I,J,K,1,3,iz] * dXdxLoc[I,J,K,2,1,iz])
    dXdxI[3,3,K,ID,Iz,IF] = dXdxLoc[I,J,K,1,1,iz] * dXdxLoc[I,J,K,2,2,iz] - 
      dXdxLoc[I,J,K,1,2,iz] * dXdxLoc[I,J,K,2,1,iz]
  end
end

@kernel inbounds = true function JacobiSphere3DGKernel!(AdaptGrid,X,dXdxI,JJ,@Const(ksi),@Const(zeta),@Const(D),
  @Const(DZ),@Const(F),@Const(z),Rad,@Const(zs))

  gi, gj, gk, gz, gF = @index(Group, NTuple)
  I, J, K, iz   = @index(Local, NTuple)
  _,_,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[4]
  N = @uniform @groupsize()[1]
  L = @uniform @groupsize()[3]
  Nz = @uniform @ndrange()[4]
  NF = @uniform @ndrange()[5]
  H = @uniform z[Nz+1]

  hR = @localmem eltype(X) (N,N,L,ColumnTilesDim)
  dXdxLoc = @localmem eltype(X) (N,N,L,3,3,ColumnTilesDim)

  if Iz <= Nz
    ID = I + (J - 1) * N
    ksi3 = zeta[K]
    z1 = z[Iz]
    z2 = z[Iz+1]
    zLoc = eltype(X)(1/2) * ((eltype(X)(1)-ksi3) * z1 + (eltype(X)(1)+ksi3) * z2)
    hR[I,J,K,iz] = AdaptGrid(zLoc,zs[I,J,IF])
  end
  @synchronize 
  if Iz <= Nz
    ID = I + (J - 1) * N
    z1 = z[Iz]
    z2 = z[Iz+1]
    ksi1 = ksi[I]
    ksi2 = ksi[J]
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
    r = sqrt(XT1 * XT1 + XT2 * XT2 + XT3 * XT3)
    X[ID,K,1,Iz,IF] = XT1 / r * (Rad + hR[I,J,K,iz])
    X[ID,K,2,Iz,IF] = XT2 / r * (Rad + hR[I,J,K,iz])
    X[ID,K,3,Iz,IF] = XT3 / r * (Rad + hR[I,J,K,iz])
    B = @SArray([F[1,1,IF] F[2,1,IF] F[3,1,IF] F[4,1,IF];
       F[1,2,IF] F[2,2,IF] F[3,2,IF] F[4,2,IF];
       F[1,3,IF] F[2,3,IF] F[3,3,IF] F[4,3,IF]])

    C = @SArray([-eltype(X)(1)+ksi2  -eltype(X)(1)+ksi1;
              eltype(X)(1)-ksi2  -eltype(X)(1)-ksi1;
              eltype(X)(1)+ksi2   eltype(X)(1)+ksi1;
             -eltype(X)(1)-ksi2   eltype(X)(1)-ksi1])
    f = (Rad + hR[I,J,K,iz]) *(XT1^2 + XT2^2 + XT3^2)^(-3/2)
    dX1dXT1 = f * (XT2^2 + XT3^2)
    dX1dXT2= -f * XT1 * XT2
    dX1dXT3= -f * XT1 * XT3
    dX2dXT1 = dX1dXT2
    dX2dXT2 = f * (XT1^2+XT3^2)
    dX2dXT3 = -f * XT2 * XT3
    dX3dXT1 = dX1dXT3
    dX3dXT2 = dX2dXT3
    dX3dXT3 = f*(XT1^2+XT2^2)

    J1  =   @SArray([dX1dXT1    dX1dXT2     dX1dXT3
                    dX2dXT1     dX2dXT2     dX2dXT3
                    dX3dXT1     dX3dXT2     dX3dXT3])
    @views dXdxLoc[I,J,K,1:3,1:2,iz] .= eltype(X)(1/4) * J1 * B * C

    DxhR = D[I,1] * hR[1,J,K,iz]
    DyhR = D[J,1] * hR[I,1,K,iz]
    DzhR = DZ[K,1] * hR[I,J,1,iz]
    for k = 2 : N
      DxhR += D[I,k] * hR[k,J,K,iz]
      DyhR += D[J,k] * hR[I,k,K,iz]
    end    
    for k = 2 : L
      DzhR += DZ[K,k] * hR[I,J,k,iz]
    end    

    dXdxLoc[I,J,K,1,1,iz] += DxhR * XT1 / r
    dXdxLoc[I,J,K,2,1,iz] += DxhR * XT2 / r
    dXdxLoc[I,J,K,3,1,iz] += DxhR * XT3 / r
    dXdxLoc[I,J,K,1,2,iz] += DyhR * XT1 / r
    dXdxLoc[I,J,K,2,2,iz] += DyhR * XT2 / r
    dXdxLoc[I,J,K,3,2,iz] += DyhR * XT3 / r
    dXdxLoc[I,J,K,1,3,iz] = DzhR * XT1 / r
    dXdxLoc[I,J,K,2,3,iz] = DzhR * XT2 / r
    dXdxLoc[I,J,K,3,3,iz] = DzhR * XT3 / r

    @views JJ[ID,K,Iz,IF] = Det3(dXdxLoc[I,J,K,:,:,iz])
    dXdxI[1,1,K,ID,Iz,IF] = dXdxLoc[I,J,K,2,2,iz] * dXdxLoc[I,J,K,3,3,iz] - 
      dXdxLoc[I,J,K,2,3,iz] * dXdxLoc[I,J,K,3,2,iz]
    dXdxI[2,1,K,ID,Iz,IF] = -(dXdxLoc[I,J,K,2,1,iz] * dXdxLoc[I,J,K,3,3,iz] - 
      dXdxLoc[I,J,K,2,3,iz] * dXdxLoc[I,J,K,3,1,iz])
    dXdxI[3,1,K,ID,Iz,IF] = dXdxLoc[I,J,K,2,1,iz] * dXdxLoc[I,J,K,3,2,iz] - 
      dXdxLoc[I,J,K,2,2,iz] * dXdxLoc[I,J,K,3,1,iz]

    dXdxI[1,2,K,ID,Iz,IF] = -(dXdxLoc[I,J,K,1,2,iz] * dXdxLoc[I,J,K,3,3,iz] - 
      dXdxLoc[I,J,K,1,3,iz] * dXdxLoc[I,J,K,3,2,iz])
    dXdxI[2,2,K,ID,Iz,IF] = dXdxLoc[I,J,K,1,1,iz] * dXdxLoc[I,J,K,3,3,iz] - 
      dXdxLoc[I,J,K,1,3,iz] * dXdxLoc[I,J,K,3,1,iz]
    dXdxI[3,2,K,ID,Iz,IF] = -(dXdxLoc[I,J,K,1,1,iz] * dXdxLoc[I,J,K,3,2,iz] - 
      dXdxLoc[I,J,K,1,2,iz] * dXdxLoc[I,J,K,3,1,iz])

    dXdxI[1,3,K,ID,Iz,IF] = dXdxLoc[I,J,K,1,2,iz] * dXdxLoc[I,J,K,2,3,iz] - 
      dXdxLoc[I,J,K,1,3,iz] * dXdxLoc[I,J,K,2,2,iz]
    dXdxI[2,3,K,ID,Iz,IF] = -(dXdxLoc[I,J,K,1,1,iz] * dXdxLoc[I,J,K,2,3,iz] - 
      dXdxLoc[I,J,K,1,3,iz] * dXdxLoc[I,J,K,2,1,iz])
    dXdxI[3,3,K,ID,Iz,IF] = dXdxLoc[I,J,K,1,1,iz] * dXdxLoc[I,J,K,2,2,iz] - 
      dXdxLoc[I,J,K,1,2,iz] * dXdxLoc[I,J,K,2,1,iz]
  end
end

@kernel inbounds = true function JacobiSphereDeep3Kernel!(AdaptGrid,X,dXdxI,JJ,@Const(ksi),@Const(zeta),@Const(D),
  @Const(F),@Const(z),Rad,@Const(zs))

  gi, gj, gk, gz, gF = @index(Group, NTuple)
  I, J, K, iz   = @index(Local, NTuple)
  _,_,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[4]
  N = @uniform @groupsize()[1]
  L = @uniform @groupsize()[3]
  Nz = @uniform @ndrange()[4]
  NF = @uniform @ndrange()[5]
  H = @uniform z[Nz+1]

  hR = @localmem eltype(X) (N,N,L,ColumnTilesDim)
  dXdxLoc = @localmem eltype(X) (N,N,L,3,3,ColumnTilesDim)

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
    hR[I,J,K,iz] = AdaptGrid(zLoc,zs[I,J,IF])
    r = sqrt(X1 * X1 + X2 * X2 + X3 * X3)
    f = Rad / r
    X[ID,K,1,Iz,IF] = X1 / r * (Rad + hR[I,J,K,iz])
    X[ID,K,2,Iz,IF] = X2 / r * (Rad + hR[I,J,K,iz])
    X[ID,K,3,Iz,IF] = X3 / r * (Rad + hR[I,J,K,iz])
    (lam,theta)=cart2sphere(X1,X2,X3)

    DD=@SArray([-sin(lam) * (Rad + hR[I,J,K,iz]) / Rad  cos(lam) * (Rad + hR[I,J,K,iz]) / Rad  eltype(X)(0);
      eltype(X)(0)       eltype(X)(0)     eltype(X)(1) * (Rad + hR[I,J,K,iz]) / Rad ])
    sinlam = sin(lam)
    coslam = cos(lam)
    sinth = sin(theta)
    costh = cos(theta)
    a11 = sinlam * sinlam * costh * costh + sinth * sinth
    a12 = -sinlam * coslam * costh * costh
    a13 = -coslam * sinth * costh
    a21 = a12
    a22 = coslam * coslam * costh * costh + sinth * sinth
    a23 = -sinlam * sinth * costh
    a31 = -coslam * sinth
    a32 = -sinlam * sinth
    a33 = costh
    A = @SArray([a11 a12 a13;
      a21 a22 a23;
      a31 a32 a33])
    B = @SArray([F[1,1,IF] F[2,1,IF] F[3,1,IF] F[4,1,IF];
       F[1,2,IF] F[2,2,IF] F[3,2,IF] F[4,2,IF];
       F[1,3,IF] F[2,3,IF] F[3,3,IF] F[4,3,IF]])

    C = @SArray([-eltype(X)(1)+ksi2  -eltype(X)(1)+ksi1;
              eltype(X)(1)-ksi2  -eltype(X)(1)-ksi1;
              eltype(X)(1)+ksi2   eltype(X)(1)+ksi1;
             -eltype(X)(1)-ksi2   eltype(X)(1)-ksi1])
    @views dXdxLoc[I,J,K,1:2,1:2,iz] .= eltype(X)(1/4) * f * DD * A * B * C
    dXdxLoc[I,J,K,1,3,iz] = eltype(X)(0)
    dXdxLoc[I,J,K,2,3,iz] = eltype(X)(0)
  end  

  @synchronize 

  if Iz <= Nz
    ID = I + (J - 1) * N
    DxhR = D[I,1] * hR[1,J,K,iz]
    DyhR = D[J,1] * hR[I,1,K,iz]
    for k = 2 : N
      DxhR += D[I,k] * hR[k,J,K,iz]
      DyhR += D[J,k] * hR[I,k,K,iz]
    end    
    DzhR = eltype(X)(0.5) * (hR[I,J,2,iz] - hR[I,J,1,iz])

    dXdxLoc[I,J,K,3,1,iz] = DxhR
    dXdxLoc[I,J,K,3,2,iz] = DyhR
    dXdxLoc[I,J,K,3,3,iz] = DzhR
    @views JJ[ID,K,Iz,IF] = Det3(dXdxLoc[I,J,K,:,:,iz])
    dXdxI[1,1,K,ID,Iz,IF] = dXdxLoc[I,J,K,2,2,iz] * dXdxLoc[I,J,K,3,3,iz] - 
      dXdxLoc[I,J,K,2,3,iz] * dXdxLoc[I,J,K,3,2,iz]
    dXdxI[2,1,K,ID,Iz,IF] = -(dXdxLoc[I,J,K,2,1,iz] * dXdxLoc[I,J,K,3,3,iz] - 
      dXdxLoc[I,J,K,2,3,iz] * dXdxLoc[I,J,K,3,1,iz])
    dXdxI[3,1,K,ID,Iz,IF] = dXdxLoc[I,J,K,2,1,iz] * dXdxLoc[I,J,K,3,2,iz] - 
      dXdxLoc[I,J,K,2,2,iz] * dXdxLoc[I,J,K,3,1,iz]

    dXdxI[1,2,K,ID,Iz,IF] = -(dXdxLoc[I,J,K,1,2,iz] * dXdxLoc[I,J,K,3,3,iz] - 
      dXdxLoc[I,J,K,1,3,iz] * dXdxLoc[I,J,K,3,2,iz])
    dXdxI[2,2,K,ID,Iz,IF] = dXdxLoc[I,J,K,1,1,iz] * dXdxLoc[I,J,K,3,3,iz] - 
      dXdxLoc[I,J,K,1,3,iz] * dXdxLoc[I,J,K,3,1,iz]
    dXdxI[3,2,K,ID,Iz,IF] = -(dXdxLoc[I,J,K,1,1,iz] * dXdxLoc[I,J,K,3,2,iz] - 
      dXdxLoc[I,J,K,1,2,iz] * dXdxLoc[I,J,K,3,1,iz])

    dXdxI[1,3,K,ID,Iz,IF] = dXdxLoc[I,J,K,1,2,iz] * dXdxLoc[I,J,K,2,3,iz] - 
      dXdxLoc[I,J,K,1,3,iz] * dXdxLoc[I,J,K,2,2,iz]
    dXdxI[2,3,K,ID,Iz,IF] = -(dXdxLoc[I,J,K,1,1,iz] * dXdxLoc[I,J,K,2,3,iz] - 
      dXdxLoc[I,J,K,1,3,iz] * dXdxLoc[I,J,K,2,1,iz])
    dXdxI[3,3,K,ID,Iz,IF] = dXdxLoc[I,J,K,1,1,iz] * dXdxLoc[I,J,K,2,2,iz] - 
      dXdxLoc[I,J,K,1,2,iz] * dXdxLoc[I,J,K,2,1,iz]
  end
end

@inline function JacobiSphere2Loc!(X,dXdx,ksi1,ksi2,F,Rad)
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
  X3 = quarter * (F[1,3] * (one-ksi1)*(one-ksi2) +
   F[2,3] * (one+ksi1)*(one-ksi2) +
   F[3,3] * (one+ksi1)*(one+ksi2) +
   F[4,3] * (one-ksi1)*(one+ksi2))

  r = sqrt(X1 * X1 + X2 * X2 + X3 * X3)
  f = Rad / r
  X1 = X1 / r
  X2 = X2 / r
  X3 = X3 / r
  (lam,theta)=cart2sphere(X1,X2,X3)

  DD=@SArray([-sin(lam) cos(lam) zero;
      zero       zero     one])

  sinlam = sin(lam)
  coslam = cos(lam)
  sinth = sin(theta)
  costh = cos(theta)
  a11 = sinlam * sinlam * costh * costh + sinth * sinth
  a12 = -sinlam * coslam * costh * costh
  a13 = -coslam * sinth * costh
  a21 = a12
  a22 = coslam * coslam * costh * costh + sinth * sinth
  a23 = -sinlam * sinth * costh
  a31 = -coslam * sinth
  a32 = -sinlam * sinth
  a33 = costh
  A = @SArray([a11 a12 a13;
      a21 a22 a23;
      a31 a32 a33])

  B = @SArray([F[1,1] F[2,1] F[3,1] F[4,1];
       F[1,2] F[2,2] F[3,2] F[4,2];
       F[1,3] F[2,3] F[3,3] F[4,3]])

  C = @SArray([-one+ksi2  -one+ksi1;
              one-ksi2  -one-ksi1;
              one+ksi2   one+ksi1;
             -one-ksi2   one-ksi1])
  D = quarter * f * DD * A * B * C
  dXdx[1,1] = D[1,1]	    
  dXdx[1,2] = D[1,2]	    
  dXdx[2,1] = D[2,1]	    
  dXdx[2,2] = D[2,2]	    
  X[1] = X1 * Rad 
  X[2] = X2 * Rad 
  X[3] = X3 * Rad

end

@inline function Det2(A)
  A[1,1] * A[2,2]  - A[1,2] * A[2,1] 
end  

@inline function Adjunct2!(Ad,A)
#   A[1,1] A[1,2] A[1,3]
#   A[2,1] A[2,2] A[2,3]
#   A[3,1] A[3,2] A[3,3]

  Ad[1,1] = A[2,2] 
  Ad[2,1] = -A[2,1]
  
  Ad[1,2] = -A[1,2] 
  Ad[2,2] = A[1,1] 

end  

@inline function Det3(A)
  A[1,1] * (A[2,2] * A[3,3] - A[2,3] * A[3,2]) -
  A[1,2] * (A[2,1] * A[3,3] - A[2,3] * A[3,1]) +
  A[1,3] * (A[2,1] * A[3,2] - A[2,2] * A[3,1]) 
end  
