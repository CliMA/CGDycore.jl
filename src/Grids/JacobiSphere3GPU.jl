#function JacobiSphere3GPU!(X,dXdxI,J,FE::CGTri,F,z,zs,Rad)
#
#  backend = get_backend(X)
#  FT = eltype(X)
#
#end
function JacobiSphere2GPU!(X,dXdx,dXdxI,J,FE,F,Rad)

  backend = get_backend(X)
  FT = eltype(X)

  NF = size(X,3)
  N = size(FE.xw,1)

  NFG = min(div(256,N*N),NF)
  group = (N, N, NFG)
  ndrange = (N, N, NF)

  KJacobiSphere2Kernel! = JacobiSphere2Kernel!(backend,group)

  KJacobiSphere2Kernel!(X,dXdx,dXdxI,J,FE.xw,FE.DS,F,Rad,ndrange=ndrange)
end

function JacobiSphere3GPU!(AdaptGrid,X,dXdx,dXdxI,J,FE,F,z,zs,Rad,Equation::Models.CompressibleShallow)

  backend = get_backend(X)
  FT = eltype(X)

  NF = size(X,5)
  N = size(FE.xw,1)
  Nz = size(X,4)
 
  NzG = min(div(256,N*N*2),Nz)
  group = (N, N, 2, NzG, 1)
  ndrange = (N, N, 2, Nz, NF)

  KJacobiSphere3Kernel! = JacobiSphere3Kernel!(backend,group)

  KJacobiSphere3Kernel!(AdaptGrid,X,dXdx,dXdxI,J,FE.xw,FE.xwZ,FE.DS,F,z,Rad,zs,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)
end

function JacobiSphere3GPU!(AdaptGrid,X,dXdx,dXdxI,J,FE,F,z,zs,Rad,Equation::Models.CompressibleDeep)

  backend = get_backend(X)
  FT = eltype(X)

  NF = size(X,5)
  N = size(FE.xw,1)
  Nz = size(X,4)

  NzG = min(div(256,N*N*2),Nz)
  group = (N, N, 2, NzG, 1)
  ndrange = (N, N, 2, Nz, NF)

  KJacobiSphereDeep3Kernel! = JacobiSphereDeep3Kernel!(backend,group)

  KJacobiSphereDeep3Kernel!(AdaptGrid,X,dXdx,dXdxI,J,FE.xw,FE.xwZ,FE.DS,F,z,Rad,zs,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)
end

@kernel function JacobiSphere2Kernel!(X,dXdx,dXdxI,JJ,@Const(ksi),@Const(D),@Const(F),Rad)

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
    @views @inbounds JacobiSphere2Loc!(X[ID,:,IF],dXdxLoc,[I,J,:,:,iF],ksi[I],eta[J],F[:,:,IF],Rad)
    @views @inbounds JJ[ID,IF] = Det2(dXdxLoc[I,J,:,:,iF])
    @views @inbounds Adjunct2!(dXdxI[:,:,ID,IF],dXdxLoc[I,J,:,:,iF])
    dXdx[1,1,K,ID,Iz,IF] = dXdxLoc[I,J,K,1,1,iz] 
    dXdx[1,2,K,ID,Iz,IF] = dXdxLoc[I,J,K,1,2,iz] 
    dXdx[2,1,K,ID,Iz,IF] = dXdxLoc[I,J,K,2,1,iz] 
    dXdx[2,2,K,ID,Iz,IF] = dXdxLoc[I,J,K,2,2,iz] 
  end
end

@kernel function JacobiSphere3Kernel!(AdaptGrid,X,dXdx,dXdxI,JJ,@Const(ksi),@Const(zeta),@Const(D),
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
    @inbounds z1 = z[Iz]
    @inbounds z2 = z[Iz+1]
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
#   hR[I,J,K,iz], D33 = GalChen(zLoc,H,zs[I,J,IF])
    hR[I,J,K,iz], D33 = AdaptGrid(zLoc,zs[I,J,IF])
    D33 = eltype(X)(1/2) * D33*(z2-z1)
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
    dXdxLoc[I,J,K,3,3,iz] = D33
  end  

  @synchronize 

  if Iz <= Nz
    ID = I + (J - 1) * N
    @inbounds DxhR = D[I,1] * hR[1,J,K,iz]
    @inbounds DyhR = D[J,1] * hR[I,1,K,iz]
    for k = 2 : N
      @inbounds DxhR += D[I,k] * hR[k,J,K,iz]
      @inbounds DyhR += D[J,k] * hR[I,k,K,iz]
    end    

    @inbounds dXdxLoc[I,J,K,3,1,iz] = DxhR
    @inbounds dXdxLoc[I,J,K,3,2,iz] = DyhR
    @views @inbounds JJ[ID,K,Iz,IF] = Det3(dXdxLoc[I,J,K,:,:,iz])
    @inbounds dXdxI[1,1,K,ID,Iz,IF] = dXdxLoc[I,J,K,2,2,iz] * dXdxLoc[I,J,K,3,3,iz] - 
      dXdxLoc[I,J,K,2,3,iz] * dXdxLoc[I,J,K,3,2,iz]
    @inbounds dXdxI[2,1,K,ID,Iz,IF] = -(dXdxLoc[I,J,K,2,1,iz] * dXdxLoc[I,J,K,3,3,iz] - 
      dXdxLoc[I,J,K,2,3,iz] * dXdxLoc[I,J,K,3,1,iz])
    @inbounds dXdxI[3,1,K,ID,Iz,IF] = dXdxLoc[I,J,K,2,1,iz] * dXdxLoc[I,J,K,3,2,iz] - 
      dXdxLoc[I,J,K,2,2,iz] * dXdxLoc[I,J,K,3,1,iz]

    @inbounds dXdxI[1,2,K,ID,Iz,IF] = -(dXdxLoc[I,J,K,1,2,iz] * dXdxLoc[I,J,K,3,3,iz] - 
      dXdxLoc[I,J,K,1,3,iz] * dXdxLoc[I,J,K,3,2,iz])
    @inbounds dXdxI[2,2,K,ID,Iz,IF] = dXdxLoc[I,J,K,1,1,iz] * dXdxLoc[I,J,K,3,3,iz] - 
      dXdxLoc[I,J,K,1,3,iz] * dXdxLoc[I,J,K,3,1,iz]
    @inbounds dXdxI[3,2,K,ID,Iz,IF] = -(dXdxLoc[I,J,K,1,1,iz] * dXdxLoc[I,J,K,3,2,iz] - 
      dXdxLoc[I,J,K,1,2,iz] * dXdxLoc[I,J,K,3,1,iz])

    @inbounds dXdxI[1,3,K,ID,Iz,IF] = dXdxLoc[I,J,K,1,2,iz] * dXdxLoc[I,J,K,2,3,iz] - 
      dXdxLoc[I,J,K,1,3,iz] * dXdxLoc[I,J,K,2,2,iz]
    @inbounds dXdxI[2,3,K,ID,Iz,IF] = -(dXdxLoc[I,J,K,1,1,iz] * dXdxLoc[I,J,K,2,3,iz] - 
      dXdxLoc[I,J,K,1,3,iz] * dXdxLoc[I,J,K,2,1,iz])
    @inbounds dXdxI[3,3,K,ID,Iz,IF] = dXdxLoc[I,J,K,1,1,iz] * dXdxLoc[I,J,K,2,2,iz] - 
      dXdxLoc[I,J,K,1,2,iz] * dXdxLoc[I,J,K,2,1,iz]
    dXdx[1,1,K,ID,Iz,IF] = dXdxLoc[I,J,K,1,1,iz] 
    dXdx[1,2,K,ID,Iz,IF] = dXdxLoc[I,J,K,1,2,iz] 
    dXdx[1,3,K,ID,Iz,IF] = dXdxLoc[I,J,K,1,3,iz] 
    dXdx[2,1,K,ID,Iz,IF] = dXdxLoc[I,J,K,2,1,iz] 
    dXdx[2,2,K,ID,Iz,IF] = dXdxLoc[I,J,K,2,2,iz] 
    dXdx[2,3,K,ID,Iz,IF] = dXdxLoc[I,J,K,2,3,iz] 
    dXdx[3,1,K,ID,Iz,IF] = dXdxLoc[I,J,K,3,1,iz] 
    dXdx[3,2,K,ID,Iz,IF] = dXdxLoc[I,J,K,3,2,iz] 
    dXdx[3,3,K,ID,Iz,IF] = dXdxLoc[I,J,K,3,3,iz] 
  end
end

@kernel function JacobiSphereDeep3Kernel!(AdaptGrid,X,dXdx,dXdxI,JJ,@Const(ksi),@Const(zeta),@Const(D),
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
    @inbounds z1 = z[Iz]
    @inbounds z2 = z[Iz+1]
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
#   hR[I,J,K,iz], D33 = GalChen(zLoc,H,zs[I,J,IF])
    hR[I,J,K,iz], D33 = AdaptGrid(zLoc,zs[I,J,IF])
    D33 = eltype(X)(1/2) * D33*(z2-z1)
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
    dXdxLoc[I,J,K,3,3,iz] = D33
  end  

  @synchronize 

  if Iz <= Nz
    ID = I + (J - 1) * N
    @inbounds DxhR = D[I,1] * hR[1,J,K,iz]
    @inbounds DyhR = D[J,1] * hR[I,1,K,iz]
    for k = 2 : N
      @inbounds DxhR += D[I,k] * hR[k,J,K,iz]
      @inbounds DyhR += D[J,k] * hR[I,k,K,iz]
    end    

    @inbounds dXdxLoc[I,J,K,3,1,iz] = DxhR
    @inbounds dXdxLoc[I,J,K,3,2,iz] = DyhR
    @views @inbounds JJ[ID,K,Iz,IF] = Det3(dXdxLoc[I,J,K,:,:,iz])
    @inbounds dXdxI[1,1,K,ID,Iz,IF] = dXdxLoc[I,J,K,2,2,iz] * dXdxLoc[I,J,K,3,3,iz] - 
      dXdxLoc[I,J,K,2,3,iz] * dXdxLoc[I,J,K,3,2,iz]
    @inbounds dXdxI[2,1,K,ID,Iz,IF] = -(dXdxLoc[I,J,K,2,1,iz] * dXdxLoc[I,J,K,3,3,iz] - 
      dXdxLoc[I,J,K,2,3,iz] * dXdxLoc[I,J,K,3,1,iz])
    @inbounds dXdxI[3,1,K,ID,Iz,IF] = dXdxLoc[I,J,K,2,1,iz] * dXdxLoc[I,J,K,3,2,iz] - 
      dXdxLoc[I,J,K,2,2,iz] * dXdxLoc[I,J,K,3,1,iz]

    @inbounds dXdxI[1,2,K,ID,Iz,IF] = -(dXdxLoc[I,J,K,1,2,iz] * dXdxLoc[I,J,K,3,3,iz] - 
      dXdxLoc[I,J,K,1,3,iz] * dXdxLoc[I,J,K,3,2,iz])
    @inbounds dXdxI[2,2,K,ID,Iz,IF] = dXdxLoc[I,J,K,1,1,iz] * dXdxLoc[I,J,K,3,3,iz] - 
      dXdxLoc[I,J,K,1,3,iz] * dXdxLoc[I,J,K,3,1,iz]
    @inbounds dXdxI[3,2,K,ID,Iz,IF] = -(dXdxLoc[I,J,K,1,1,iz] * dXdxLoc[I,J,K,3,2,iz] - 
      dXdxLoc[I,J,K,1,2,iz] * dXdxLoc[I,J,K,3,1,iz])

    @inbounds dXdxI[1,3,K,ID,Iz,IF] = dXdxLoc[I,J,K,1,2,iz] * dXdxLoc[I,J,K,2,3,iz] - 
      dXdxLoc[I,J,K,1,3,iz] * dXdxLoc[I,J,K,2,2,iz]
    @inbounds dXdxI[2,3,K,ID,Iz,IF] = -(dXdxLoc[I,J,K,1,1,iz] * dXdxLoc[I,J,K,2,3,iz] - 
      dXdxLoc[I,J,K,1,3,iz] * dXdxLoc[I,J,K,2,1,iz])
    @inbounds dXdxI[3,3,K,ID,Iz,IF] = dXdxLoc[I,J,K,1,1,iz] * dXdxLoc[I,J,K,2,2,iz] - 
      dXdxLoc[I,J,K,1,2,iz] * dXdxLoc[I,J,K,2,1,iz]
    dXdx[1,1,K,ID,Iz,IF] = dXdxLoc[I,J,K,1,1,iz] 
    dXdx[1,2,K,ID,Iz,IF] = dXdxLoc[I,J,K,1,2,iz] 
    dXdx[1,3,K,ID,Iz,IF] = dXdxLoc[I,J,K,1,3,iz] 
    dXdx[2,1,K,ID,Iz,IF] = dXdxLoc[I,J,K,2,1,iz] 
    dXdx[2,2,K,ID,Iz,IF] = dXdxLoc[I,J,K,2,2,iz] 
    dXdx[2,3,K,ID,Iz,IF] = dXdxLoc[I,J,K,2,3,iz] 
    dXdx[3,1,K,ID,Iz,IF] = dXdxLoc[I,J,K,3,1,iz] 
    dXdx[3,2,K,ID,Iz,IF] = dXdxLoc[I,J,K,3,2,iz] 
    dXdx[3,3,K,ID,Iz,IF] = dXdxLoc[I,J,K,3,3,iz] 
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


abstract type AdaptGrid end

Base.@kwdef struct GalChen <: AdaptGrid end

function (::GalChen)(H)
  function AdaptHeight(zRef,zs)
    z = zRef + (H - zRef) * zs / H
    DzDzRef  = eltype(zRef)(1) - zs / H
    return (z, DzDzRef)
  end
  return AdaptHeight
end

Base.@kwdef struct Sleve{T} <: AdaptGrid 
  etaH::T = .7
  s::T = 8/10

end

function (F::Sleve)(H)
  (;s,etaH) = F
  function AdaptHeight(zRef,zs)
    eta = zRef / H
    if eta <= etaH
      z = eta * H + zs * sinh((etaH - eta) / s / etaH) / sinh(1 / s) 
      DzDzRef  = eltype(zRef)(1) - zs / H / s / etaH * cosh((etaH - eta) / s / etaH) / sinh(1 / s) 
    else
      z = eta * H
      DzDzRef  = eltype(zRef)(1) 
    end  
    return (z, DzDzRef)
  end
  return AdaptHeight
end

function AdaptGrid(FT,Type,H)

  if Type == "GalChen"
    AdaptGridFunction = Grids.GalChen()(H)
  elseif Type == "Sleve"
    AdaptGridFunction = Grids.Sleve{FT}()(H)
  end
end
