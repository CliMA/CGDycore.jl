function JacobiSphere3GPU!(X,dXdxI,J,FE,F,z,zs,Rad)

  backend = get_backend(X)
  FT = eltype(X)

  NF = size(X,5)
  N = size(FE.xw,1)
  Nz = size(X,4)
 
  NzG = min(div(512,N*N*2),Nz)
  group = (N, N, 2, NzG, 1)
  ndrange = (N, N, 2, Nz, NF)

  KJacobiSphere3Kernel! = JacobiSphere3Kernel!(backend,group)

  H = z[end]
  KJacobiSphere3Kernel!(X,dXdxI,J,FE.xw,FE.xwZ,FE.DS,F,z,Rad,H,zs,ndrange=ndrange)
end

@kernel function JacobiSphere3Kernel!(X,dXdxI,JJ,@Const(ksi),@Const(zeta),@Const(D),
  @Const(F),@Const(z),Rad,H,@Const(zs))

  gi, gj, gk, gz, gF = @index(Group, NTuple)
  I, J, K, iz   = @index(Local, NTuple)
  _,_,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[4]
  N = @uniform @groupsize()[1]
  L = @uniform @groupsize()[3]
  Nz = @uniform @ndrange()[4]
  NF = @uniform @ndrange()[5]

  hR = @localmem eltype(X) (N,N,L,ColumnTilesDim)
  dXdx = @localmem eltype(X) (N,N,L,3,3,ColumnTilesDim)

  eta = ksi
  if Iz <= Nz
    ID = I + (J - 1) * N
    @inbounds z1 = z[Iz]
    @inbounds z2 = z[Iz+1]
    hRLoc = eltype(X)(0) 
    @views @inbounds JacobiSphere3Loc!(X[ID,K,:,Iz,IF],dXdx[I,J,K,:,:,iz],hRLoc,ksi[I],eta[J],zeta[K],F[:,:,IF],z1,z2,Rad,H,zs[I,J,IF])
    @inbounds hR[I,J,K,iz] = hRLoc
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

    @inbounds dXdx[I,J,K,3,1,iz] = DxhR
    @inbounds dXdx[I,J,K,3,2,iz] = DyhR
    @inbounds JJ[ID,K,Iz,IF] = det(reshape(dXdx[I,J,K,:,:,iz],3,3))
    @inbounds dXdxI[:,:,K,ID,Iz,IF] = inv(reshape(dXdx[I,J,K,:,:,iz],3,3))*JJ[ID,K,Iz,IF]
  end
end

@inline function JacobiSphere3Loc!(X,dXdx,hR,ksi1,ksi2,ksi3,F,z1,z2,Rad,H,zs)
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
  zLoc = half * ((one-ksi3) * z1 + (one+ksi3) * z2)
  hR = zLoc + (H - zLoc) * zs / H
  D33  = one - zs / H;
  D33 = half * D33*(z2-z1)

  r = sqrt(X1 * X1 + X2 * X2 + X3 * X3)
  f = Rad / r
  X1 = X1 / r
  X2 = X2 / r
  X3 = X3 / r
  (lam,theta)=cart2sphere(X1,X2,X3)

  DD=[-sin(lam) cos(lam) zero
      zero       zero     one]

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
  A = [a11 a12 a13
      a21 a22 a23
      a31 a32 a33]

  B = [F[1,1] F[2,1] F[3,1] F[4,1]
       F[1,2] F[2,2] F[3,2] F[4,2]
       F[1,3] F[2,3] F[3,3] F[4,3]]

  C = quarter * [-one+ksi2  -one+ksi1
              one-ksi2  -one-ksi1
              one+ksi2   one+ksi1
             -one-ksi2   one-ksi1]
  D = f * DD * A * B * C
  dXdx .= [D [zero; zero]
           zero zero D33]
  X[1] = X1 * (Rad + hR)
  X[2] = X2 * (Rad + hR)
  X[3] = X3 * (Rad + hR)

end

