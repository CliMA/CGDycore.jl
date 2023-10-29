function JacobiSphere3GPU!(X,dXdx,dXdxI,J,FE,F,z,zs)

  NF = size(X,4)
  N = size(FE.xw,1))
  Nz = size(X,2)
 
  NzG = min(div(512,N*N),Nz)
  group = (N, N, NzG, 1)
  ndrange = (N, N, Nz, NF)

  KJacobiSphere3Kernel! = JacobiSphere3Kernel!(backend,group)

  KJacobiSphere3Kernel!(X,dXdx,dXdxI,J,FE.xw,FE.DS,
end

@kernel function JacobiSphere3Kernel!(X,dXdx,dXdxI,J,ksi,zeta,DS,F,z,H,zs)

  gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  eta = ksi
  n3 = 2
  ID = I + (J - 1) * N
  @inbounds for k = 1 : n3
    JacobiSphere3Loc!(X[ID,k,:],dXdx[ID,k,:,:],ksi[I],eta[J],zeta[k],F,z,Topography.Rad,Topography,zs[I,J])
  end

  @inbounds for k=1:n3
    dXdx[I,:,k,3,1]=DS*hR[I,:,k]
    dXdx[I,:,k,3,2]=reshape(hR[I,:,k],N,N)*DS'
  end
  @inbounds for k=1:n3
    J[i,j,k]=det(reshape(dXdx[i,j,k,:,:],3,3))
    dXdxI[:,:,k,i,j]=inv(reshape(dXdx[i,j,k,:,:],3,3))*J[i,j,k]
  end

  X = reshape(X,n*n,n3,3)
  J = reshape(J,n*n,n3)
  dXdx = reshape(dXdx,n*n,n3,3,3)
  dXdxI = reshape(dXdxI,3,3,n3,n*n)
end

function JacobiSphere3Loc(X,dXdx,ksi1,ksi2,ksi3,F,z,Rad,H,zs)

  X1 = 0.25 * (F.P[1].x .* (1-ksi1)*(1-ksi2)+
    F.P[2].x .* (1+ksi1)* (1-ksi2)+
    F.P[3].x .* (1+ksi1)* (1+ksi2)+
    F.P[4].x .* (1-ksi1)* (1+ksi2))
  X2 = 0.25 * (F.P[1].y .*(1-ksi1)*(1-ksi2)+
    F.P[2].y .*(1+ksi1)*(1-ksi2)+
    F.P[3].y .*(1+ksi1)*(1+ksi2)+
    F.P[4].y .*(1-ksi1)*(1+ksi2))
  X3 = 0.25 * (F.P[1].z .*(1-ksi1)*(1-ksi2)+
    F.P[2].z .*(1+ksi1)*(1-ksi2)+
    F.P[3].z .*(1+ksi1)*(1+ksi2)+
    F.P[4].z .*(1-ksi1)*(1+ksi2))
  zLoc = 0.5 * ((1-ksi3) * z[1] + (1+ksi3) * z[2])
  hR = zLoc + (H - zLoc) * zs / H
  D33  = 1 - zs / Topography.H;
  D33 = 0.5 * D33*(z[2]-z[1])

  r = sqrt(X1^2 + X2^2 + X3^2)
  f = Rad / r
  X1 = X1 / r
  X2 = X2 / r
  X3 = X3 / r
  (lam,theta)=cart2sphere(X1,X2,X3)

  DD=[-sin(lam) cos(lam) 0
      0          0       1]

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

  B = [F.P[1].x F.P[2].x F.P[3].x F.P[4].x
      F.P[1].y F.P[2].y F.P[3].y F.P[4].y
      F.P[1].z F.P[2].z F.P[3].z F.P[4].z]

  C =0.25 * [-1+ksi2  -1+ksi1
              1-ksi2  -1-ksi1
              1+ksi2   1+ksi1
             -1-ksi2   1-ksi1]
  D = f * DD * A * B * C
  dXdx .= [D [0; 0]
           0 0 D33]
  X .= [X1 X2 X3]*(Rad+hR)

end

