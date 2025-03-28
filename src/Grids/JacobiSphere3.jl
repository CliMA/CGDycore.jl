function JacobiSphere3CPU!(X,dXdxI,J,CG,Grid,zs)
  FT = eltype(X)

  NF = size(X,5)
  N = size(CG.xw,1)
  nz = size(X,4)
  XCPU = zeros(FT,size(X))
  JCPU = zeros(FT,size(J))
  dXdxICPU = zeros(FT,size(dXdxI))
  zCPU = zeros(FT,size(Grid.z))
  zsCPU = zeros(FT,size(zs))
  copyto!(zCPU,Grid.z)
  copyto!(zsCPU,zs)
  H = zCPU[end]
  Rad = Grid.Rad

  for iF = 1 : NF
    for iz = 1 : nz
      zI = [zCPU[iz],zCPU[iz+1]]
      @views (X_Fz,J_Fz,dXdx_Fz,dXdxI_Fz) = JacobiSphere3(CG,Grid.Faces[iF],zI,Rad,H,zsCPU[:,:,iF],iz,iF)
      @views @. XCPU[:,:,:,iz,iF] = X_Fz
      @views @. JCPU[:,:,iz,iF] = J_Fz
      @views @. dXdxICPU[:,:,:,:,iz,iF] = dXdxI_Fz
    end
  end
  copyto!(X,XCPU)
  copyto!(J,JCPU)
  copyto!(dXdxI,dXdxICPU)
end

function JacobiSphere3(CG,F,z,Rad,H,zs,iz,iF)
ksi=CG.xwCPU
eta=CG.xwCPU
zeta=CG.xwZCPU
n=CG.OrdPoly+1
n3=CG.OrdPolyZ+1
X=zeros(n,n,n3,3)
dXdx=zeros(n,n,n3,3,3)
dXdxI=zeros(3,3,n3,n,n)
J=zeros(n,n,n3)
hR=zeros(n,n,n3)
theta=zeros(n,n)
(_,DS)=DG.DerivativeMatrixSingle(CG.OrdPoly)
DST = DS'
(_,DSZ)=DG.DerivativeMatrixSingle(CG.OrdPolyZ)
@inbounds for j=1:n
  @inbounds for i=1:n
    @inbounds for k=1:n3
      (X[i,j,k,:],dXdx[i,j,k,:,:],hR[i,j,k]) =
        JacobiSphere3Loc(ksi[i],eta[j],zeta[k],F,z,Rad,H,zs[i,j])
    end
  end
end
@inbounds for k=1:n3
  dXdx[:,:,k,3,1]=DS*hR[:,:,k]
  dXdx[:,:,k,3,2]=reshape(hR[:,:,k],n,n)*DST
end
@inbounds for j=1:n
  @inbounds for i=1:n
    @inbounds for k=1:n3
      J[i,j,k]=det(reshape(dXdx[i,j,k,:,:],3,3))
      dXdxI[:,:,k,i,j]=inv(reshape(dXdx[i,j,k,:,:],3,3))*J[i,j,k]
    end
  end
end

X = reshape(X,n*n,n3,3)
J = reshape(J,n*n,n3)
dXdx = reshape(dXdx,n*n,n3,3,3)
dXdxI = reshape(dXdxI,3,3,n3,n*n)
return (X,J,dXdx,dXdxI)
end

function JacobiSphere3Loc(ksi1,ksi2,ksi3,F,z,Rad,H,zs)

X1=0.25*(F.P[1].x .*(1-ksi1)*(1-ksi2)+
  F.P[2].x .*(1+ksi1)*(1-ksi2)+
  F.P[3].x .*(1+ksi1)*(1+ksi2)+
  F.P[4].x .*(1-ksi1)*(1+ksi2))
X2=0.25*(F.P[1].y .*(1-ksi1)*(1-ksi2)+
  F.P[2].y .*(1+ksi1)*(1-ksi2)+
  F.P[3].y .*(1+ksi1)*(1+ksi2)+
  F.P[4].y .*(1-ksi1)*(1+ksi2))
X3=0.25*(F.P[1].z .*(1-ksi1)*(1-ksi2)+
  F.P[2].z .*(1+ksi1)*(1-ksi2)+
  F.P[3].z .*(1+ksi1)*(1+ksi2)+
  F.P[4].z .*(1-ksi1)*(1+ksi2))
zLoc=0.5*((1-ksi3)*z[1]+(1+ksi3)*z[2])
(hR,D33)=Topo1(X1,X2,X3,zLoc,H,zs)
D33=0.5*D33*(z[2]-z[1])

r=sqrt(X1^2+X2^2+X3^2)
f=Rad/r
X1=X1/r
X2=X2/r
X3=X3/r
(lam,theta)=cart2sphere(X1,X2,X3)


DD=[-sin(lam) cos(lam) 0
  0        0     1]

sinlam=sin(lam)
coslam=cos(lam)
sinth=sin(theta)
costh=cos(theta)
a11=sinlam*sinlam*costh*costh+sinth*sinth
a12=-sinlam*coslam*costh*costh
a13=-coslam*sinth*costh
a21=a12
a22=coslam*coslam*costh*costh+sinth*sinth
a23=-sinlam*sinth*costh
a31=-coslam*sinth
a32=-sinlam*sinth
a33=costh
A=[a11 a12 a13
  a21 a22 a23
  a31 a32 a33]

#B=[Grid.Nodes[F.N[1]].P[1] Grid.Nodes[F.N[2]].P[1] Grid.Nodes[F.N[3]].P[1] Grid.Nodes[F.N[4]].P[1]
#  Grid.Nodes[F.N[1]].P[2] Grid.Nodes[F.N[2]].P[2] Grid.Nodes[F.N[3]].P[2] Grid.Nodes[F.N[4]].P[2]
#  Grid.Nodes[F.N[1]].P[3] Grid.Nodes[F.N[2]].P[3] Grid.Nodes[F.N[3]].P[3] Grid.Nodes[F.N[4]].P[3]]
B=[F.P[1].x F.P[2].x F.P[3].x F.P[4].x
   F.P[1].y F.P[2].y F.P[3].y F.P[4].y
   F.P[1].z F.P[2].z F.P[3].z F.P[4].z]

C=0.25*[-1+ksi2  -1+ksi1
  1-ksi2  -1-ksi1
  1+ksi2   1+ksi1
  -1-ksi2   1-ksi1]
D=f*DD*A*B*C
D=[D [0; 0]
  0 0 D33]
X=[X1 X2 X3]*(Rad+hR)


return (X,D,hR)

end

function Topo1(x,y,z,zeta,H,zs)
  h = zs
  Z=zeta+(H-zeta)*h/H;
  dZ=1-h/H;
 return (Z,dZ)
end
