function JacobiSphere3(CG,F,z,Topo,Topography,zs)
ksi=CG.xw
eta=CG.xw
zeta=CG.xwZ
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
        JacobiSphere3Loc(ksi[i],eta[j],zeta[k],F,z,Topography.Rad,Topography,zs[i,j])
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

function JacobiSphere3Loc(ksi1,ksi2,ksi3,F,z,Rad,Topography,zs)

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
(hR,D33)=Topo(X1,X2,X3,zLoc,Topography,zs)
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

