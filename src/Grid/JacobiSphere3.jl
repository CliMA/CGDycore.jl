function JacobiSphere3(DG,F,z,Topo,Topography)
ksi=DG.xw;
eta=DG.xw;
zeta=DG.xwZ;
n=DG.OrdPoly+1;
nz=DG.OrdPolyZ+1;
X=zeros(n,n,nz,3);
dXdx=zeros(n,n,nz,3,3);
dXdxI=zeros(n,n,nz,3,3);
J=zeros(n,n,nz);
theta=zeros(n,n);
for j=1:n
  for i=1:n
    for k=1:nz
      (X[i,j,k,:],J[i,j,k],dXdx[i,j,k,:,:],dXdxI[i,j,k,:,:],theta[i,j]) =
        JacobiSphere3Loc(ksi[i],eta[j],zeta[k],F,z,Topography.Rad);
    end
  end
end
return (X,J,dXdx,dXdxI,theta)
end

function JacobiSphere3Loc(ksi1,ksi2,ksi3,F,z,Rad)

X1=0.25*(F.P[1].x .*(1-ksi1)*(1-ksi2)+
  F.P[2].x .*(1+ksi1)*(1-ksi2)+
  F.P[3].x .*(1+ksi1)*(1+ksi2)+
  F.P[4].x .*(1-ksi1)*(1+ksi2));
X2=0.25*(F.P[1].y .*(1-ksi1)*(1-ksi2)+
  F.P[2].y .*(1+ksi1)*(1-ksi2)+
  F.P[3].y .*(1+ksi1)*(1+ksi2)+
  F.P[4].y .*(1-ksi1)*(1+ksi2));
X3=0.25*(F.P[1].z .*(1-ksi1)*(1-ksi2)+
  F.P[2].z .*(1+ksi1)*(1-ksi2)+
  F.P[3].z .*(1+ksi1)*(1+ksi2)+
  F.P[4].z .*(1-ksi1)*(1+ksi2));
zLoc=0.5*((1-ksi3)*z[1]+(1+ksi3)*z[2]);

r=sqrt(X1^2+X2^2+X3^2);
f=Rad/r;
X1=X1/r;
X2=X2/r;
X3=X3/r;
(lam,theta)=cart2sphere(X1,X2,X3);


DD=[-sin(lam) cos(lam) 0
  0        0     1];

sinlam=sin(lam);
coslam=cos(lam);
sinth=sin(theta);
costh=cos(theta);
a11=sinlam*sinlam*costh*costh+sinth*sinth;
a12=-sinlam*coslam*costh*costh;
a13=-coslam*sinth*costh;
a21=a12;
a22=coslam*coslam*costh*costh+sinth*sinth;
a23=-sinlam*sinth*costh;
a31=-coslam*sinth;
a32=-sinlam*sinth;
a33=costh;
A=[a11 a12 a13
  a21 a22 a23
  a31 a32 a33];

#B=[Grid.Nodes[F.N[1]].P[1] Grid.Nodes[F.N[2]].P[1] Grid.Nodes[F.N[3]].P[1] Grid.Nodes[F.N[4]].P[1]
#  Grid.Nodes[F.N[1]].P[2] Grid.Nodes[F.N[2]].P[2] Grid.Nodes[F.N[3]].P[2] Grid.Nodes[F.N[4]].P[2]
#  Grid.Nodes[F.N[1]].P[3] Grid.Nodes[F.N[2]].P[3] Grid.Nodes[F.N[3]].P[3] Grid.Nodes[F.N[4]].P[3]];
B=[F.P[1].x F.P[2].x F.P[3].x F.P[4].x
   F.P[1].y F.P[2].y F.P[3].y F.P[4].y
   F.P[1].z F.P[2].z F.P[3].z F.P[4].z]

C=0.25*[-1+ksi2  -1+ksi1
  1-ksi2  -1-ksi1
  1+ksi2   1+ksi1
  -1-ksi2   1-ksi1];
D=f*DD*A*B*C;
D=[D [0;0]
  0 0 0.5*(z[2]-z[1])];
J=abs(det(D));
X=[X1 X2 X3]*(Rad+zLoc);

DI=inv(D)*J;

return (X,J,D,DI,theta)

end

