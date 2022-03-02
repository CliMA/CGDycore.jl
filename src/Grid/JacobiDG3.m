function [X,J,dXdx,dXdxI,theta]=JacobiDG3(DG,F,z,Topo,Param)
ksi=DG.xw;
eta=DG.xw;
zeta=DG.xwZ;
z1=z(1);
z2=z(2);
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
    X(i,j,k,1:2)=0.25*((1-ksi(i))*(1-eta(j))*F.P(1:2,1)...
                    +(1+ksi(i))*(1-eta(j))*F.P(1:2,2)...
                    +(1+ksi(i))*(1+eta(j))*F.P(1:2,3)...
                    +(1-ksi(i))*(1+eta(j))*F.P(1:2,4));
    z=0.5*((1-zeta(k))*z1+(1+zeta(k))*z2);              
    [X(i,j,k,3),dXdx(i,j,k,3,3)]=Topo(X(i,j,1),X(i,j,2),z,Param);
    dXdx(i,j,k,3,3)=0.5*dXdx(i,j,k,3,3)*(z2-z1);
    end
  end
end

for k=1:nz
  dXdx(:,:,k,1,1)=DG.DS*X(:,:,k,1);
  dXdx(:,:,k,2,1)=DG.DS*X(:,:,k,2);
  dXdx(:,:,k,3,1)=DG.DS*X(:,:,k,3);
  dXdx(:,:,k,1,2)=reshape(X(:,:,k,1),n,n)*DG.DS';
  dXdx(:,:,k,2,2)=reshape(X(:,:,k,2),n,n)*DG.DS';
  dXdx(:,:,k,3,2)=reshape(X(:,:,k,3),n,n)*DG.DS';
end
for j=1:n
  for i=1:n
    for k=1:nz
      J(i,j,k)=det(reshape(dXdx(i,j,k,:,:),3,3));
      dXdxI(i,j,k,:,:)=inv(reshape(dXdx(i,j,k,:,:),3,3))*J(i,j,k);
    end
  end
end
end

