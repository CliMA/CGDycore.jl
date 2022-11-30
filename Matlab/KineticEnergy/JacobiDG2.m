function [X,J,dXdx,dXdxI]=JacobiDG2(xP,zP,Dx,x,Dz,z)
ksi=x;
eta=z;
nX=size(ksi,1);
nY=size(eta,1);
X=zeros(nX,nY,2);
dXdx=zeros(nX,nY,2,2);
dXdxI=zeros(nX,nY,2,2);
J=zeros(nX,nY);
for j=1:nY
  for i=1:nX
    X(i,j,1:2)=[xP(i),zP(i,j)];
  end
end

dXdx(:,:,1,1)=Dx*X(:,:,1);
dXdx(:,:,2,1)=Dx*X(:,:,2);
dXdx(:,:,1,2)=X(:,:,1)*Dz';
dXdx(:,:,2,2)=X(:,:,2)*Dz';

for i=1:nX
  for j=1:nY
    J(i,j)=det(reshape(dXdx(i,j,:,:),2,2));
    dXdxI(i,j,:,:)=inv(reshape(dXdx(i,j,:,:),2,2))*J(i,j);
  end
end

end

