function [X,J,dXdx,dXdxI]=JacobiDG2(DG,F,Topo,Param)
ksi=DG.xwX;
eta=DG.xwY;
nX=DG.OrdPolyX+1;
nY=DG.OrdPolyY+1;
X=zeros(nX,nY,3);
dXdx=zeros(nX,nY,2,2);
dXdxI=zeros(nX,nY,2,2);
J=zeros(nX,nY);
for j=1:nY
  for i=1:nX
    X(i,j,1:2)=0.25*((1-ksi(i))*(1-eta(j))*F.P(1:2,1)...
                    +(1+ksi(i))*(1-eta(j))*F.P(1:2,2)...
                    +(1+ksi(i))*(1+eta(j))*F.P(1:2,3)...
                    +(1-ksi(i))*(1+eta(j))*F.P(1:2,4));
    [X(i,j,2),dXdx(i,j,2,2)]=Topo(X(i,j,1),X(i,j,2),Param);
  end
end

dXdx(:,:,1,1)=DG.DSX*X(:,:,1);
dXdx(:,:,2,1)=DG.DSX*X(:,:,2);
dXdx(:,:,1,2)=X(:,:,1)*DG.DSY';
dXdx(:,:,2,2)=X(:,:,2)*DG.DSY';
for j=1:nY
  for i=1:nX
    J(i,j)=det(reshape(dXdx(i,j,:,:),2,2));
    dXdxI(i,j,:,:)=inv(reshape(dXdx(i,j,:,:),2,2))*J(i,j);
  end
end
end

