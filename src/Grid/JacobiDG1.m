function [X,J,dXdx]=JacobiDG1(DG,F,Topo,Param)
ksi=DG.xwX;
eta=DG.xwY;
nX=DG.OrdPolyX+1;
nY=DG.OrdPolyY+1;
XE1=Topo(DG,F.P(:,1),F.P(:,2),Param.Grid.Edges(F.E(1)).Type,Param);
Topo(X(i,j,1),X(i,j,2),z,Param)
XE2=Topo(DG,F.P(:,2),F.P(:,3),Param.Grid.Edges(F.E(2)).Type,Param);
XE3=Topo(DG,F.P(:,4),F.P(:,3),Param.Grid.Edges(F.E(3)).Type,Param);
XE4=Topo(DG,F.P(:,1),F.P(:,4),Param.Grid.Edges(F.E(4)).Type,Param);
X=zeros(nX,nY,3);
for j=1:nY
  for i=1:nX
    X(i,j,1:2)=0.5*((1-ksi(i))*XE4(j,:)+(1+ksi(i))*XE2(j,:)+(1-eta(j))*XE1(i,:)+(1+eta(j))*XE3(i,:))-...
      0.25*((1-ksi(i))*((1-eta(j))*F.P(1:2,1)'+(1+eta(j))*F.P(1:2,4)')+...
      (1+ksi(i))*((1-eta(j))'*F.P(1:2,2)'+(1+eta(j))*F.P(1:2,3)'));
  end
end
dXdx=zeros(nX,nY,2,2);
dXdx(:,:,1,1)=DG.DSX*X(:,:,1);
dXdx(:,:,2,1)=DG.DSX*X(:,:,2);
dXdx(:,:,1,2)=X(:,:,1)*(DG.DSY');
dXdx(:,:,2,2)=X(:,:,2)*(DG.DSY');
J=dXdx(:,:,1,1).*dXdx(:,:,2,2)-dXdx(:,:,1,2).*dXdx(:,:,2,1);
end

