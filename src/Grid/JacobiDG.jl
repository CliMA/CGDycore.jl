function [J,dXdx,X]=JacobiDG(DG,F,Grid,Trans)
ksi=DG.xw;
eta=DG.xw;
n=DG.OrdPoly+1;
XE1=Trans(DG,Grid.Nodes(F.N(1)).P,Grid.Nodes(F.N(2)).P);
XE2=Trans(DG,Grid.Nodes(F.N(2)).P,Grid.Nodes(F.N(3)).P);
XE3=Trans(DG,Grid.Nodes(F.N(4)).P,Grid.Nodes(F.N(3)).P);
XE4=Trans(DG,Grid.Nodes(F.N(1)).P,Grid.Nodes(F.N(4)).P);
X=zeros(n,n,2);
for j=1:n
  for i=1:n
    X(i,j,:)=0.5*((1-ksi(i))*XE4(j,:)+(1+ksi(i))*XE2(j,:)+(1-eta(j))*XE1(i,:)+(1+eta(j))*XE3(i,:))-...
      0.25*((1-ksi(i))*((1-eta(j))*Grid.Nodes(F.N(1)).P(1:2)+(1+eta(j))*Grid.Nodes(F.N(4)).P(1:2))+...
      (1+ksi(i))*((1-eta(j))'*Grid.Nodes(F.N(2)).P(1:2)+(1+eta(j))*Grid.Nodes(F.N(3)).P(1:2)));
  end
end
dXdx=zeros(n,n,2,2);
dXdx(:,:,1,1)=DG.DS*X(:,:,1);
dXdx(:,:,2,1)=DG.DS*X(:,:,2);
dXdx(:,:,1,2)=X(:,:,1)*DG.DS';
dXdx(:,:,2,2)=X(:,:,2)*DG.DS';
J=dXdx(:,:,1,1).*dXdx(:,:,2,2)-dXdx(:,:,1,2).*dXdx(:,:,2,1);
end

