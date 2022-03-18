function [X,J,dXdx,theta]=JacobiSphereDG(DG,F,Grid)
Rad=Grid.Rad;
ksi=DG.xw;
eta=DG.xw;
n=DG.OrdPoly+1;
lam=zeros(n,n);
theta=zeros(n,n);
X=zeros(n,n,3);
for j=1:n
  ksi2=eta(j);
  for i=1:n
    ksi1=ksi(i);
    X1=0.25*(F.P(1,1).*(1-ksi1)*(1-ksi2)+...
      F.P(1,2).*(1+ksi1)*(1-ksi2)+...
      F.P(1,3).*(1+ksi1)*(1+ksi2)+...
      F.P(1,4).*(1-ksi1)*(1+ksi2));
    X2=0.25*(F.P(2,1).*(1-ksi1)*(1-ksi2)+...
      F.P(2,2).*(1+ksi1)*(1-ksi2)+...
      F.P(2,3).*(1+ksi1)*(1+ksi2)+...
      F.P(2,4).*(1-ksi1)*(1+ksi2));
    X3=0.25*(F.P(3,1).*(1-ksi1)*(1-ksi2)+...
      F.P(3,2).*(1+ksi1)*(1-ksi2)+...
      F.P(3,3).*(1+ksi1)*(1+ksi2)+...
      F.P(3,4).*(1-ksi1)*(1+ksi2));
    X(i,j,1)=X1;
    X(i,j,2)=X2;
    X(i,j,3)=X3;
    r=sqrt(X1^2+X2^2+X3^2);
    f=Rad/r;
    X1=X1/r;
    X2=X2/r;
    X3=X3/r;
    [lam(i,j),theta(i,j)]=cart2sphere(X1,X2,X3);
    
    
    % DD=[-sin(lam) cos(lam) 0
    %   0        0     1];
    %
    % sinlam=sin(lam);
    % coslam=cos(lam);
    % sinth=sin(theta);
    % costh=cos(theta);
    % a11=sinlam*sinlam*costh*costh+sinth*sinth;
    % a12=-sinlam*coslam*costh*costh;
    % a13=-coslam*sinth*costh;
    % a21=a12;
    % a22=coslam*coslam*costh*costh+sinth*sinth;
    % a23=-sinlam*sinth*costh;
    % a31=-coslam*sinth;
    % a32=-sinlam*sinth;
    % a33=costh;
    % A=[a11 a12 a13
    %    a21 a22 a23
    %    a31 a32 a33];
    %
    % B=[Grid.Nodes(F.N(1)).P(1) Grid.Nodes(F.N(2)).P(1) Grid.Nodes(F.N(3)).P(1) Grid.Nodes(F.N(4)).P(1)
    %    Grid.Nodes(F.N(1)).P(2) Grid.Nodes(F.N(2)).P(2) Grid.Nodes(F.N(3)).P(2) Grid.Nodes(F.N(4)).P(2)
    %    Grid.Nodes(F.N(1)).P(3) Grid.Nodes(F.N(2)).P(3) Grid.Nodes(F.N(3)).P(3) Grid.Nodes(F.N(4)).P(3)];
    % C=0.25*[-1+ksi2  -1+ksi1
    %          1-ksi2  -1-ksi1
    %          1+ksi2   1+ksi1
    %         -1-ksi2   1-ksi1];
    % J=f*DD*A*B*C;
    % D=abs(det(J));
    % X=[X1 X2 X3]*Rad;
  end
end
minlam=min(min(lam));
maxlam=max(max(lam));
if maxlam-minlam>pi
  for i=1:n
    for j=1:n
      if lam(i,j)>pi
        lam(i,j)=lam(i,j)-pi;
      else
        lam(i,j)=lam(i,j)+pi;
      end
    end
  end
end
lam=Rad*lam; %.*cos(theta);
costheta=cos(theta);
theta=theta*Rad;
dXdx=zeros(n,n,2,2);
dXdx(:,:,1,1)=(DG.DS*lam).*costheta;
dXdx(:,:,2,1)=DG.DS*theta;
dXdx(:,:,1,2)=(lam*DG.DS').*costheta;
dXdx(:,:,2,2)=theta*DG.DS';
J=dXdx(:,:,1,1).*dXdx(:,:,2,2)-dXdx(:,:,1,2).*dXdx(:,:,2,1);

end

