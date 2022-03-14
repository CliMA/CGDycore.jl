function [D,J,JInv,X]=JacobiSphere(ksi,F,Grid)
Rad=Grid.Rad;
X1=Grid.Nodes(F.N(1)).P(1)...
  +(Grid.Nodes(F.N(2)).P(1)-Grid.Nodes(F.N(1)).P(1))*ksi(1)...
  +(Grid.Nodes(F.N(4)).P(1)-Grid.Nodes(F.N(1)).P(1))*ksi(2)...
  +(Grid.Nodes(F.N(3)).P(1)-Grid.Nodes(F.N(4)).P(1)-Grid.Nodes(F.N(2)).P(1)+Grid.Nodes(F.N(1)).P(1))*ksi(1)*ksi(2);
X2=Grid.Nodes(F.N(1)).P(2)...
  +(Grid.Nodes(F.N(2)).P(2)-Grid.Nodes(F.N(1)).P(2))*ksi(1)...
  +(Grid.Nodes(F.N(4)).P(2)-Grid.Nodes(F.N(1)).P(2))*ksi(2)...
  +(Grid.Nodes(F.N(3)).P(2)-Grid.Nodes(F.N(4)).P(2)-Grid.Nodes(F.N(2)).P(2)+Grid.Nodes(F.N(1)).P(2))*ksi(1)*ksi(2);
X3=Grid.Nodes(F.N(1)).P(3)...
  +(Grid.Nodes(F.N(2)).P(3)-Grid.Nodes(F.N(1)).P(3))*ksi(1)...
  +(Grid.Nodes(F.N(4)).P(3)-Grid.Nodes(F.N(1)).P(3))*ksi(2)...
  +(Grid.Nodes(F.N(3)).P(3)-Grid.Nodes(F.N(4)).P(3)-Grid.Nodes(F.N(2)).P(3)+Grid.Nodes(F.N(1)).P(3))*ksi(1)*ksi(2);
f=Rad*(X1^2+X2^2+X3^2)^(-3/2);

dx1dX1=f*(X2^2+X3^2);
dx1dX2=-f*X1*X2; 
dx1dX3=-f*X1*X3;
dx2dX1=dx1dX2; 
dx2dX2=f*(X1^2+X3^2); 
dx2dX3=-f*X2*X3;
dx3dX1=dx1dX3; 
dx3dX2=dx2dX3; 
dx3dX3=f*(X1^2+X2^2);

% dx1dX1=1;
% dx1dX2=-0; 
% dx1dX3=0;
% dx2dX1=0; 
% dx2dX2=1; 
% dx2dX3=0;
% dx3dX1=0; 
% dx3dX2=0; 
% dx3dX3=1;

dX1dksi1=(Grid.Nodes(F.N(2)).P(1)-Grid.Nodes(F.N(1)).P(1))...
  +(Grid.Nodes(F.N(3)).P(1)-Grid.Nodes(F.N(4)).P(1)...
  -Grid.Nodes(F.N(2)).P(1)+Grid.Nodes(F.N(1)).P(1))*ksi(2);
dX1dksi2=(Grid.Nodes(F.N(4)).P(1)-Grid.Nodes(F.N(1)).P(1))...
  +(Grid.Nodes(F.N(3)).P(1)-Grid.Nodes(F.N(4)).P(1)...
  -Grid.Nodes(F.N(2)).P(1)+Grid.Nodes(F.N(1)).P(1))*ksi(1);
dX2dksi1=(Grid.Nodes(F.N(2)).P(2)-Grid.Nodes(F.N(1)).P(2))...
  +(Grid.Nodes(F.N(3)).P(2)-Grid.Nodes(F.N(4)).P(2)...
  -Grid.Nodes(F.N(2)).P(2)+Grid.Nodes(F.N(1)).P(2))*ksi(2);
dX2dksi2=(Grid.Nodes(F.N(4)).P(2)-Grid.Nodes(F.N(1)).P(2))...
  +(Grid.Nodes(F.N(3)).P(2)-Grid.Nodes(F.N(4)).P(2)...
  -Grid.Nodes(F.N(2)).P(2)+Grid.Nodes(F.N(1)).P(2))*ksi(1);
dX3dksi1=(Grid.Nodes(F.N(2)).P(3)-Grid.Nodes(F.N(1)).P(3))...
  +(Grid.Nodes(F.N(3)).P(3)-Grid.Nodes(F.N(4)).P(3)...
  -Grid.Nodes(F.N(2)).P(3)+Grid.Nodes(F.N(1)).P(3))*ksi(2);
dX3dksi2=(Grid.Nodes(F.N(4)).P(3)-Grid.Nodes(F.N(1)).P(3))...
  +(Grid.Nodes(F.N(3)).P(3)-Grid.Nodes(F.N(4)).P(3)...
  -Grid.Nodes(F.N(2)).P(3)+Grid.Nodes(F.N(1)).P(3))*ksi(1);

J=zeros(3,2);

J(1,1)=dx1dX1*dX1dksi1+dx1dX2*dX2dksi1+dx1dX3*dX3dksi1;
J(2,1)=dx2dX1*dX1dksi1+dx2dX2*dX2dksi1+dx2dX3*dX3dksi1;
J(3,1)=dx3dX1*dX1dksi1+dx3dX2*dX2dksi1+dx3dX3*dX3dksi1;
J(1,2)=dx1dX1*dX1dksi2+dx1dX2*dX2dksi2+dx1dX3*dX3dksi2;
J(2,2)=dx2dX1*dX1dksi2+dx2dX2*dX2dksi2+dx2dX3*dX3dksi2;
J(3,2)=dx3dX1*dX1dksi2+dx3dX2*dX2dksi2+dx3dX3*dX3dksi2;

D=norm(cross(J(:,1),J(:,2)),2);
if nargout > 2
  X=[X1 X2 X3]*(Rad/sqrt((X1^2+X2^2+X3^2)));   X=[X1 X2 X3];
end
JInv=inv(J'*J)*J'*sqrt(det(J'*J));
end

