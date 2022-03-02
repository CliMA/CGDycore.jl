function [X,D,J,theta]=JacobiSphere2(ksi1,ksi2,F,Grid)
Rad=Grid.Rad;

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
% if sqrt(X1^2+X2^2)<=1.e-12*Rad
%   X1=X1+1.e-6*Rad;
%   X2=X2+1.e-6*Rad;
% end
r=sqrt(X1^2+X2^2+X3^2);
f=Rad/r;
X1=X1/r;
X2=X2/r;
X3=X3/r;
[lam,theta]=cart2sphere(X1,X2,X3);


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
 
B=[Grid.Nodes(F.N(1)).P(1) Grid.Nodes(F.N(2)).P(1) Grid.Nodes(F.N(3)).P(1) Grid.Nodes(F.N(4)).P(1)
   Grid.Nodes(F.N(1)).P(2) Grid.Nodes(F.N(2)).P(2) Grid.Nodes(F.N(3)).P(2) Grid.Nodes(F.N(4)).P(2)
   Grid.Nodes(F.N(1)).P(3) Grid.Nodes(F.N(2)).P(3) Grid.Nodes(F.N(3)).P(3) Grid.Nodes(F.N(4)).P(3)];
C=0.25*[-1+ksi2  -1+ksi1
         1-ksi2  -1-ksi1
         1+ksi2   1+ksi1
        -1-ksi2   1-ksi1];
J=f*DD*A*B*C;
D=abs(det(J));
X=[X1 X2 X3]*Rad;
end

