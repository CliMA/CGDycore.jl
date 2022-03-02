function TestJacobi()
dlambda=pi/8;
dtheta=pi/8;
s1=[0 0];
s2=[dlambda 0];
s3=[dlambda dtheta];
s4=[0 dtheta];

c1=SphereToCart(s1(1),s1(2));
c2=SphereToCart(s2(1),s2(2));
c3=SphereToCart(s3(1),s3(2));
c4=SphereToCart(s4(1),s4(2));

xi1=.3;
xi2=.4;
cTilde=rTilde(xi1,xi2,c1,c2,c3,c4);
c=cTilde/norm(cTilde);
[lambda,theta]=CartToSphereToCart(c);

cNeu=SphereToCart(lambda,theta);

DV1=0.25/norm(cTilde)*JacS1(lambda,theta)*JacS2(lambda,theta);

D0=[1/norm(cTilde) 0
   0 norm(cTilde)^2];
 
 DV2=0.25*D0*JacC1(cTilde)*JacC2(cTilde);
 
 DV3=0.25/norm(cTilde)^(3/2)*JacCC1(c)*JacCC2(cTilde);
 ee=dd();
end

function x=SphereToCart(lambda,theta)
x=zeros(3,1);
x(1)=cos(lambda)*cos(theta);
x(2)=sin(lambda)*cos(theta);
x(3)=sin(theta);
end

function [lambda,theta]=CartToSphereToCart(x)
r=norm(x);
lambda=atan2(x(2),x(1));
if x(2)<0
  lambda=lambda+2*pi;
end
theta=asin(x(3)/r);
end

function x=rTilde(xi1,xi2,c1,c2,c3,c4)
x=0.25*((1-xi1)*(1-xi2)*c1...
  +(1+xi1)*(1-xi2)*c2...
  +(1+xi1)*(1+xi2)*c3...
  +(1-xi1)*(1+xi2)*c4);
end

function D=JacCC1(x)
D=zeros(2,2);
D(1,1)=-x(2)/(x(1)*x(1)+x(2)*x(2));
D(1,2)=x(1)/(x(1)*x(1)+x(2)*x(2));

D(2,1)=x(1)/((x(1)^2 + x(2)^2)^(1/2)*(- x(1)^2 - x(2)^2 + 1)^(1/2));
D(2,2)=x(2)/((x(1)^2 + x(2)^2)^(1/2)*(- x(1)^2 - x(2)^2 + 1)^(1/2));
end

function D=JacCC2(x)
D=zeros(2,3);
D(1,1)=x(2)*x(2)+x(3)*x(3);
D(1,2)=-x(1)*x(2);
D(1,3)=-x(1)*x(3);

D(2,2)=-x(1)*x(2);
D(2,2)=x(1)*x(1)+x(3)*x(3);
D(2,3)=-x(2)*x(3);


end


function D=JacC1(x)
D=zeros(2,3);
D(1,1)=-x(2)/(x(1)*x(1)+x(2)*x(2));
D(1,2)=x(1)/(x(1)*x(1)+x(2)*x(2));

D(2,1)=x(1)*x(3)/sqrt(x(1)*x(1)+x(2)*x(2));
D(2,2)=x(2)*x(3)/sqrt(x(1)*x(1)+x(2)*x(2));
D(2,3)=-sqrt(x(1)*x(1)+x(2)*x(2));
end

function D=JacC2(x)
D=zeros(3,3);
D(1,1)=x(2)*x(2)+x(3)*x(3);
D(1,2)=-x(1)*x(2);
D(1,3)=-x(1)*x(3);

D(2,2)=-x(1)*x(2);
D(2,2)=x(1)*x(1)+x(3)*x(3);
D(2,3)=-x(2)*x(3);

D(3,1)=-x(1)*x(3);
D(3,2)=-x(2)*x(3);
D(3,3)=x(1)*x(1)+x(2)*x(2);
end

function D=JacS1(lambda,theta)
D=zeros(2,3);
D(1,1)=-sin(lambda);
D(1,2)=cos(lambda);

D(2,3)=1;
end

function D=JacS2(lambda,theta)
D=zeros(3,3);
D(1,1)=sin(lambda)^2*cos(theta)^2+sin(theta)^2;
D(1,2)=-0.5*sin(2*lambda)*cos(theta)^2;
D(1,3)=0.5*cos(lambda)*sin(2*theta);

D(2,1)=-0.5*sin(2*lambda)*cos(theta)^2;
D(2,2)=cos(lambda)^2*cos(theta)^2+sin(theta)^2;
D(2,3)=-0.5*sin(lambda)*sin(2*theta);

D(3,1)=-0.5*cos(lambda)*sin(theta);
D(3,2)=-0.5*sin(lambda)*sin(theta);
D(3,3)=cos(theta);
end

function erg=dd()
syms x
syms y
f=acos(sqrt(1-x*x-y*y));
erg=diff(f,x);
end



