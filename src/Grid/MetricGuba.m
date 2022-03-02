syms lon lat r
syms x1 y1 z1
syms xx1 yy1 zz1
syms lam1 theta1 r1

% xx1=x1/sqrt(x1^2+y1^2+z1^2);
% yy1=y1/sqrt(x1^2+y1^2+z1^2);
% zz1=z1/sqrt(x1^2+y1^2+z1^2);

% xx1=r*cos(lam)*cos(theta)
% yy1=r*sin(lam)*cos(theta)
% zz1=r*sin(theta);
% yy1/xx1=sin(lam)/cos(lam)
% lam=atan(yy1/xx1)
% theta=asin(zz1/sqrt(xx1^2+yy1^2+zz1^2))
% r=sqrt(xx1^2+yy1^2+zz1^2);
%
% d atan(s) / ds = 1/(1+s^2)
% d asin(s) / ds = 1/sqrt(1-s^2)


lam=atan(yy1/xx1);
theta=asin(zz1);

D1=[diff(lam,xx1) diff(lam,yy1) diff(lam,zz1)
    diff(theta,xx1) diff(theta,yy1) diff(theta,zz1)];
D11=subs(D1,{xx1,yy1,zz1},{cos(lam1)*cos(theta1),sin(lam1)*cos(theta1),sin(theta1)}); 
D111=simplify(D11);

C=[cos(theta1) 0 
   0           1 ];
D1111=C*D111; 

lam=atan(yy1/xx1);
theta=asin(zz1/sqrt(xx1^2+yy1^2+zz1^2));
r=sqrt(xx1^2+yy1^2+zz1^2);
E1=[diff(lam,xx1) diff(lam,yy1) diff(lam,zz1)
    diff(theta,xx1) diff(theta,yy1) diff(theta,zz1)
    diff(r,xx1) diff(r,yy1) diff(r,zz1)];
E11=subs(E1,{xx1,yy1,zz1},{r1*cos(lam1)*cos(theta1),r1*sin(lam1)*cos(theta1),r1*sin(theta1)}); 
E111=simplify(E11);

C=[r1*cos(theta1) 0  0
   0              r1 0
   0              0  1];
 
E1111=C*E111;

% [-sin(lam1)                cos(lam1)                     0 
%  -cos(lam1)*sin(theta1)   -sin(lam1)*sin(theta1)    cos(theta1)
%   cos(lam1)*cos(theta1)    sin(lam1)*cos(theta1)    sin(theta1)] 
               



A(1,1)=diff(xx1,x1);
A(1,2)=diff(xx1,y1);
A(1,3)=diff(xx1,z1);
A(2,1)=diff(yy1,x1);
A(2,2)=diff(yy1,y1);
A(2,3)=diff(yy1,z1);
A(3,1)=diff(zz1,x1);
A(3,2)=diff(zz1,y1);
A(3,3)=diff(zz1,z1);



AA=simplify(subs(A,{x1,y1,z1},{r*cos(lon)*cos(lat),r*sin(lon)*cos(lat),r*sin(lat)}));
B=[r*cos(lat) 0 0
  0         r 0
  0         0 r];
C=simplify(B*AA);
bb=77;
