clear all
syms lon lat r
syms x1 y1 z1

x=cos(lon)*cos(lat);
y=sin(lon)*cos(lat);
z=sin(lat);


D(1,1)=diff(x,lon);
D(2,1)=diff(y,lon);
D(3,1)=diff(z,lon);
D(1,2)=diff(x,lat);
D(2,2)=diff(y,lat);
D(3,2)=diff(z,lat);
DPen=simplify(inv(transpose(D)*D)*transpose(D));

DD=simplify(inv(D));

D1(1,1)=diff(lon1,x1);
D1(1,2)=diff(lon1,y1);
D1(1,3)=diff(lon1,z1);
D1(2,1)=diff(lat1,x1);
D1(2,2)=diff(lat1,y1);
D1(2,3)=diff(lat1,z1);

DPen=inv(transpose(D)*D)*transpose(D);
DPen=simplify(DPen);
A(1,1)=cos(lat);
A(1,2)=0;
A(2,1)=0;
A(2,2)=1;
DD=A*DPen;
aa=3;