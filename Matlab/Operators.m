syms J11 J12 
syms J21 J22
syms J31 J32 J33
syms JI11 JI12 
syms JI21 JI22
syms JI31 JI32 JI33
syms u v w
syms uu vv ww
syms dx dy dz
syms dxc dyc dzc

J=[J11 J12 0
   J21 J22 0
   J31 J32 J33];
% J=[J11 J12 0
%    J21 J22 0
%    0  0 1]; 
invSJ=inv(J)*det(J);
invJ=simplify(invSJ);
invJ=[JI11 JI12 0
       JI21 JI22 0
       JI31 JI32 JI33];
U=[u;v;0];
W=[0;0;w];
% Curl operator
cU1=cross(invJ(1,:),U);
cU2=cross(invJ(2,:),U);
cU3=cross(invJ(3,:),U);
grad=transpose(invJ(1,:)*dxc+invJ(2,:)*dyc+invJ(3,:)*dzc);

cU=transpose(dx*cU1+dy*cU2+dz*cU3);
cUV=transpose(dx*cU1+dy*cU2);
cUH=transpose(dz*cU3);
cW1=cross(invJ(1,:),W);
cW2=cross(invJ(2,:),W);
cW3=cross(invJ(3,:),W);
cW=transpose(dx*cW1+dy*cW2+dz*cW3);
UU=[uu;vv;0];
WW=[0;0;ww];
c=cU+cW;
cH=subs(c,{w},{0});

Adv=cross(UU+WW,cU+cW);

AdvH=subs(Adv,{w,ww},{0,0});
Adv1uv=subs(Adv(1),{w},{0});
Adv2uv=subs(Adv(2),{w},{0});
Adv3uv=subs(Adv(3),{w},{0});
Adv1w=subs(Adv(1),{u,v},{0,0});
Adv2w=subs(Adv(2),{u,v},{0,0});
Adv3w=subs(Adv(3),{u,v},{0,0});
%Div operator
UContr=invJ*(UU+WW);
%UContrI=invJE*(UU+WW);
div=dx*UContr(1)+dy*UContr(2)+dz*UContr(3);


aa=3.