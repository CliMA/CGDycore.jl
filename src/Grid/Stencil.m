syms h1 h2
syms x a0 a1
syms c1 c2
x1 = 0;
x0 = x1 - h1;
x2 = x1 + h1;
x3 = x2 + h2;

f1 = int(1,x);
A11 = subs(f1,x,x2)-subs(f1,x,x1);
A21 = subs(f1,x,x3)-subs(f1,x,x2);
f2 = int(x,x);
A12 = subs(f2,x,x2)-subs(f2,x,x1);
A22 = subs(f2,x,x3)-subs(f2,x,x2);

A=[A11 A12
   A21 A22];
r=[c1*h1
   c2*h2];
s=A\r; 
g0 = int(s(1)+s(2)*x,x);
g0 = simplify((subs(g0,x,x1)-subs(g0,x,x0))/(x1-x0));

val_g0=subs(g0,{c2,c1,h1,h2},{5,4,1,1});


aa=3;