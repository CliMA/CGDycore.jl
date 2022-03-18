function [w,x] = GaussLegendreQuad(n)
w=zeros(n+1,1);
x=zeros(n+1,1);
switch n
  case 0
    w(1)=2;
    x(1)=0;
  case 1
    w(1)=1;
    w(2)=1;
    x(1)=-1/sqrt(3);
    x(2)=1/sqrt(3);
  case 2
    w(2)=8/9;
    w(1)=5/9;
    w(3)=5/9;
    
    x(2)=0;
    x(1)=-sqrt(3/5);
    x(3)=sqrt(3/5);
    
  case 3
    w(2)=(18+sqrt(30))/36;
    w(3)=(18+sqrt(30))/36;
    w(1)=(18-sqrt(30))/36;
    w(4)=(18-sqrt(30))/36;
    
    x(2)=-sqrt(3/7-2/7*sqrt(6/5));
    x(3)=sqrt(3/7-2/7*sqrt(6/5));
    x(1)=-sqrt(3/7+2/7*sqrt(6/5));
    x(4)=sqrt(3/7+2/7*sqrt(6/5));
    
  case 4
    w(3)=128/225;
    w(2)=(322+13*sqrt(70))/900;
    w(4)=(322+13*sqrt(70))/900;
    w(1)=(322-13*sqrt(70))/900;
    w(5)=(322-13*sqrt(70))/900;
    
    x(3)=0;
    x(2)=-1/3*sqrt(5-2*sqrt(10/7));
    x(4)=1/3*sqrt(5-2*sqrt(10/7));
    x(1)=-1/3*sqrt(5+2*sqrt(10/7));
    x(5)=1/3*sqrt(5+2*sqrt(10/7));
    
  otherwise
    error('ord1+ord2 zu groﬂ')
    
end
end
