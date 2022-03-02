function [wVec,xBarVec] = GetQuadMeth(n,QuadType,Type)
switch Type
  case 'Quad'
    if strcmp(QuadType,'GaussLobattoQuad')==1
      [w,x]=GaussLobattoQuad(n);
    elseif strcmp(QuadType,'GaussLegendreQuad')==1
      [w,x]=GaussLegendreQuad(n);
    end
    x=x;
    w=w;
    wVec=zeros(size(w,2)^2,1);
    xBarVec=zeros(size(w,2)^2,2);
    iVec=1;
    for i=1:size(w,2)
      for j=1:size(w,2)
        wVec(iVec)=w(i)*w(j);
        xBarVec(iVec,1)=x(i);
        xBarVec(iVec,2)=x(j);
        iVec=iVec+1;
      end
    end
  case 'Edge'
    if strcmp(QuadType,'GaussLobattoQuad')==1
      [w,x]=GaussLobattoQuad(n);
    elseif strcmp(QuadType,'GaussLegendreQuad')==1
      [w,x]=GaussLegendreQuad(n);
    end
    xBarVec=x';
    wVec=w';
  case 'Triangle'
    switch n
      case 1
        xBarVec=zeros(1,2);
        wVec=zeros(1,1);
        xBarVec(1,1)=1/3;
        xBarVec(1,2)=1/3;
        wVec(1)=0.5;
      case 2
        xBarVec=zeros(3,2);
        wVec=zeros(3,1);
        xBarVec(1,1)=1/6;
        xBarVec(1,2)=1/6;
        wVec(1)=1/6;
        xBarVec(2,1)=1/6;
        xBarVec(2,2)=2/3;
        wVec(2)=1/6;
        xBarVec(3,1)=2/3;
        xBarVec(3,2)=1/6;
        wVec(3)=1/6;
      case 3
        xBarVec=zeros(4,2);
        wVec=zeros(4,1);
        xBarVec(1,1)=1/3;
        xBarVec(1,2)=1/3;
        wVec(1)=-27/96;
        xBarVec(2,1)=1/5;
        xBarVec(2,2)=1/5;
        wVec(2)=25/96;
        xBarVec(3,1)=1/5;
        xBarVec(3,2)=3/5;
        wVec(3)=25/96;
        xBarVec(4,1)=3/5;
        xBarVec(4,2)=1/5;
        wVec(4)=25/96;
    end
end
end

function [w,x] = GaussLegendreQuad(n)
w=zeros(n,1);
x=zeros(n,1);
switch n
  case 0
  case 1
    w(1)=2;
    x(1)=0;
  case 2
    w(1)=1;
    w(2)=1;
    x(1)=-1/sqrt(3);
    x(2)=1/sqrt(3);
    
  case 3
    w(1)=8/9;
    w(2)=5/9;
    w(3)=5/9;
    
    x(1)=0;
    x(2)=-sqrt(3/5);
    x(3)=sqrt(3/5);
    
  case 4
    w(1)=(18+sqrt(30))/36;
    w(2)=(18+sqrt(30))/36;
    w(3)=(18-sqrt(30))/36;
    w(4)=(18-sqrt(30))/36;
    
    x(1)=-sqrt(3/7-2/7*sqrt(6/5));
    x(2)=sqrt(3/7-2/7*sqrt(6/5));
    x(3)=-sqrt(3/7+2/7*sqrt(6/5));
    x(4)=sqrt(3/7+2/7*sqrt(6/5));
    
  case 5
    w(1)=128/225;
    w(2)=(322+13*sqrt(70))/900;
    w(3)=(322+13*sqrt(70))/900;
    w(4)=(322-13*sqrt(70))/900;
    w(5)=(322-13*sqrt(70))/900;
    
    x(1)=0;
    x(2)=-1/3*sqrt(5-2*sqrt(10/7));
    x(3)=1/3*sqrt(5-2*sqrt(10/7));
    x(4)=-1/3*sqrt(5+2*sqrt(10/7));
    x(5)=1/3*sqrt(5+2*sqrt(10/7));
    
  otherwise
    error('ord1+ord2 zu groﬂ')
    
end
w=w';
x=x';
end

function [w,x] = GaussLobattoQuad(n)
w=zeros(n,1);
x=zeros(n,1);
switch n
  case 1
    w(1)=1;
    x(1)=0;
  case 2
    w(1)=1;
    w(2)=1;
    x(1)=-1;
    x(2)=1;
    
  case 3
    w(1)=1/3;
    w(2)=4/3;
    w(3)=1/3;
    
    x(1)=-1;
    x(2)=0;
    x(3)=1;
    
  case 4
    w(1)=1/6;
    w(2)=5/6;
    w(3)=5/6;
    w(4)=1/6;
    
    x(1)=-1;
    x(2)=-1/sqrt(5);
    x(3)=1/sqrt(5);
    x(4)=1;
  case 5
    w(1)=1/10;
    w(2)=49/90;
    w(3)=32/45;
    w(4)=49/90;
    w(5)=1/10;
    
    x(1)=-1;
    x(2)=-sqrt(3/7);
    x(3)=0;
    x(4)=sqrt(3/7);
    x(5)=1;
    
  case 6
    w(1)=(14+sqrt(7))/30;
    w(2)=(14+sqrt(7))/30;
    w(3)=(14-sqrt(7))/30;
    w(4)=(14-sqrt(7))/30;
    w(5)=1/15;
    w(6)=1/15;
    
    x(1)=-sqrt(1/3-2*sqrt(7)/21);
    x(2)=sqrt(1/3-2*sqrt(7)/21);
    x(3)=-sqrt(1/3+2*sqrt(7)/21);
    x(4)=sqrt(1/3+2*sqrt(7)/21);
    x(5)=-1;
    x(6)=1;
    
  case 7
    w(1)=256/525;
    w(2)=(124+7*sqrt(15))/350;
    w(3)=(124+7*sqrt(15))/350;
    w(4)=(124-7*sqrt(15))/350;
    w(5)=(124-7*sqrt(15))/350;
    w(6)=1/21;
    w(7)=1/21;
    
    x(1)=0;
    x(2)=-sqrt(5/11-2/11*sqrt(5/3));
    x(3)=sqrt(5/11-2/11*sqrt(5/3));
    x(4)=-sqrt(5/11+2/11*sqrt(5/3));
    x(5)=sqrt(5/11+2/11*sqrt(5/3));
    x(6)=-1;
    x(7)=1;
    
  otherwise
    error('ord1+ord2 zu groﬂ oder zu klein')
    
end
w=w';
x=x';

end
