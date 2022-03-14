function [X]=Cart(DG,P1,P2,TypeE,Param)

if strcmp(TypeE,'Y')
  ksi=DG.xwY;
  n=DG.OrdPolyY+1;
  X=zeros(n,2);
  X(:,1)=.5*(P1(1)*(1-ksi)+...
    P2(1)*(1+ksi));
  X(:,2)=.5*(P1(2)*(1-ksi)+...
    P2(2)*(1+ksi));
else
  ksi=DG.xwX;
  n=DG.OrdPolyX+1;
  X=zeros(n,2);
  X(:,1)=.5*(P1(1)*(1-ksi)+...
    P2(1)*(1+ksi));
  X(:,2)=.5*(P1(2)*(1-ksi)+...
    P2(2)*(1+ksi));
  X(:,1)=.5*(P1(1)*(1-ksi)+...
    P2(1)*(1+ksi));
  hSX1=hS(P1(1),Param);
  hSX2=hS(P2(1),Param);
  if abs(hSX2-hSX1)>0
    for i=1:n
      X(i,2)=P1(2)+(P2(2)-P1(2))/(hSX2-hSX1)*(hS(X(i,1),Param)-hSX1);
    end
  else
    X(:,2)=.5*(P1(2)*(1-ksi)+...
      P2(2)*(1+ksi));
  end
end
end

