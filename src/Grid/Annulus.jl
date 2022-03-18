function [X]=Annulus(DG,P1,P2)
ksi=DG.xw;
n=DG.OrdPoly+1;
X=zeros(n,2);
Rad1=sqrt(P1(1)^2+P1(2)^2);
Rad2=sqrt(P2(1)^2+P2(2)^2);
Rad=0.5*(Rad1*(1-ksi)+...
  Rad2*(1+ksi));
X(:,1)=.5*((1-ksi)*P1(1)+(1+ksi)*P2(1));
X(:,2)=.5*((1-ksi)*P1(2)+(1+ksi)*P2(2));
X=X./sqrt((X(:,1).^2+X(:,2).^2));
X=X.*Rad;
end

