function [X]=TransSphere(ksi,eta,z,F,Topo,Param)
X=zeros(1,3);
X(1:3)=0.25*((1-ksi)*(1-eta)*F.P(1:3,1)...
  +(1+ksi)*(1-eta)*F.P(1:3,2)...
  +(1+ksi)*(1+eta)*F.P(1:3,3)...
  +(1-ksi)*(1+eta)*F.P(1:3,4));
r=norm(X);
X=X/r;
X=X*(Param.RadPrint+z);
end

