function [X,J,dXdx]=JacobiAnnulus(ksi,eta,F,Grid)
n1=size(ksi,1);
n2=size(eta,2);
X=zeros(n1,n2,2);
Phi1=atan2(F.P(2,1),F.P(1,1));
Rad1=sqrt(F.P(1,1)^2+F.P(2,1)^2);
Phi2=atan2(F.P(2,2),F.P(1,2));
Rad2=sqrt(F.P(1,2)^2+F.P(2,2)^2);
Phi3=atan2(F.P(2,3),F.P(1,3));
Rad3=sqrt(F.P(1,3)^2+F.P(2,3)^2);
Phi4=atan2(F.P(2,4),F.P(1,4));
Rad4=sqrt(F.P(1,4)^2+F.P(2,4)^2);
if Phi1<0
  Phi1=Phi1+2*pi;
end
if Phi2<0
  Phi2=Phi2+2*pi;
end  
if Phi3<0
  Phi3=Phi3+2*pi;
end
if Phi4<0
  Phi4=Phi4+2*pi;
end
if Phi1>Phi3
  Phi3=Phi3+2*pi;
end
if Phi2>Phi4
  Phi4=Phi4+2*pi;
end
Phi=0.25*(Phi1*(1-ksi)*(1-eta)+...
  Phi2*(1+ksi)*(1-eta)+...
  Phi3*(1+ksi)*(1+eta)+...
  Phi4*(1-ksi)*(1+eta));
Rad=0.25*(Rad1*(1-ksi)*(1-eta)+...
  Rad2*(1+ksi)*(1-eta)+...
  Rad3*(1+ksi)*(1+eta)+...
  Rad4*(1-ksi)*(1+eta));
X(:,:,1)=cos(Phi).*Rad ;
X(:,:,2)=sin(Phi).*Rad;

if nargout > 2
  dX1dPhi=-sin(Phi).*Rad;
  dX1dRad=cos(Phi);
  dX2dPhi=cos(Phi).*Rad;
  dX2dRad=sin(Phi);
  
  dPhidksi=0.25*((-Phi1+Phi2+...
    Phi3-Phi4)*ones(n1,1)*ones(1,n2)+...
    ( Phi1-Phi2+...
    Phi3-Phi4)*ones(n1,1)*eta);
  dPhideta=0.25*((-Phi1-Phi2+...
    Phi3+Phi4)*ones(n1,1)*ones(1,n2)+...
    ( Phi1-Phi2+...
    Phi3-Phi4)*ksi*ones(1,n2));
  dRaddksi=0.25*((-Rad1+Rad2+...
    Rad3-Rad4)*ones(n1,1)*ones(1,n2)+...
    ( Rad1-Rad2+...
    Rad3-Rad4)*ones(n1,1)*eta);
  dRaddeta=0.25*((-Rad1-Rad2+...
    Rad3+Rad4)*ones(n1,1)*ones(1,n2)+...
    ( Rad1-Rad2+...
    Rad3-Rad4)*ksi*ones(1,n2));
  
  dXdx=zeros(n1,n2,2,2);
  dXdx(:,:,1,1)=dX1dPhi.*dPhidksi+dX1dRad.*dRaddksi;
  dXdx(:,:,2,1)=dX2dPhi.*dPhidksi+dX2dRad.*dRaddksi;
  dXdx(:,:,1,2)=dX1dPhi.*dPhideta+dX1dRad.*dRaddeta;
  dXdx(:,:,2,2)=dX2dPhi.*dPhideta+dX2dRad.*dRaddeta;
  J=dXdx(:,:,1,1).*dXdx(:,:,2,2)-dXdx(:,:,1,2).*dXdx(:,:,2,1);
end
end

