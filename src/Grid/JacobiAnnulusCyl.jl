function [X,J,dXdx]=JacobiAnnulusCyl(ksi,eta,F,Grid)
n1=size(ksi,1);
n2=size(eta,2);
X=zeros(n1,n2,2);
Phi1=atan2(Grid.Nodes(F.N(1)).P(2),Grid.Nodes(F.N(1)).P(1));
Rad1=sqrt(Grid.Nodes(F.N(1)).P(1)^2+Grid.Nodes(F.N(1)).P(2)^2);
Phi2=atan2(Grid.Nodes(F.N(2)).P(2),Grid.Nodes(F.N(2)).P(1));
Rad2=sqrt(Grid.Nodes(F.N(2)).P(1)^2+Grid.Nodes(F.N(2)).P(2)^2);
Phi3=atan2(Grid.Nodes(F.N(3)).P(2),Grid.Nodes(F.N(3)).P(1));
Rad3=sqrt(Grid.Nodes(F.N(3)).P(1)^2+Grid.Nodes(F.N(3)).P(2)^2);
Phi4=atan2(Grid.Nodes(F.N(4)).P(2),Grid.Nodes(F.N(4)).P(1));
Rad4=sqrt(Grid.Nodes(F.N(4)).P(1)^2+Grid.Nodes(F.N(4)).P(2)^2);
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
X(:,:,1)=Phi ;
X(:,:,2)=Rad;

if nargout > 2 
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
  dXdx(:,:,1,1)=dPhidksi;
  dXdx(:,:,2,1)=dRaddksi;
  dXdx(:,:,1,2)=dPhideta;
  dXdx(:,:,2,2)=dRaddeta;
  J=dXdx(:,:,1,1).*dXdx(:,:,2,2)-dXdx(:,:,1,2).*dXdx(:,:,2,1);
end
end

