function TransSphere(ksi,eta,z,F,Topo,Param)
X=zeros(3);
X[1]=0.25*((1-ksi)*(1-eta)*F.P[1].x+
  (1+ksi)*(1-eta)*F.P[2].x+
  (1+ksi)*(1+eta)*F.P[3].x+
  (1-ksi)*(1+eta)*F.P[4].x);
X[2]=0.25*((1-ksi)*(1-eta)*F.P[1].y+
  (1+ksi)*(1-eta)*F.P[2].y+
  (1+ksi)*(1+eta)*F.P[3].y+
  (1-ksi)*(1+eta)*F.P[4].y);
X[3]=0.25*((1-ksi)*(1-eta)*F.P[1].z+
  (1+ksi)*(1-eta)*F.P[2].z+
  (1+ksi)*(1+eta)*F.P[3].z+
  (1-ksi)*(1+eta)*F.P[4].z);
r=norm(X);
X=X/r;
X=X*(Param.RadPrint+z);
return X
end

