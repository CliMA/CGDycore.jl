function TransCart(ksi,eta,z,F,Topo,Grid)
X=zeros(3);
X[1]=0.25*((1-ksi)*(1-eta)*F.P[1].x+
  (1+ksi)*(1-eta)*F.P[2].x+
  (1+ksi)*(1+eta)*F.P[3].x+
  (1-ksi)*(1+eta)*F.P[4].x);
X[2]=0.25*((1-ksi)*(1-eta)*F.P[1].y+
  (1+ksi)*(1-eta)*F.P[2].y+
  (1+ksi)*(1+eta)*F.P[3].y+
  (1-ksi)*(1+eta)*F.P[4].y);
X[3]=first(Topo(X[1],X[2],z,Grid.Topography));
return X
end

function TransCartX(ksi,eta,zeta,X,CG,Global)
OrdPoly=CG.OrdPoly
  OrdPolyZ=CG.OrdPolyZ
  XP=zeros(3);
  for j=1:OrdPoly+1
    for i=1:OrdPoly+1
      for k=1:OrdPolyZ+1
        XP[:]=XP[:]+Lagrange(ksi,CG.xw,i)*
          Lagrange(eta,CG.xw,j)*Lagrange(zeta,CG.xwZ,k)*X[i,j,k,:];
      end
    end
  end
  if abs(XP[1])<1.e-20
    XP[1] = 0.0  
  end  
  if abs(XP[2])<1.e-20
    XP[2] = 0.0  
  end  
  if abs(XP[3])<1.e-20
    XP[3] = 0.0  
  end  
  return XP
end  



