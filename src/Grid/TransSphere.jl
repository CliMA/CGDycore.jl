function TransSphere(ksi,eta,z,F,Topo,Param)
X=zeros(3)
X[1]=0.25*((1-ksi)*(1-eta)*F.P[1].x+
  (1+ksi)*(1-eta)*F.P[2].x+
  (1+ksi)*(1+eta)*F.P[3].x+
  (1-ksi)*(1+eta)*F.P[4].x)
X[2]=0.25*((1-ksi)*(1-eta)*F.P[1].y+
  (1+ksi)*(1-eta)*F.P[2].y+
  (1+ksi)*(1+eta)*F.P[3].y+
  (1-ksi)*(1+eta)*F.P[4].y)
X[3]=0.25*((1-ksi)*(1-eta)*F.P[1].z+
  (1+ksi)*(1-eta)*F.P[2].z+
  (1+ksi)*(1+eta)*F.P[3].z+
  (1-ksi)*(1+eta)*F.P[4].z)
r=norm(X)
X=X/r
X=X*(Param.RadPrint+z)
return X
end

function TransSphereS(ksi,eta,F,Rad)
X=zeros(3)
X[1]=0.25*((1-ksi)*(1-eta)*F.P[1].x+
  (1+ksi)*(1-eta)*F.P[2].x+
  (1+ksi)*(1+eta)*F.P[3].x+
  (1-ksi)*(1+eta)*F.P[4].x)
X[2]=0.25*((1-ksi)*(1-eta)*F.P[1].y+
  (1+ksi)*(1-eta)*F.P[2].y+
  (1+ksi)*(1+eta)*F.P[3].y+
  (1-ksi)*(1+eta)*F.P[4].y)
X[3]=0.25*((1-ksi)*(1-eta)*F.P[1].z+
  (1+ksi)*(1-eta)*F.P[2].z+
  (1+ksi)*(1+eta)*F.P[3].z+
  (1-ksi)*(1+eta)*F.P[4].z)
r=norm(X)
X=X/r
X=X*Rad
end

function TransSphereX(ksi,eta,zeta,X,CG,Global)
  OrdPoly=CG.OrdPoly
  OrdPolyZ=CG.OrdPolyZ
  XP=zeros(3)
  for j=1:OrdPoly+1
    for i=1:OrdPoly+1
      for k=1:OrdPolyZ+1
        XP[:]=XP[:]+Lagrange(ksi,CG.xw,i)*
          Lagrange(eta,CG.xw,j)*Lagrange(zeta,CG.xwZ,k)*X[i,j,k,:]
      end
    end
  end
  r=norm(XP)
  z=max(r-Global.Grid.Rad,0.0)
  XP=XP/r
  XP=XP*(Global.Output.RadPrint+z)
  return XP
end

