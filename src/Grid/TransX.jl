function TransSphereX(ksi,eta,zeta,X,CG,Global)
  XP=zeros(3);
  for j=1:OrdPoly+1
    for i=1:OrdPoly+1
      for k=1:OrdPolyZ+1  
        XP[:]=XP[:]+reshape(Lagrange(ksi,CG.xw,i)*
          Lagrange(eta,CG.xw,j)*Lagrange(eta,CG.xwZ,k)*X[i,j,k,:];
      end
    end
  end
  r=norm(X);
  z=max(r-Global.Grid.RadEarth,0.0)
  X=X/r;
  X=X*(Global.Output.RadPrint+z);
  return XP
end

