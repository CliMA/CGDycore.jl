function f=GradGrad(OrdPoly)
DG.OrdPoly=OrdPoly;
[DG.w,DG.xw]=GaussLobattoQuad(OrdPoly);
[DG.DW,DG.DS,DG.DV,DG.DVT]=DerivativeMatrix(DG);
c=zeros(OrdPoly+1,1);
for i=1:OrdPoly+1
  c(i)=(i-OrdPoly/2-1)^2;
end
f=zeros(OrdPoly+1,1);
for i=1:OrdPoly+1
  for l=1:OrdPoly+1
    grad_phi_i=DG.DS(l,i);
    w=DG.w(l);
    t=0;
    for j=1:OrdPoly+1
      t=t+c(j)*DG.DS(l,j);
    end
    f(i)=f(i)+w*t*grad_phi_i;
  end
  f(i)=f(i)/DG.w(i);
end
f2=DG.DS*c;
f2=DG.DW*f2;
end