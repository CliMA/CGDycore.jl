OrdPoly=3;
[wX,xw]=GaussLobattoQuad(OrdPoly);
Dx=zeros(OrdPoly+1,OrdPoly+1);
  for i=1:OrdPoly+1
    for j=1:OrdPoly+1
      Dx(i,j)=DLagrange(xw(i),xw,j);
    end
  end
a=rand(OrdPoly+1,1);
b=rand(OrdPoly+1,1);

r=sum(wX.*(a.*(Dx*b)+b.*(Dx*a)));
s=sum(wX.*(Dx*(a.*b)));

ee=4;