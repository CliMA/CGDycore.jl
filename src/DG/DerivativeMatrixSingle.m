function [DW,DS,DV,DVT,B]=DerivativeMatrixSingle(OrdPoly)
[w,xw]=GaussLobattoQuad(OrdPoly);
DS=zeros(OrdPoly+1,OrdPoly+1);
for i=1:OrdPoly+1
  for j=1:OrdPoly+1
    DS(i,j)=DLagrange(xw(i),xw,j);
  end
end
DW=-inv(diag(w))*DS'*diag(w);
DV=2*DS;
DV(1,1)=2*DS(1,1)+1/w(1);
DV(OrdPoly+1,OrdPoly+1)=2*DS(OrdPoly+1,OrdPoly+1)-1/w(OrdPoly+1);
DVT=DV';
B=zeros(OrdPoly+1,OrdPoly+1);
B(1,1)=-1/w(1);
B(OrdPoly+1,OrdPoly+1)=1/w(OrdPoly+1);
end
