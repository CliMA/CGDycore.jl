function DerivativeMatrixSingle(OrdPoly)
(xw,wmat)=gausslobatto(OrdPoly+1);
DS=zeros(OrdPoly+1,OrdPoly+1);
for i=1:OrdPoly+1
  for j=1:OrdPoly+1
    DS[i,j]=DLagrange(xw[i],xw,j);
  end
end
w = vec(wmat)
DW=-inv(diagm(w))*DS'*diagm(w);
DV=2*DS;
DV[1,1]=2*DS[1,1]+1/w[1];
DV[OrdPoly+1,OrdPoly+1]=2*DS[OrdPoly+1,OrdPoly+1]-1/w[OrdPoly+1];
DVT=DV';
B=zeros(OrdPoly+1,OrdPoly+1);
B[1,1]=-1/w[1];
B[OrdPoly+1,OrdPoly+1]=1/w[OrdPoly+1];
return (DW,DS,DV,DVT,B)
end
