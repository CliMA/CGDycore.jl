function ProjectW(Fun,CG,Param)
OrdPoly=CG.OrdPoly;
nz=Param.Grid.nz;
p=zeros(CG.NumG,nz,1);
fLoc=zeros(OrdPoly+1,OrdPoly+1);
X = Param.cache.X
J = Param.cache.J
for iz=1:nz-1
  for iF=1:Param.Grid.NumFaces
    for j=1:OrdPoly+1
      for i=1:OrdPoly+1
        x=0.5*(X[i,j,2,:,iF,iz]+X[i,j,1,:,iF,iz+1]);
        det=0.5*(J[i,j,2,iF,iz]+J[i,j,1,iF,iz+1]);
        fLoc[i,j]=Fun(x,Param)*det;
      end
    end
    p[CG.Glob[:,iF],iz]=p[CG.Glob[:,iF],iz]+reshape(fLoc,(OrdPoly+1)*(OrdPoly+1),1);
  end
end
@views p[:,1:nz-1] .= p[:,1:nz-1]./CG.MW;
return p
end
