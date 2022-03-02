function [p]=Project(Fun,CG,Param)
OrdPoly=CG.OrdPoly;
nz=Param.Grid.nz;
p=zeros(CG.NumG,nz,1);
fLoc=zeros(OrdPoly+1,OrdPoly+1);
for iz=1:nz
  for iF=1:Param.Grid.NumFaces
    for j=1:OrdPoly+1
      for i=1:OrdPoly+1
        x=0.5*(Param.X(i,j,1,:,iF,iz)+Param.X(i,j,2,:,iF,iz));
        det=Param.JC(i,j,iF,iz);
        fLoc(i,j)=Fun(x,Param)*det;
      end
    end
    p(CG.Faces(iF).Glob,iz)=p(CG.Faces(iF).Glob,iz)+reshape(fLoc,(OrdPoly+1)*(OrdPoly+1),1);
  end
end
p=p./CG.M;
end
