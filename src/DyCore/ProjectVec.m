function [uS,vS]=ProjectVec(Fun,CG,Param)
OrdPoly=CG.OrdPoly;
OrdPolyZ=CG.OrdPolyZ;
nz=Param.Grid.nz;
uS=zeros(CG.NumG,nz,1);
vS=zeros(CG.NumG,nz,1);
uLoc=zeros(OrdPoly+1,OrdPoly+1);
vLoc=zeros(OrdPoly+1,OrdPoly+1);
for iz=1:nz
  for iF=1:Param.Grid.NumFaces
    for j=1:OrdPoly+1
      for i=1:OrdPoly+1
        x=0.5*(Param.X(i,j,1,:,iF,iz)+Param.X(i,j,2,:,iF,iz));
        det=Param.JC(i,j,iF,iz);
        [uLoc(i,j),vLoc(i,j)]=Fun(x,Param);
        uLoc(i,j)=uLoc(i,j)*det;
        vLoc(i,j)=vLoc(i,j)*det;
      end
    end
    uS(CG.Faces(iF).Glob,iz)=uS(CG.Faces(iF).Glob,iz)+reshape(uLoc,(OrdPoly+1)*(OrdPoly+1),1);
    vS(CG.Faces(iF).Glob,iz)=vS(CG.Faces(iF).Glob,iz)+reshape(vLoc,(OrdPoly+1)*(OrdPoly+1),1);
  end
end
uS=uS./CG.M;
vS=vS./CG.M;
end
