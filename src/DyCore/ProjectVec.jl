function ProjectVec(Fun,CG,Global)
OrdPoly=CG.OrdPoly;
OrdPolyZ=CG.OrdPolyZ;
nz=Global.Grid.nz;
uS=zeros(nz,CG.NumG);
vS=zeros(nz,CG.NumG);
uLoc=zeros(OrdPoly+1,OrdPoly+1);
vLoc=zeros(OrdPoly+1,OrdPoly+1);
JC = Global.Metric.JC
X = Global.Metric.X
for iz=1:nz
  for iF=1:Global.Grid.NumFaces
    for j=1:OrdPoly+1
      for i=1:OrdPoly+1
        x=0.5*(X[i,j,1,:,iz,iF]+X[i,j,2,:,iz,iF]);
        (uLoc,vLoc)=Fun(x,Global);
        uS[iz,CG.Glob[i,j,iF]]+=uLoc*JC[i,j,iz,iF]
        vS[iz,CG.Glob[i,j,iF]]+=vLoc*JC[i,j,iz,iF]
      end
    end
  end
end
uS=uS./CG.M;
vS=vS./CG.M;
return (uS,vS)
end
