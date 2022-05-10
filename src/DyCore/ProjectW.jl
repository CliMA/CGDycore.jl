function ProjectW(Fun,time,CG,Global)
OrdPoly=CG.OrdPoly;
nz=Global.Grid.nz;
p=zeros(nz,CG.NumG);
fLoc=zeros(OrdPoly+1,OrdPoly+1);
X = Global.Metric.X
JF = Global.Metric.JF
for iz=1:nz-1
  for iF=1:Global.Grid.NumFaces
    for j=1:OrdPoly+1
      for i=1:OrdPoly+1
        x=0.5*(X[i,j,2,:,iz,iF]+X[i,j,1,:,iz+1,iF]);
        p[iz,CG.Glob[i,j,iF]]+=Fun(x,time,Global)*JF[i,j,iz+1,iF]
      end
    end
  end
end
@views p[1:nz-1,:] .= p[1:nz-1,:]./CG.MW;
return p
end
