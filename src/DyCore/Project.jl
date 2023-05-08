function Project(Fun,time,CG,Global,Param)
OrdPoly=CG.OrdPoly;
nz=Global.Grid.nz;
p=zeros(nz,CG.NumG);
X = Global.Metric.X
JC = Global.Metric.JC
for iF=1:Global.Grid.NumFaces
  for j=1:OrdPoly+1
    for i=1:OrdPoly+1
      ind = CG.Glob[i,j,iF]
      for iz=1:nz
        x=0.5*(X[i,j,1,:,iz,iF]+X[i,j,2,:,iz,iF]);
        p[iz,ind] = Fun(x,time,Global,Param) 
      end
    end
  end
end
return p
end

function ProjectSurf(Fun,time,CG,Global,Param)
OrdPoly=CG.OrdPoly;
NF = Global.Grid.NumFaces
p=zeros(OrdPoly+1,OrdPoly+1,NF);
X = Global.Metric.X
for iF=1:Global.Grid.NumFaces
  for j=1:OrdPoly+1
    for i=1:OrdPoly+1
      x=X[i,j,1,:,1,iF]
      p[i,j,iF] = Fun(x,time,Global,Param)
    end
  end
end
return p
end
