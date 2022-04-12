function MassCG(CG,Global)
OrdPoly=CG.OrdPoly;
nz=Global.Grid.nz;
M=zeros(nz,CG.NumG);
J = Global.Metric.J
for iF=1:Global.Grid.NumFaces
  for iz=1:nz  
    iG=0
    for j=1:OrdPoly+1
      for i=1:OrdPoly+1
        iG += 1  
        ind = CG.Glob[iG,iF]  
          M[iz,ind] += 0.5*(J[i,j,1,iz,iF]+J[i,j,2,iz,iF])
      end
    end
  end
end
MW=0.5*(M[1:end-1,:]+M[2:end,:]);
return (M,MW)
end
