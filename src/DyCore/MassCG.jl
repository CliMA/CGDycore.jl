function MassCG(CG,Global)
OrdPoly=CG.OrdPoly;
w=CG.w
nz=Global.Grid.nz;
M=zeros(nz,CG.NumG);
MMass=zeros(Float64,nz,CG.NumG);
J = Global.Metric.J
for iF=1:Global.Grid.NumFaces
  for iz=1:nz  
    for j=1:OrdPoly+1
      for i=1:OrdPoly+1
        ind = CG.Glob[i,j,iF]  
          M[iz,ind] += 0.5 * (J[i,j,1,iz,iF]+J[i,j,2,iz,iF])
          MMass[iz,ind] += 0.5 * (J[i,j,1,iz,iF]+J[i,j,2,iz,iF])*w[i]*w[j]
      end
    end
  end
end

ExchangeData!(M,Global.Exchange)
ExchangeData!(MMass,Global.Exchange)
MW=(M[1:end-1,:]+M[2:end,:]);
return (M,MW,MMass)
end
