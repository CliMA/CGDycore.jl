function MassCG(CG,J,Glob,Global)
  OrdPoly=CG.OrdPoly;
  w=CG.w
  nz=Global.Grid.nz;
  M=zeros(nz,CG.NumG);
  MMass=zeros(nz,CG.NumG);
  MW =zeros(nz-1,CG.NumG);
  @inbounds for iF = 1 : Global.Grid.NumFaces
    @inbounds for iz = 1 : nz  
      @inbounds for j = 1 : OrdPoly + 1
        @inbounds for i = 1 : OrdPoly + 1
          ind = Glob[i,j,iF]  
          M[iz,ind] += (J[i,j,1,iz,iF]+J[i,j,2,iz,iF])
          MMass[iz,ind] += 0.5 * (J[i,j,1,iz,iF]+J[i,j,2,iz,iF])*w[i]*w[j]
        end
      end
    end
    @inbounds for iz = 1 : nz - 1  
      @inbounds for j = 1 : OrdPoly + 1
        @inbounds for i = 1 : OrdPoly+1
          ind = Glob[i,j,iF]  
          MW[iz,ind] += (J[i,j,2,iz,iF]+J[i,j,1,iz+1,iF])
        end
      end
    end
  end
  ExchangeData!(M,Global.Exchange)
  ExchangeData!(MMass,Global.Exchange)
  ExchangeData!(MW,Global.Exchange)
  return (M,MW,MMass)
end
