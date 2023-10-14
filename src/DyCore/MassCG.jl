function MassCG(CG,J,Glob,Exchange)
  OrdPoly = CG.OrdPoly
  DoF = CG.DoF
  w = CG.w
  nz = size(J,3)
  NF = size(Glob,2)
  M = zeros(nz,CG.NumG)
  MMass = zeros(nz,CG.NumG)
  MW = zeros(nz-1,CG.NumG)
  @inbounds for iF = 1 : NF
    iD = 0
    @inbounds for j = 1 : OrdPoly + 1
      @inbounds for i = 1 : OrdPoly + 1
        iD += 1
        ind = Glob[iD,iF]  
        @inbounds for iz = 1 : nz  
          M[iz,ind] += (J[iD,1,iz,iF] + J[iD,2,iz,iF])
          MMass[iz,ind] += 0.5 * (J[iD,1,iz,iF] + J[iD,2,iz,iF]) * w[i] * w[j]
        end  
        @inbounds for iz = 1 : nz - 1  
          MW[iz,ind] += (J[iD,2,iz,iF] + J[iD,1,iz+1,iF])
        end
      end
    end
  end
  ExchangeData!(M,Exchange)
  ExchangeData!(MMass,Exchange)
  ExchangeData!(MW,Exchange)
  return (M,MW,MMass)
end
