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
  Parallels.ExchangeData!(M,Exchange)
  Parallels.ExchangeData!(MMass,Exchange)
  Parallels.ExchangeData!(MW,Exchange)
  return (M,MW,MMass)
end

function MassCGGPU!(CG,J,Glob,Exchange,Global)
  backend = get_backend(J)
  FT = eltype(J)
  N = CG.OrdPoly + 1
  DoF = CG.DoF
  w = CG.w
  Nz = size(J,3)
  NF = size(Glob,2)
  M = CG.M
  MMass = CG.MMass
  MW = CG.MW

  @. M = FT(0)
  @. MMass = FT(0)
  @. MW = FT(0)

  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU

  NzG = min(div(NumberThreadGPU,N*N),Nz)
  group = (N, N, NzG, 1)
  ndrange = (N, N, Nz, NF)

  KMassCGKernel! = MassCGKernel!(backend,group)
  KMassCGKernel!(M,MMass,MW,J,w,Glob,ndrange=ndrange)

  Parallels.ExchangeData!(M,Exchange)
  Parallels.ExchangeData!(MMass,Exchange)
  Parallels.ExchangeData!(MW,Exchange)
end

@kernel function MassCGKernel!(M,MMass,MW,@Const(JJ),@Const(w),@Const(Glob))
  I,J,Iz,IF = @index(Global, NTuple)

  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N  
  @inbounds ind = Glob[ID,IF]  
  if Iz <= Nz && IF <= NF
    @inbounds @atomic M[Iz,ind] += (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF])
    @inbounds @atomic MMass[Iz,ind] += 0.5 * (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]) * w[I] * w[J]
  end  
  if Iz < Nz && IF <= NF
    @inbounds @atomic MW[Iz,ind] += (JJ[ID,2,Iz,IF] + JJ[ID,1,Iz+1,IF])
  end  
end


