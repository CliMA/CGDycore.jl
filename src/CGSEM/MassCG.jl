function MassCG(CG,J,Glob,Exchange)
  OrdPoly = CG.OrdPoly
  DoF = CG.DoF
  w = CG.w
  nz = size(J,3)
  NF = size(Glob,2)
  M = zeros(nz,CG.NumG)
  MMass = zeros(nz,CG.NumG)
  MW = zeros(nz-1,CG.NumG)
  for iF = 1 : NF
    iD = 0
    for j = 1 : OrdPoly + 1
      for i = 1 : OrdPoly + 1
        iD += 1
        ind = Glob[iD,iF]  
        for iz = 1 : nz  
          M[iz,ind] += (J[iD,1,iz,iF] + J[iD,2,iz,iF])
          MMass[iz,ind] += 0.5 * (J[iD,1,iz,iF] + J[iD,2,iz,iF]) * w[i] * w[j]
        end  
        for iz = 1 : nz - 1  
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

  M .= FT(0)
  MMass .= FT(0)

  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU

  NzG = min(div(NumberThreadGPU,N*N),Nz)
  group = (N, N, NzG, 1)
  ndrange = (N, N, Nz, NF)

  KMassCGKernel! = MassCGKernel!(backend,group)
  KMassCGKernel!(M,MMass,J,w,Glob,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)

  @views Parallels.ExchangeData!(M[:,:,1],Exchange)
  @views Parallels.ExchangeData!(M[:,:,2],Exchange)
  Parallels.ExchangeData!(MMass,Exchange)
end

@kernel inbounds = true function MassCGKernel!(M,MMass,@Const(JJ),@Const(w),@Const(Glob))
  I,J,Iz,IF = @index(Global, NTuple)

  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  if Iz <= Nz && IF <= NF
    ID = I + (J - 1) * N  
    ind = Glob[ID,IF]  
    @atomic :monotonic M[Iz,ind,1] += JJ[ID,1,Iz,IF]
    @atomic :monotonic M[Iz,ind,2] += JJ[ID,2,Iz,IF]
    @atomic :monotonic MMass[Iz,ind] += eltype(M)(0.5) * (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]) * w[I] * w[J]
  end  
end


