function Source!(F,U,Metric,DG,Grid,Phys)
  uPos = 2
  vPos = 3
  Omega = Phys.Omega
  NX = DG.OrdPoly + 1
  iDG = 0
  X = Metric.X
  @inbounds for iF = 1 : Grid.NumFaces
    @inbounds for iQ = 1 : NX * NX  
      iDG += 1  
      fac = 2.0 * Omega * X[iQ,1,3,1,iF] / sqrt(X[iQ,1,1,1,iF]^2 + 
        X[iQ,1,2,1,iF]^2 + X[iQ,1,3,1,iF]^2)
      F[1,1,iDG,uPos] += fac * U[1,1,iDG,vPos]  
      F[1,1,iDG,vPos] += -fac * U[1,1,iDG,uPos]  
    end
  end  
end


@kernel inbounds = true function SourceKernel!(F,U,X,Phys)

  iQ,  = @index(Local, NTuple)
  _,IF = @index(Global, NTuple)

  NQ = @uniform @ndrange()[1]
  NF = @uniform @ndrange()[2]
  FT = eltype(F)

  if IF <= NF
    uPos = 2  
    vPos = 3  
    iD = iQ + (IF - 1) * NQ
    fac = FT(2.0) * Phys.Omega * X[iQ,1,3,1,IF] / sqrt(X[iQ,1,1,1,IF]^2 + 
      X[iQ,1,2,1,IF]^2 + X[iQ,1,3,1,IF]^2)
    F[1,1,iD,uPos] += fac * U[1,1,iD,vPos]  
    F[1,1,iD,vPos] += -fac * U[1,1,iD,uPos]  
  end
end
