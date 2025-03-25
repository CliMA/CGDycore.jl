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

@kernel inbounds = true function BuoyancyKernel!(F,U,Phys)

  _,_,iD  = @index(Local, NTuple)
  Iz,K,ID = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  if ID <= ND
    RhoPos = 1
    wPos = 4
    F[Iz,K,ID,wPos] += -Phys.Grav * U[Iz,K,ID,RhoPos]
  end
end
