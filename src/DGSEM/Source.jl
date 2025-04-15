@kernel inbounds = true function CoriolisKernel!(F,U,X,Glob,Phys)

  _,_,iD,  = @index(Local, NTuple)
  Iz,K,ID,IF = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  if ID <= ND
    uPos = 2  
    vPos = 3  
    ind = Glob[ID,IF]
    fac = eltype(F)(2.0) * Phys.Omega * X[ID,K,3,Iz,IF] / sqrt(X[ID,K,1,Iz,IF]^2 + 
      X[ID,K,2,Iz,IF]^2 + X[ID,K,3,Iz,IF]^2)
    F[Iz,K,ind,uPos] += fac * U[Iz,K,ind,vPos]  
    F[Iz,K,ind,vPos] += -fac * U[Iz,K,ind,uPos]
  end  
end

@kernel inbounds = true function BuoyancyKernel!(Fun,F,@Const(U),@Const(X),@Const(Glob))

  _,_,iD,  = @index(Local, NTuple)
  Iz,K,ID,IF = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  if ID <= ND
    ind = Glob[ID,IF]
    @views Fun(F[Iz,K,ind,:],U[Iz,K,ind,:],X[ID,K,:,Iz,IF])
  end
end

@kernel inbounds = true function BuoyancySphereKernel!(F,@Const(U),@Const(X),@Const(Glob),Phys)

  _,_,iD,  = @index(Local, NTuple)
  Iz,K,ID,IF = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  if ID <= ND
    RhoPos = 1
    uPos = 2
    vPos = 3
    wPos = 4
    ind = Glob[ID,IF]
    r = sqrt(X[ID,K,1,Iz,IF]^2 +
      X[ID,K,2,Iz,IF]^2 + X[ID,K,3,Iz,IF]^2) 

    fac = -Phys.Grav * (Phys.RadEarth / r)^2 / r * U[Iz,K,ind,RhoPos]
    F[Iz,K,ind,uPos] += fac * X[ID,K,1,Iz,IF]
    F[Iz,K,ind,vPos] += fac * X[ID,K,2,Iz,IF]
    F[Iz,K,ind,wPos] += fac * X[ID,K,3,Iz,IF]
  end
end

@kernel inbounds = true function BuoyancyCartKernel!(F,@Const(U),@Const(Glob),Phys)

  _,_,iD,  = @index(Local, NTuple)
  Iz,K,ID,IF = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  if ID <= ND
    RhoPos = 1
    wPos = 4
    ind = Glob[ID,IF]

    fac = -Phys.Grav * U[Iz,K,ind,RhoPos]
    F[Iz,K,ind,wPos] += fac
  end
end

@kernel inbounds = true function DampSphereKernel!(Damp,F,@Const(U),@Const(X),@Const(Glob),Phys)
  _,_,iD,  = @index(Local, NTuple)
  Iz,K,ID,IF = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  if ID <= ND
    RhoPos = 1
    uPos = 2
    vPos = 3
    wPos = 4
    ind = Glob[ID,IF]
    h = sqrt(X[ID,K,1,Iz,IF]^2 +
      X[ID,K,2,Iz,IF]^2 + X[ID,K,3,Iz,IF]^2) - Phys.RadEarth
    Fu,Fv,Fw = Damp(h,view(U,Iz,K,ind,1:5))
    F[Iz,K,ind,2] += Fu
    F[Iz,K,ind,3] += Fv
    F[Iz,K,ind,4] += Fw
  end
end

@kernel inbounds = true function DampCartKernel!(Damp,F,@Const(U),@Const(X),@Const(Glob))
  _,_,iD,  = @index(Local, NTuple)
  Iz,K,ID,IF = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  if ID <= ND
    RhoPos = 1
    uPos = 2
    vPos = 3
    wPos = 4
    ind = Glob[ID,IF]
    h = X[ID,K,3,Iz,IF]
    Fu,Fv,Fw = Damp(h,view(U,Iz,K,ind,1:5))
    F[Iz,K,ind,2] += Fu
    F[Iz,K,ind,3] += Fv
    F[Iz,K,ind,4] += Fw
  end
end

@kernel inbounds = true function GeoPotentialKernel!(GPF,GP,@Const(X),@Const(Glob))

  _,_,iD,  = @index(Local, NTuple)
  Iz,K,ID,IF = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  if ID <= ND
    ind = Glob[ID,IF]
    GP[Iz,K,ind] = GPF(X[ID,K,1,Iz,IF],X[ID,K,2,Iz,IF],X[ID,K,3,Iz,IF])
  end
end
