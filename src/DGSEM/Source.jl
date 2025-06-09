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
    F[K,Iz,ind,uPos] += fac * U[K,Iz,ind,vPos]  
    F[K,Iz,ind,vPos] += -fac * U[K,Iz,ind,uPos]
  end  
end

@kernel inbounds = true function BuoyancyKernel!(Fun,F,@Const(U),@Const(X),@Const(Glob))

  _,_,iD,  = @index(Local, NTuple)
  Iz,K,ID,IF = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  if ID <= ND
    ind = Glob[ID,IF]
    @views Fun(F[K,Iz,ind,:],U[K,Iz,ind,:],X[ID,K,:,Iz,IF])
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

    fac = -Phys.Grav * (Phys.RadEarth / r)^2 / r * U[K,Iz,ind,RhoPos]
    F[K,Iz,ind,uPos] += fac * X[ID,K,1,Iz,IF]
    F[K,Iz,ind,vPos] += fac * X[ID,K,2,Iz,IF]
    F[K,Iz,ind,wPos] += fac * X[ID,K,3,Iz,IF]
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

    fac = -Phys.Grav * U[K,Iz,ind,RhoPos]
    F[K,Iz,ind,wPos] += fac
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
    Fu,Fv,Fw = Damp(h,view(U,K,Iz,ind,1:5))
    F[K,Iz,ind,2] += Fu
    F[K,Iz,ind,3] += Fv
    F[K,Iz,ind,4] += Fw
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
    Fu,Fv,Fw = Damp(h,view(U,K,Iz,ind,1:5))
    F[K,Iz,ind,2] += Fu
    F[K,Iz,ind,3] += Fv
    F[K,Iz,ind,4] += Fw
  end
end

@kernel inbounds = true function GeoPotentialKernel!(GPF,GP,@Const(X),@Const(Glob))

  _,_,iD,  = @index(Local, NTuple)
  Iz,K,ID,IF = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  if ID <= ND
    ind = Glob[ID,IF]
    GP[K,Iz,ind] = GPF(X[ID,K,1,Iz,IF],X[ID,K,2,Iz,IF],X[ID,K,3,Iz,IF])
  end
end
