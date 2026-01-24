function Coriolis!(F,U,DG,NumberThreadGPU,Phys)
  backend = get_backend(F)
  M = size(U,1)
  Nz = size(U,2)
  DoF = size(DG.Glob,1)
  NF = size(DG.Glob,2)
  DoFG = min(div(NumberThreadGPU,Nz*M),DoF)
  group = (Nz,M,DoFG,1)
  ndrange = (Nz,M,DoF,NF)
  KCoriolisKernel! = CoriolisKernel!(backend,group)
  KCoriolisKernel!(F,U,DG.Glob,Phys;ndrange=ndrange)
end    

@kernel inbounds = true function CoriolisKernel!(F,U,Glob,Phys)

  _,_,iD,  = @index(Local, NTuple)
  Iz,K,ID,IF = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]
  uPos = 2
  vPos = 3
  if ID <= ND
    ind = Glob[ID,IF]
    fac = eltype(F)(2.0) * Phys.Omega
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

  uPos = 2
  vPos = 3
  wPos = 4
  if ID <= ND
    x1 = X[ID,K,1,Iz,IF]  
    x2 = X[ID,K,1,Iz,IF]  
    x3 = X[ID,K,1,Iz,IF]  
    ind = Glob[ID,IF]
    Rad = sqrt(x1^2 + x2^2 + x3^2) - Phys.RadEarth
    h = Rad - Phys.RadEarth
    Damp = -DampF(h) / Rad^2 * 
      (U[K,Iz,ind,uPos] * x1 + U[K,Iz,ind,vPos] * x2 + U[K,Iz,ind,wPos] * x3)
    F[K,Iz,ind,uPos] += Damp * x1
    F[K,Iz,ind,vPos] += Damp * x2
    F[K,Iz,ind,wPos] += Damp * x3
  end
end

#=
Radius = SQRT(SUM(Elem_xGP(:,i,j,k,iElem)*Elem_xGP(:,i,j,k,iElem)))
      zPloc=Radius-RadEarth
      Ut(2:4,i,j,k,iElem) = Ut(2:4,i,j,k,iElem) - Grav*(RadEarth/Radius)**2*Elem_xGP(:,i,j,k,iElem)/Radius*U(1,i,j,k,iElem)
      Ut(2:4,i,j,k,iElem) = Ut(2:4,i,j,k,iElem) - 2.0*OmegaEarth*CROSS((/0.0,0.0,1.0/),U(2:4,i,j,k,iElem))
      IF (zPLoc>=Height-StrideDamp) THEN
        Damp = Relax*DampF((1.0 - (Height - zPLoc)/StrideDamp))
        Ut(2:4,i,j,k,iElem) = Ut(2:4,i,j,k,iElem) - &
             Damp*SUM(U(2:4,i,j,k,iElem)*Elem_xGP(:,i,j,k,iElem))/Radius*Elem_xGP(:,i,j,k,iElem)/Radius
=#             

@kernel inbounds = true function DampCartKernel!(Damp,F,@Const(U),@Const(X),@Const(Glob))
  _,_,iD,  = @index(Local, NTuple)
  Iz,K,ID,IF = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  if ID <= ND
    ind = Glob[ID,IF]
    h = X[ID,K,3,Iz,IF]
    Fu,Fv,Fw = Damp(h,view(U,K,Iz,ind,1:5))
    F[K,Iz,ind,1] += Fu
    F[K,Iz,ind,2] += Fv
    F[K,Iz,ind,3] += Fw
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

@kernel inbounds = true function ForceKernel!(Force,F,U,p,xS)
  _,_,iD,  = @index(Local, NTuple)
  Iz,K,ID = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  if ID <= ND
    Force(view(F,K,Iz,ID,:),view(U,K,Iz,ID,:),p[K,Iz,ID],xS[2,ID])
  end
end
