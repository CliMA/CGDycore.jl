function Buoyancy!(Buo,F,U,Glob,X,NumberThreadGPU)
  backend = get_backend(F)
  M = size(U,1)
  Nz = size(U,2)
  DoF = size(Glob,1)
  NF = size(Glob,2)
  DoFG = min(div(NumberThreadGPU,Nz*M),DoF)
  group = (M,Nz,DoFG,1)
  ndrange = (M,Nz,DoF,NF)
  KBuoyancyKernel! = BuoyancyKernel!(backend,group)
  KBuoyancyKernel!(Buo,F,U,X,Glob,;ndrange=ndrange)
end

@kernel inbounds = true function BuoyancyKernel!(Buo,F,@Const(U),@Const(X),@Const(Glob))
  _,_,iD,  = @index(Local, NTuple)
  K,Iz,ID,IF = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  if ID <= ND
    ind = Glob[ID,IF]
    @views Buo(F[K,Iz,ind,:],U[K,Iz,ind,:],X[ID,K,:,Iz,IF])
  end
end

function GeoPotential!(GPF,Aux,Glob,X,NumberThreadGPU)
  backend = get_backend(Aux)
  M = size(Aux,1)
  Nz = size(Aux,2)
  DoF = size(Glob,1)
  NF = size(Glob,2)
  DoFG = min(div(NumberThreadGPU,Nz*M),DoF)
  group = (Nz,M,DoFG,1)
  ndrange = (Nz,M,DoF,NF)
  KGeoPotentialKernel! = GeoPotentialKernel!(backend,group)
  KGeoPotentialKernel!(GPF,Aux,X,Glob;ndrange=ndrange)
 end

 @kernel inbounds = true function GeoPotentialKernel!(GPF,Aux,@Const(X),@Const(Glob))

  _,_,iD,  = @index(Local, NTuple)
  Iz,K,ID,IF = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  if ID <= ND
    ind = Glob[ID,IF]
    @views GPF(Aux[K,Iz,ind,:],X[ID,K,:,Iz,IF])
  end
end

