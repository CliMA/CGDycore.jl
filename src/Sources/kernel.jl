function Forcing!(Force,F,U,Aux,FE,Metric,NumberThreadGPU)
  backend = get_backend(F)
  M = size(U,1)
  Nz = size(U,2)
  NumG = size(U,3)
  NumV = size(U,4)
  NumAux = size(Aux,4)
  NumGG = min(div(NumberThreadGPU,Nz*M),NumG)
  group = (M,Nz,NumGG)
  ndrange = (M,Nz,NumG)
  KForceKernel! = ForceKernel!(backend, group)
  KForceKernel!(Force,F,U,Aux,Metric.xS,Val(NumV),Val(NumAux);ndrange=ndrange)
end

@kernel inbounds = true function ForceKernel!(Force,F,@Const(U),@Const(Aux),@Const(xS),
  ::Val{NUMV}, ::Val{NUMAUX}) where {NUMV,NUMAUX}

  _,_,iD  = @index(Local, NTuple)
  K,Iz,ID = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  FLoc = @private eltype(F) (NUMV,)
  ULoc = @private eltype(F) (NUMV,)
  AuxLoc = @private eltype(F) (NUMAUX,)
  xSLoc = @private eltype(F) (2,)

  if ID <= ND
    xSLoc[1] = xS[1,ID]
    xSLoc[2] = xS[2,ID]
    @unroll for iv = 1 : NUMV
      ULoc[iv] = U[K,Iz,ID,iv]
      FLoc[iv] = eltype(F)(0)
    end  
    @unroll for iAux = 1 : NUMAUX
      AuxLoc[iAux] = Aux[K,Iz,ID,iAux]
    end  

    @views Force(FLoc,ULoc,AuxLoc,xSLoc)
    @unroll for iv = 1 : NUMV
      F[K,Iz,ID,iv] += FLoc[iv]
    end  
  end
end 

function Damping!(Damp,F,U,FE::FiniteElements.CGElement,Metric,NumberThreadGPU)
  backend = get_backend(F)
  M = size(U,1)
  Nz = size(U,2)
  NumG = size(U,3)
  NumGG = min(div(NumberThreadGPU,Nz*M),NumG)
  group = (M,Nz,NumGG)
  ndrange = (M,Nz,NumG)
  KDampKernel! = DampCGKernel!(backend, group)
  KDampKernel!(Damp,F,U,Metric.zP;ndrange=ndrange)
end

@kernel inbounds = true function DampCGKernel!(Damp,F,@Const(U),@Const(zP))
  _,_,iD  = @index(Local, NTuple)
  K,Iz,ID = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  if ID <= ND
    @views Fu,Fv,Fw = Damp(zP[Iz,ID],U[K,Iz,ID,:])
    F[K,Iz,ID,2] += Fu
    F[K,Iz,ID,3] += Fv
    F[K,Iz,ID,4] += Fw
  end
end

function Damping!(Damp,F,U,FE::FiniteElements.DGElement,Metric,NumberThreadGPU)
  backend = get_backend(F)
  M = size(U,1)
  Nz = size(U,2)
  DoF = size(FE.Glob,1)
  NF = size(FE.Glob,2)
  DoFG = min(div(NumberThreadGPU,Nz*M),DoF)
  group = (M,Nz,DoFG,1)
  ndrange = (M,Nz,DoF,NF)
  KDampKernel! = DampDGKernel!(backend, group)
  KDampKernel!(Damp,F,U,Metric.X,FE.Glob;ndrange=ndrange)
end

@kernel inbounds = true function DampDGKernel!(Damp,F,@Const(U),@Const(X),@Const(Glob))
  _,_,iD,  = @index(Local, NTuple)
  K,Iz,ID,IF = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  if ID <= ND
    ind = Glob[ID,IF]
    @views Fu,Fv,Fw = Damp(X[ID,K,:,Iz,IF],U[K,Iz,ind,:])
    F[K,Iz,ind,2] += Fu
    F[K,Iz,ind,3] += Fv
    F[K,Iz,ind,4] += Fw
  end
end

function Coriolis!(Cor,F,U,FE::FiniteElements.DGElement,Metric,NumberThreadGPU)
  backend = get_backend(F)
  M = size(U,1)
  Nz = size(U,2)
  DoF = size(FE.Glob,1)
  NF = size(FE.Glob,2)
  DoFG = min(div(NumberThreadGPU,Nz*M),DoF)
  group = (M,Nz,DoFG,1)
  ndrange = (M,Nz,DoF,NF)
  KCoriolisKernel! = CoriolisKernel!(backend,group)
  KCoriolisKernel!(Cor,F,U,Metric.X,FE.Glob,;ndrange=ndrange)
end

@kernel inbounds = true function CoriolisKernel!(Cor,F,@Const(U),@Const(X),@Const(Glob))
  _,_,iD,  = @index(Local, NTuple)
  K,Iz,ID,IF = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  if ID <= ND
    ind = Glob[ID,IF]
    @views Cor(F[K,Iz,ind,:],U[K,Iz,ind,:],X[ID,K,:,Iz,IF])
  end
end

function Buoyancy!(Buo,F,U,Glob,X,NumberThreadGPU)
  backend = get_backend(F)
  M = size(U,1)
  Nz = size(U,2)
  DoF = size(Glob,1)
  NF = size(Glob,2)
  DoFG = min(div(NumberThreadGPU,Nz*M),DoF)
  group = (Nz,M,DoFG,1)
  ndrange = (Nz,M,DoF,NF)
  KBuoyancyKernel! = BuoyancyKernel!(backend,group)
  KBuoyancyKernel!(Buo,F,U,X,Glob,;ndrange=ndrange)
end

@kernel inbounds = true function BuoyancyKernel!(Buo,F,@Const(U),@Const(X),@Const(Glob))
  _,_,iD,  = @index(Local, NTuple)
  Iz,K,ID,IF = @index(Global, NTuple)

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

