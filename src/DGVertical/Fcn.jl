function FcnGPUVert!(F,U,DG1,X,dXdxI,J,CacheU,Pressure,Phys,FluxAverage,RiemannSolver)
  backend = get_backend(F)

  NumberThreadGPU = 256
  NV = 3
  NAUX = 2
  Nz = size(U,2)
  M = DG1.OrdPolyZ + 1
  @views Aux = CacheU[:,:,:]
  @views @. Aux[:,:,1] = Pressure(U[:,:,3])
  @views @. Aux[:,:,2] = Phys.Grav * X

  @. F = 0

  NzG = min(div(NumberThreadGPU,M),Nz)
  group = (M,NzG)
  ndrange = (M,Nz)
  KFluxSplitVolumeNonLinVKernel! = FluxSplitVolumeNonLinVertKernel!(backend,group)
  KFluxSplitVolumeNonLinVKernel!(FluxAverage,F,U,Aux,dXdxI,DG1.DVZT,
    Val(NV),Val(NAUX);ndrange=ndrange)

  group = (Nz+1,)
  ndrange = (Nz+1,)
  KRiemanNonLinVKernel! = RiemanNonLinVertKernel!(backend,group)
  KRiemanNonLinVKernel!(RiemannSolver,F,U,Aux,
   DG1.wZ,Val(M),Val(NV),Val(NAUX);ndrange=ndrange) 

  for iv = 1  : NV
    @views @. F[:,:,iv] /= J
  end

end

function FcnSplitAccousticGPUVert!(F,U,DG1,X,dXdxI,J,CacheU,Pressure,Phys,FluxAverage,RiemannSolver)
  backend = get_backend(F)

  NumberThreadGPU = 256
  NV = 2
  NAUX = 0
  Nz = size(U,2)
  M = DG1.OrdPolyZ + 1
  @views Aux = CacheU[:,:,:]

  @. F = 0

  NzG = min(div(NumberThreadGPU,M),Nz)
  group = (M,NzG)
  ndrange = (M,Nz)
  KFluxSplitVolumeNonLinVKernel! = FluxSplitVolumeNonLinVertKernel!(backend,group)
  KFluxSplitVolumeNonLinVKernel!(FluxAverage,F,U,Aux,dXdxI,DG1.DVZT,
    Val(NV),Val(NAUX);ndrange=ndrange)

  group = (Nz+1,)
  ndrange = (Nz+1,)
  KRiemanNonLinVKernel! = RiemanNonLinVertKernel!(backend,group)
  KRiemanNonLinVKernel!(RiemannSolver,F,U,Aux,
   DG1.wZ,Val(M),Val(NV),Val(NAUX);ndrange=ndrange)

  for iv = 1  : NV
    @views @. F[:,:,iv] /= J
  end

end

function FcnAccousticGPUVert!(F,U,DG1,X,dXdxI,J,CacheU,Pressure,Phys,Flux,RiemannSolver)
  backend = get_backend(F)

  NumberThreadGPU = 256
  NV = 2
  NAUX = 0
  Nz = size(U,2)
  M = DG1.OrdPolyZ + 1
  @views Aux = CacheU[:,:,:]

  @. F = 0

  NzG = min(div(NumberThreadGPU,M),Nz)
  group = (M,NzG)
  ndrange = (M,Nz)
  KFluxVolumeNonLinVKernel! = FluxVolumeNonLinVertKernel!(backend,group)
  KFluxVolumeNonLinVKernel!(Flux,F,U,Aux,DG1.DWZ,
    Val(NV),Val(NAUX);ndrange=ndrange)

  group = (Nz+1,)
  ndrange = (Nz+1,)
  KRiemanNonLinVKernel! = RiemanNonLinVertKernel!(backend,group)
  KRiemanNonLinVKernel!(RiemannSolver,F,U,Aux,
   DG1.wZ,Val(M),Val(NV),Val(NAUX);ndrange=ndrange)

  for iv = 1  : NV
    @views @. F[:,:,iv] /= J
  end

end



