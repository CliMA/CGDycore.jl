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

  NzG = min(div(NumberThreadGPU,M),Nz)
  group = (M,NzG)
  ndrange = (M,Nz)
  KFluxSplitVolumeNonLinVKernel! = FluxSplitVolumeNonLinVertKernel!(backend,group)
  KFluxSplitVolumeNonLinVKernel!(FluxAverage,F,U,Aux,dXdxI,DG1.DVZT,
    Val(NV),Val(NAUX);ndrange=ndrange)

  group = (Nz+1,)
  ndrange = (Nz+1,)
  @show group,ndrange
  KRiemanNonLinVKernel! = RiemanNonLinVertKernel!(backend,group)
  KRiemanNonLinVKernel!(RiemannSolver,F,U,Aux,
   DG1.wZ,Val(M),Val(NV),Val(NAUX);ndrange=ndrange) 

#=
  Damp = Model.Damp
  GeoPotential = Model.GeoPotential
  NonConservativeFlux = Model.NonConservativeFlux
  RiemannSolver = Model.RiemannSolver
  FT = eltype(F)
  N = DG.OrdPoly + 1
  M = DG.OrdPolyZ + 1
  NF = Grid.NumFaces
  NE = Grid.NumEdges
  Nz = Grid.nz
  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU
  Proc = Global.ParallelCom.Proc
  NV = 5
  NAUX = 2
  @views V = CacheU[:,:,:,1:NV]
  @views VI = V[:,:,1:DG.NumI,:]
  @views FV = CacheF[:,:,:,1:NV]
  @views Aux = CacheU[:,:,:,NV+1:NV+NAUX]
  @views AuxI = Aux[:,:,1:DG.NumI,:]
  @. FV = 0


  NQ = N * N
  NQG = min(div(NumberThreadGPU,Nz*M),NQ)
  group = (Nz,M,NQG,1)
  ndrange = (Nz,M,NQ,NF)
  KVSp2VCart3Kernel! = VSp2VCart3Kernel!(backend,group)
  KVSp2VCart3Kernel!(VI,U,Metric.Rotate,DG.Glob;ndrange=ndrange)
  @views @. AuxI[:,:,:,1] = Model.Pressure(U[:,:,1:DG.NumI,5])
  NQ = N * N
  NQG = min(div(NumberThreadGPU,Nz*M),NQ)
  group = (Nz,M,NQG,1)
  ndrange = (Nz,M,NQ,NF)
  KGeoPotentialKernel! = GeoPotentialKernel!(backend,group)
  @views KGeoPotentialKernel!(GeoPotential,AuxI[:,:,:,2],Metric.X,DG.Glob;ndrange=ndrange)
    

  @views Parallels.ExchangeData3DSendGPU(reshape(CacheU[:,:,:,1:NV+NAUX],
    Nz*M,size(CacheU,3),NV+NAUX),Exchange)

  NzG = min(div(NumberThreadGPU,N*N),M*Nz)
  group = (N,N,NzG,1)
  ndrange = (N,N,M*Nz,NF)
  KFluxSplitVolumeNonLinHKernel! = FluxSplitVolumeNonLinHQuadKernel!(backend,group)
  KFluxSplitVolumeNonLinHKernel!(Model.FluxAverage,FV,V,Aux,Metric.dXdxI,DG.DVT,DG.Glob,
    Val(NV),Val(NAUX);ndrange=ndrange)

  NDG = min(div(NumberThreadGPU,M*Nz),N*N)
  group = (M,Nz,NDG,1)
  ndrange = (M,Nz,N*N,NF)
  KFluxSplitVolumeNonLinV3Kernel! = FluxSplitVolumeNonLinV3Kernel!(backend,group)
  KFluxSplitVolumeNonLinV3Kernel!(Model.FluxAverage,FV,V,Aux,Metric.dXdxI,DG.DVZT,DG.Glob,
    Val(NV),Val(NAUX);ndrange=ndrange)

  @views Parallels.ExchangeData3DRecvSetGPU!(reshape(CacheU[:,:,:,1:NV+NAUX],
    Nz*M,size(CacheU,3),NV+NAUX),Exchange)

  NEG = min(div(NumberThreadGPU,N*M),Nz)
  group = (N,M,NzG,1)
  ndrange = (N,M,Nz,NE)
  KRiemanNonLinH3Kernel! = RiemanNonLinH3Kernel!(backend,group)
  KRiemanNonLinH3Kernel!(Model.RiemannSolver,FV,V,Aux,DG.GlobE,Grid.EF,Grid.FE,Metric.NH,
    Metric.VolSurfH,DG.wF,Grid.NumFaces,Val(NV),Val(NAUX);ndrange=ndrange) 

  NQ = N * N
  NQG = min(div(NumberThreadGPU,Nz+1),NQ)
  group = (Nz+1,NQG,1)
  ndrange = (Nz+1,NQ,NF)
  KRiemanNonLinV3Kernel! = RiemanNonLinV3Kernel!(backend,group)
  KRiemanNonLinV3Kernel!(RiemannSolver,NonConservativeFlux,FV,V,Aux,DG.Glob,Metric.NV,
    Metric.VolSurfV,DG.wZ,Val(M),Val(NV),Val(NAUX);ndrange=ndrange) 

  NQ = N * N
  NQG = min(div(NumberThreadGPU,Nz*M),NQ)
  group = (Nz,M,NQG,1)
  ndrange = (Nz,M,NQ,NF)
  KScaleMassMatrixKernel! = ScaleMassMatrixKernel!(backend,group)
  KScaleMassMatrixKernel!(FV,Metric.J,DG.Glob,Val(NV);ndrange=ndrange)

  if Model.Buoyancy
    NQ = N * N
    NQG = min(div(NumberThreadGPU,Nz*M),NQ)
    group = (Nz,M,NQG,1)
    ndrange = (Nz,M,NQ,NF)
    if Grid.Form == "Sphere"
      KBuoyancyKernel! = BuoyancySphereKernel!(backend,group)
      KBuoyancyKernel!(FV,V,Metric.X,DG.Glob,Phys;ndrange=ndrange)
    else  
      KBuoyancyKernel! = BuoyancyCartKernel!(backend,group)
      KBuoyancyKernel!(FV,V,DG.Glob,Phys;ndrange=ndrange)
    end  
  end  

  NQ = N * N
  NQG = min(div(NumberThreadGPU,Nz*M),NQ)
  group = (Nz,M,NQG,1)
  ndrange = (Nz,M,NQ,NF)
  KVCart2VSp3Kernel! = VCart2VSp3Kernel!(backend,group)
  KVCart2VSp3Kernel!(F,FV,Metric.Rotate,DG.Glob;ndrange=ndrange)

  if Model.Damping
    NQ = N * N
    NQG = min(div(NumberThreadGPU,Nz*M),NQ)
    group = (Nz,M,NQG,1)
    ndrange = (Nz,M,NQ,NF)
    if Grid.Form == "Sphere"
      KDampKernel! = DampSphereKernel!(backend, group)
      KDampKernel!(Damp,F,U,Metric.X,DG.Glob,Phys;ndrange=ndrange)
    else  
      KDampKernel! = DampCartKernel!(backend, group)
      KDampKernel!(Damp,F,U,Metric.X,DG.Glob;ndrange=ndrange)
    end  
  end  

  if Model.Coriolis
    NQ = N * N
    NQG = min(div(NumberThreadGPU,Nz*M),NQ)
    group = (Nz,M,NQG,1)
    ndrange = (Nz,M,NQ,NF)
    KCoriolisKernel! = CoriolisKernel!(backend,group)
    KCoriolisKernel!(F,U,Metric.X,DG.Glob,Phys;ndrange=ndrange)
  end  
  if Model.Buoyancy
    NDG = min(div(NumberThreadGPU,Nz*M),size(F,3))
    group = (Nz,M,NDG)
    ndrange = (Nz,M,size(F,3))
    KBuoyancyKernel! = BuoyancyKernel!(backend,group)
    KBuoyancyKernel!(F,U,Phys;ndrange=ndrange)
  end    
=#  

end

