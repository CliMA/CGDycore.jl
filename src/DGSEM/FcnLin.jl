function FcnGPULin!(F,U,DG,Model,Metric,Exchange,Grid,Cache,CacheS,Phys,Global,::Grids.Quad)
  backend = get_backend(F)
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
  NV = Model.NumV
  NAUX = Model.NumAux
  @views UI = U[:,:,1:DG.NumI,:]
  @views Aux = Cache[:,:,:,1:NAUX]
  FS = CacheS
  @. F = 0
  @. FS = 0

# Auftrieb

  if Model.Buoyancy
    NQ = N * N
    NQG = min(div(NumberThreadGPU,Nz*M),NQ)
    group = (Nz,M,NQG,1)
    ndrange = (Nz,M,NQ,NF)
    if Grid.Form == "Sphere"
      KBuoyancyKernel! = BuoyancySphereKernel!(backend,group)
      KBuoyancyKernel!(FS,U,Metric.X,Metric.J,DG.Glob,Phys;ndrange=ndrange)
    else  
      KBuoyancyKernel! = BuoyancyCartKernel!(backend,group)
      KBuoyancyKernel!(FS,U,DG.Glob,Phys;ndrange=ndrange)
    end  
  end  
  @show "Buo"
  @show FS[end,1,DG.Glob[:,1],3]


  NQ = N * N
  NQG = min(div(NumberThreadGPU,Nz*M),NQ)
  group = (Nz,M,NQG,1)
  ndrange = (Nz,M,NQ,NF)
  KVSp2VCart3Kernel! = VSp2VCart3Kernel!(backend,group)
  @views KVSp2VCart3Kernel!(UI[:,:,:,2:4],Metric.Rotate,DG.Glob;ndrange=ndrange)

  @views Parallels.ExchangeData3DSendGPU(reshape(U,
    Nz*M,size(U,3),NV),Exchange)

  NzG = min(div(NumberThreadGPU,N*N),M*Nz)
  group = (N,N,NzG,1)
  ndrange = (N,N,M*Nz,NF)
  KFluxVolumeNonLinH3Kernel! = FluxVolumeNonLinH3Kernel!(backend,group)
  KFluxVolumeNonLinH3Kernel!(Model.FluxLin,F,U,Aux,Metric.dXdxI,DG.DW,DG.Glob,
    Val(NV),Val(NAUX);ndrange=ndrange)
  @show "VolH",sum(abs.(F))
  @show F[1,1,DG.Glob[:,1],3]

  NDG = min(div(NumberThreadGPU,M*Nz),N*N)
  group = (M,Nz,NDG,1)
  ndrange = (M,Nz,N*N,NF)
  KFluxVolumeNonLinV3Kernel! = FluxVolumeNonLinV3Kernel!(backend,group)
# KFluxVolumeNonLinV3Kernel!(Model.FluxLin,F,U,Aux,Metric.dXdxI,DG.DWT,DG.Glob,
#   Val(NV),Val(NAUX);ndrange=ndrange)
  @show "VolV",sum(abs.(F))
  @show F[1,1,DG.Glob[:,1],3]

  @views Parallels.ExchangeData3DRecvSetGPU!(reshape(U,
    Nz*M,size(U,3),NV),Exchange)

  NEG = min(div(NumberThreadGPU,N*M),Nz)
  group = (N,M,NEG,1)
  ndrange = (N,M,Nz,NE)
  KRiemanNonLinH3Kernel! = RiemanNonLinH3Kernel!(backend,group)
  KRiemanNonLinH3Kernel!(Model.RiemannSolverLin,F,U,Aux,DG.GlobE,Grid.EF,Grid.FE,Metric.NH,
    Metric.VolSurfH,DG.wF,Grid.NumFaces,Val(NV),Val(NAUX);ndrange=ndrange) 
  @show "RieH",sum(abs.(F))
  @show F[1,1,DG.Glob[:,1],3]

  NQ = N * N
  NQG = min(div(NumberThreadGPU,Nz+1),NQ)
  group = (Nz+1,NQG,1)
  ndrange = (Nz+1,NQ,NF)
  KRiemanNonLinV3Kernel! = RiemanNonLinV3Kernel!(backend,group)
  KRiemanNonLinV3Kernel!(Model.RiemannSolverLin,NonConservativeFlux,F,U,Aux,DG.Glob,Metric.NV,
    Metric.VolSurfV,DG.wZ,Val(M),Val(NV),Val(NAUX);ndrange=ndrange) 
  @show "RieV",sum(abs.(F))
  @show F[1,1,DG.Glob[:,1],3]

  NQ = N * N
  NQG = min(div(NumberThreadGPU,Nz*M),NQ)
  group = (Nz,M,NQG,1)
  ndrange = (Nz,M,NQ,NF)
  KScaleMassMatrixKernel! = ScaleMassMatrixKernel!(backend,group)
  KScaleMassMatrixKernel!(F,Metric.J,DG.Glob,Val(NV);ndrange=ndrange)

  NQ = N * N
  NQG = min(div(NumberThreadGPU,Nz*M),NQ)
  group = (Nz,M,NQG,1)
  ndrange = (Nz,M,NQ,NF)
  KVCart2VSp3Kernel! = VCart2VSp3Kernel!(backend,group)
  @views KVCart2VSp3Kernel!(F[:,:,:,2:4],Metric.Rotate,DG.Glob;ndrange=ndrange)
  @show F[1,1,DG.Glob[:,1],3]

  @. F += FS

end

function FcnGPU!(F,U,DG,Model,Metric,Exchange,Grid,CacheU,CacheS,Phys,Global,::Grids.Tri)
  backend = get_backend(F)
  Damp = Model.Damp
  GeoPotential = Model.GeoPotential
  NonConservativeFlux = Model.NonConservativeFlux
  RiemannSolver = Model.RiemannSolver
  FT = eltype(F)
  DoF = DG.DoF
  DoFE = DG.DoFE
  M = DG.OrdPolyZ + 1
  NF = Grid.NumFaces
  NE = Grid.NumEdges
  Nz = Grid.nz
  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU
  Proc = Global.ParallelCom.Proc
  NV = Model.NumV
  NAUX = Model.NumAux
  @views UI = U[:,:,1:DG.NumI,:]
  @views Aux = CacheU[:,:,:,NV+1:NV+NAUX]
  @views p = Aux[:,:,1:DG.NumI,1]
  @views FS = CacheS
  @. F = 0
  @. FS = 0

  @views @. p = Model.Pressure(U[:,:,1:DG.NumI,5])
  if NAUX > 1
    @views GeoPot = Aux[:,:,1:DG.NumI,2]
    DoF = DoF
    DoFG = min(div(NumberThreadGPU,Nz*M),DoF)
    group = (Nz,M,DoFG,1)
    ndrange = (Nz,M,DoF,NF)
    KGeoPotentialKernel! = GeoPotentialKernel!(backend,group)
    KGeoPotentialKernel!(GeoPotential,GeoPot,Metric.X,DG.Glob;ndrange=ndrange)
  end

  if Model.Damping
    DoFG = min(div(NumberThreadGPU,Nz*M),DoF)
    group = (Nz,M,DoFG,1)
    ndrange = (Nz,M,DoF,NF)
    if Grid.Form == "Sphere"
      KDampKernel! = DampSphereKernel!(backend, group)
      KDampKernel!(Damp,FS,U,Metric.X,DG.Glob,Phys;ndrange=ndrange)
    else
      KDampKernel! = DampCartKernel!(backend, group)
      KDampKernel!(Damp,FS,U,Metric.X,DG.Glob;ndrange=ndrange)
    end
  end

  if Model.Coriolis
    DoFG = min(div(NumberThreadGPU,Nz*M),DoF)
    group = (Nz,M,DoFG,1)
    ndrange = (Nz,M,DoF,NF)
    KCoriolisKernel! = CoriolisKernel!(backend,group)
    KCoriolisKernel!(FS,U,Metric.X,DG.Glob,Phys;ndrange=ndrange)
  end

  if Model.Forcing
    NQ = N * N
    NQG = min(div(NumberThreadGPU,Nz*M),NQ)
    group = (Nz,M,NQG,1)
    ndrange = (Nz,M,NQ,NF)
    KForceKernel! = ForceKernel!(backend, group)
    KForceKernel!(Model.Force,FS,U,p,Metric.xS,ndrange=ndrange)
  end


  DoFG = min(div(NumberThreadGPU,Nz*M),DoF)
  group = (Nz,M,DoFG,1)
  ndrange = (Nz,M,DoF,NF)
  KVSp2VCart3Kernel! = VSp2VCart3Kernel!(backend,group)
  @views  KVSp2VCart3Kernel!(UI[:,:,:,2:4],Metric.Rotate,DG.Glob;ndrange=ndrange)
    

  @views Parallels.ExchangeData3DSendGPU(reshape(CacheU[:,:,:,1:NV+NAUX],
    Nz*M,size(CacheU,3),NV+NAUX),Exchange)

  NzG = min(div(NumberThreadGPU,DoF),M*Nz)
  group = (DoF,NzG,1)
  ndrange = (DoF,M*Nz,NF)
  KFluxVolumeNonLinHKernel! = FluxVolumeNonLinHTriKernel!(backend,group)
  KFluxVolumeNonLinHKernel!(Model.FluxAverage,F,U,Aux,Metric.dXdxI,DG.DSx1,DG.DSx2,DG.Glob,
    Val(NV),Val(NAUX);ndrange=ndrange)

  NDG = min(div(NumberThreadGPU,M*Nz),DoF)
  group = (M,Nz,NDG,1)
  ndrange = (M,Nz,DoF,NF)
  KFluxVolumeNonLinV3Kernel! = FluxVolumeNonLinV3Kernel!(backend,group)
  KFluxVolumeNonLinV3Kernel!(Model.FluxAverage,F,U,Aux,Metric.dXdxI,DG.DVZT,DG.Glob,
    Val(NV),Val(NAUX);ndrange=ndrange)

  @views Parallels.ExchangeData3DRecvSetGPU!(reshape(CacheU[:,:,:,1:NV+NAUX],
    Nz*M,size(CacheU,3),NV+NAUX),Exchange)

  NEG = min(div(NumberThreadGPU,DoFE*M),Nz)
  group = (DoFE,M,NzG,1)
  ndrange = (DoFE,M,Nz,NE)
  KRiemanNonLinH3Kernel! = RiemanNonLinH3Kernel!(backend,group)
  KRiemanNonLinH3Kernel!(Model.RiemannSolver,F,U,Aux,DG.GlobE,Grid.EF,Grid.FE,Metric.NH,
    Metric.VolSurfH,DG.wF,Grid.NumFaces,Val(NV),Val(NAUX);ndrange=ndrange) 
  DoFG = min(div(NumberThreadGPU,Nz+1),DoF)
  group = (Nz+1,DoFG,1)
  ndrange = (Nz+1,DoF,NF)
  KRiemanNonLinV3Kernel! = RiemanNonLinV3Kernel!(backend,group)
  KRiemanNonLinV3Kernel!(RiemannSolver,NonConservativeFlux,F,U,Aux,DG.Glob,Metric.NV,
    Metric.VolSurfV,DG.wZ,Val(M),Val(NV),Val(NAUX);ndrange=ndrange) 

  DoFG = min(div(NumberThreadGPU,Nz*M),DoF)
  group = (Nz,M,DoFG,1)
  ndrange = (Nz,M,DoF,NF)
  KScaleMassMatrixKernel! = ScaleMassMatrixKernel!(backend,group)
  KScaleMassMatrixKernel!(F,Metric.J,DG.Glob,Val(NV);ndrange=ndrange)

  DoFG = min(div(NumberThreadGPU,Nz*M),DoF)
  group = (Nz,M,DoFG,1)
  ndrange = (Nz,M,DoF,NF)
  KVCart2VSp3Kernel! = VCart2VSp3Kernel!(backend,group)
  @views KVCart2VSp3Kernel!(F[:,:,:,2:4],Metric.Rotate,DG.Glob;ndrange=ndrange)

  @. F += FS

end

