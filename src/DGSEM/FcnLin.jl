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
  KFluxVolumeLinH3Kernel! = FluxVolumeLinH3Kernel!(backend,group)
  KFluxVolumeLinH3Kernel!(F,U,Aux,Metric.dXdxI,DG.DW,DG.Glob,
    Val(NV),Val(NAUX);ndrange=ndrange)

  NDG = min(div(NumberThreadGPU,M*Nz),N*N)
  group = (M,Nz,NDG,1)
  ndrange = (M,Nz,N*N,NF)
  KFluxVolumeLinV3Kernel! = FluxVolumeLinV3Kernel!(backend,group)
  KFluxVolumeLinV3Kernel!(F,U,Aux,Metric.dXdxI,DG.DWZ,DG.Glob,
    Val(NV),Val(NAUX);ndrange=ndrange)

  @views Parallels.ExchangeData3DRecvSetGPU!(reshape(U,
    Nz*M,size(U,3),NV),Exchange)

  NEG = min(div(NumberThreadGPU,N*M),Nz)
  group = (N,M,NEG,1)
  ndrange = (N,M,Nz,NE)
  KRiemanNonLinH3Kernel! = RiemanNonLinH3Kernel!(backend,group)
  KRiemanNonLinH3Kernel!(Model.RiemannSolverLin,F,U,Aux,DG.GlobE,Grid.EF,Grid.FE,Metric.NH,
    Metric.VolSurfH,DG.wF,Grid.NumFaces,Val(NV),Val(NAUX);ndrange=ndrange) 

  NQ = N * N
  NQG = min(div(NumberThreadGPU,Nz+1),NQ)
  group = (Nz+1,NQG,1)
  ndrange = (Nz+1,NQ,NF)
  KRiemanNonLinV3Kernel! = RiemanNonLinV3Kernel!(backend,group)
  KRiemanNonLinV3Kernel!(Model.RiemannSolverLin,NonConservativeFlux,F,U,Aux,DG.Glob,Metric.NV,
    Metric.VolSurfV,DG.wZ,Val(M),Val(NV),Val(NAUX);ndrange=ndrange) 

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

  @. F += FS

end

function FcnGPULin1!(F,U,DG,Model,Metric,Exchange,Grid,Cache,CacheS,Phys,Global,::Grids.Quad)
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
# KFluxVolumeLinH3Kernel! = FluxVolumeNonLinH3Kernel!(backend,group)
# KFluxVolumeLinH3Kernel!(Model.FluxLin,F,U,Aux,Metric.dXdxI,DG.DW,DG.Glob,
#   Val(NV),Val(NAUX);ndrange=ndrange)
  KFluxVolumeLinH3Kernel! = FluxVolumeLinH3Kernel!(backend,group)
  KFluxVolumeLinH3Kernel!(F,U,Aux,Metric.dXdxI,DG.DW,DG.Glob,
    Val(NV),Val(NAUX);ndrange=ndrange)

  NDG = min(div(NumberThreadGPU,M*Nz),N*N)
  group = (M,Nz,NDG,1)
  ndrange = (M,Nz,N*N,NF)
# KFluxVolumeLinV3Kernel! = FluxVolumeNonLinV3Kernel!(backend,group)
# KFluxVolumeLinV3Kernel!(Model.FluxLin,F,U,Aux,Metric.dXdxI,DG.DWZ,DG.Glob,
#   Val(NV),Val(NAUX);ndrange=ndrange)
  KFluxVolumeLinV3Kernel! = FluxVolumeLinV3Kernel!(backend,group)
  KFluxVolumeLinV3Kernel!(F,U,Aux,Metric.dXdxI,DG.DWZ,DG.Glob,
    Val(NV),Val(NAUX);ndrange=ndrange)

  @views Parallels.ExchangeData3DRecvSetGPU!(reshape(U,
    Nz*M,size(U,3),NV),Exchange)

  NEG = min(div(NumberThreadGPU,N*M),Nz)
  group = (N,M,NEG,1)
  ndrange = (N,M,Nz,NE)
  KRiemanNonLinH3Kernel! = RiemanNonLinH3Kernel!(backend,group)
  KRiemanNonLinH3Kernel!(Model.RiemannSolverLin,F,U,Aux,DG.GlobE,Grid.EF,Grid.FE,Metric.NH,
    Metric.VolSurfH,DG.wF,Grid.NumFaces,Val(NV),Val(NAUX);ndrange=ndrange) 

  NQ = N * N
  NQG = min(div(NumberThreadGPU,Nz+1),NQ)
  group = (Nz+1,NQG,1)
  ndrange = (Nz+1,NQ,NF)
  KRiemanNonLinV3Kernel! = RiemanNonLinV3Kernel!(backend,group)
  KRiemanNonLinV3Kernel!(Model.RiemannSolverLin,NonConservativeFlux,F,U,Aux,DG.Glob,Metric.NV,
    Metric.VolSurfV,DG.wZ,Val(M),Val(NV),Val(NAUX);ndrange=ndrange) 

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

  @. F += FS

end

