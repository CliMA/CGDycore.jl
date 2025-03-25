function Fcn!(F,U,DG,Metric,Grid,Cache,Phys)
  @views V = Cache[:,:,:,1:4]
  @views FV = Cache[:,:,:,5:8]
  @. FV = 0
  VSp2VCart!(V,U,Metric.Rotate)
  FluxVolumeNonLin!(FV,V,DG,Metric.dXdxI,Grid,Phys)
  RiemanNonLin!(FV,V,DG,Metric,Grid,Phys)
  VCart2VSp!(F,FV,Metric.Rotate)
  Source!(F,U,Metric,DG,Grid,Phys)
end  

function FcnSplit!(F,U,DG,Metric,Grid,Cache,Phys)
  @views V = Cache[:,:,:,1:4]
  @views FV = Cache[:,:,:,5:8]
  @. FV = 0
  VSp2VCart!(V,U,Metric.Rotate)
  FluxSplitVolumeNonLin!(FluxNonLinAver!,FV,V,DG,Metric.dXdxI,Grid,Phys)
  RiemanNonLin!(FV,V,DG,Metric,Grid,Phys)
  VCart2VSp!(F,FV,Metric.Rotate)
  Source!(F,U,Metric,DG,Grid,Phys)
end  


function FcnGPUSplit!(F,U,DG,Model,Metric,Grid,Cache,Phys,Global)
  backend = get_backend(F)
  FT = eltype(F)
  N = DG.OrdPoly + 1
  NF = Grid.NumFaces
  NE = Grid.NumEdges
  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU
  @views V = Cache[:,:,:,1:4]
  @views FV = Cache[:,:,:,5:8]
  @views Aux = Cache[:,:,:,9:9]
  @. FV = 0

  hPos = Model.RhoPos
  fac = FT(0.5)
  @views @. Aux[:,:,:,1] = fac * Phys.Grav * U[:,:,:,hPos]^2
  NV = 4
  NAUX = 1

  NFG = min(div(NumberThreadGPU,N*N),NF)
  group = (N*N,NFG)
  ndrange = (N*N,NF)
  KVSp2VCartKernel! = VSp2VCartKernel!(backend,group)
  KVSp2VCartKernel!(V,U,Metric.Rotate;ndrange=ndrange)

  NFG = min(div(NumberThreadGPU,N*N),NF)
  group = (N,N,NFG)
  ndrange = (N,N,NF)
  KFluxSplitVolumeNonLinKernel! = FluxSplitVolumeNonLinKernel!(backend,group)
  KFluxSplitVolumeNonLinKernel!(Model.FluxAverage,FV,V,Aux,Metric.dXdxI,Metric.J,DG.DVT,DG.Glob,
    Val(NV),Val(NAUX);ndrange=ndrange)

  NEG = min(div(NumberThreadGPU,N),NE)
  group = (N,NEG)
  ndrange = (N,NE)
  KRiemanNonLinKernel! = RiemanNonLinKernel!(backend,group)
  KRiemanNonLinKernel!(Model.RiemannSolver,FV,V,Aux,DG.Glob,DG.IndE,Grid.EF,Grid.FE,Metric.NH,Metric.T1H,
    Metric.T2H,Metric.VolSurfH,Metric.J,DG.w[1],Val(NV),Val(NAUX);ndrange=ndrange) 

  NFG = min(div(NumberThreadGPU,N*N),NF)
  group = (N*N,NFG)
  ndrange = (N*N,NF)
  KVCart2VSpKernel! = VCart2VSpKernel!(backend,group)
  KVCart2VSpKernel!(F,FV,Metric.Rotate;ndrange=ndrange)

  NFG = min(div(NumberThreadGPU,N*N),NF)
  group = (N*N,NFG)
  ndrange = (N*N,NF)
  KSourceKernel! = SourceKernel!(backend,group)
  KSourceKernel!(F,U,Metric.X,Phys;ndrange=ndrange)

end

function FcnGPUSplitPar!(F,U,DG,Model,Metric,Exchange,Grid,CacheU,CacheF,Phys,Global)
  backend = get_backend(F)
  FT = eltype(F)
  N = DG.OrdPoly + 1
  NF = Grid.NumFaces
  NE = Grid.NumEdges
  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU
  Proc = Global.ParallelCom.Proc
  @views V = CacheU[:,:,:,1:4]
  @views VI = V[:,:,1:DG.NumI,:]
  @views FV = CacheF[:,:,:,1:4]
  @views Aux = CacheU[:,:,:,5:5]
  @views AuxI = Aux[:,:,1:DG.NumI,:]
  @. FV = 0

  NV = 4
  NAUX = 1

  NFG = min(div(NumberThreadGPU,N*N),NF)
  group = (N*N,NFG)
  ndrange = (N*N,NF)
  KVSp2VCartKernel! = VSp2VCartKernel!(backend,group)
  KVSp2VCartKernel!(VI,U,Metric.Rotate;ndrange=ndrange)
  hPos = Model.RhoPos
  fac = FT(0.5)
  @views @. AuxI[:,:,:,1] = fac * Phys.Grav * VI[:,:,:,hPos]^2

  @views Parallels.ExchangeData3DSendGPU(CacheU[:,1,:,1:NV+NAUX],Exchange)

  NFG = min(div(NumberThreadGPU,N*N),NF)
  group = (N,N,NFG)
  ndrange = (N,N,NF)
  KFluxSplitVolumeNonLinHKernel! = FluxSplitVolumeNonLinHKernel!(backend,group)
  KFluxSplitVolumeNonHLinKernel!(Model.FluxAverage,FV,VI,Aux,Metric.dXdxI,DG.DVT,DG.Glob,
    Val(NV),Val(NAUX);ndrange=ndrange)

  @views Parallels.ExchangeData3DRecvSetGPU!(CacheU[:,1,:,1:NV+NAUX],Exchange)

  NEG = min(div(NumberThreadGPU,N),NE)
  group = (N,NEG)
  ndrange = (N,NE)
  KRiemanNonLinHKernel! = RiemanNonLinHKernel!(backend,group)
  KRiemanNonLinHKernel!(Model.RiemannSolver,FV,V,Aux,DG.GlobE,Grid.EF,Grid.FE,Metric.NH,Metric.T1H,
    Metric.T2H,Metric.VolSurfH,DG.w[1],Grid.NumFaces,Val(NV),Val(NAUX);ndrange=ndrange) 

  NFG = min(div(NumberThreadGPU,N*N),NF)
  group = (N*N,NFG)
  ndrange = (N*N,NF)
  KVCart2VSpKernel! = VCart2VSpKernel!(backend,group)
  KVCart2VSpKernel!(F,FV,Metric.Rotate,Metric.J,;ndrange=ndrange)

  if Model.Coriolis
    NFG = min(div(NumberThreadGPU,N*N),NF)
    group = (N*N,NFG)
    ndrange = (N*N,NF)
    KSourceKernel! = SourceKernel!(backend,group)
    KSourceKernel!(F,U,Metric.X,Phys;ndrange=ndrange)
  end  
end  

function FcnGPUSplitPar3!(F,U,DG,Model,Metric,Exchange,Grid,CacheU,CacheF,Phys,Global)
  backend = get_backend(F)
  FT = eltype(F)
  N = DG.OrdPoly + 1
  M = DG.OrdPolyZ + 1
  NF = Grid.NumFaces
  NE = Grid.NumEdges
  Nz = Grid.nz
  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU
  Proc = Global.ParallelCom.Proc
  NV = 5
  NAUX = 1
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

  @views Parallels.ExchangeData3DSendGPU(reshape(CacheU[:,:,:,1:NV+NAUX],
    Nz*M,size(CacheU,3),NV+NAUX),Exchange)

  NzG = min(div(NumberThreadGPU,N*N),M*Nz)
  group = (N,N,NzG,1)
  ndrange = (N,N,M*Nz,NF)
  KFluxSplitVolumeNonLinH3Kernel! = FluxSplitVolumeNonLinH3Kernel!(backend,group)
  KFluxSplitVolumeNonLinH3Kernel!(Model.FluxAverage,FV,V,Aux,Metric.dXdxI,DG.DVT,DG.Glob,
    Val(NV),Val(NAUX);ndrange=ndrange)


  NDG = min(div(NumberThreadGPU,M*Nz),N*N)
  group = (M,Nz,NDG,1)
  ndrange = (M,Nz,N*N,NF)
  KFluxSplitVolumeNonLinV3Kernel! = FluxSplitVolumeNonLinV3Kernel!(backend,group)
  KFluxSplitVolumeNonLinV3Kernel!(Model.FluxAverage,FV,V,Aux,Metric.dXdxI,DG.DVT,DG.Glob,
    Val(NV),Val(NAUX);ndrange=ndrange)

  @views Parallels.ExchangeData3DRecvSetGPU!(reshape(CacheU[:,:,:,1:NV+NAUX],
    Nz*M,size(CacheU,3),NV+NAUX),Exchange)

  NEG = min(div(NumberThreadGPU,N*M),Nz)
  group = (N,M,NzG,1)
  ndrange = (N,M,Nz,NE)
  KRiemanNonLinH3Kernel! = RiemanNonLinH3Kernel!(backend,group)
  KRiemanNonLinH3Kernel!(Model.RiemannSolver,FV,V,Aux,DG.GlobE,Grid.EF,Grid.FE,Metric.NH,Metric.T1H,
    Metric.T2H,Metric.VolSurfH,DG.w[1],Grid.NumFaces,Val(NV),Val(NAUX);ndrange=ndrange) 

  NQ = N * N
  NQG = min(div(NumberThreadGPU,Nz+1),NQ)
  group = (Nz+1,NQG,1)
  ndrange = (Nz+1,NQ,NF)
  KRiemanNonLinV3Kernel! = RiemanNonLinV3Kernel!(backend,group)
  KRiemanNonLinV3Kernel!(Model.RiemannSolver,FV,V,Aux,DG.Glob,Metric.NV,Metric.T1V,
    Metric.T2V,Metric.VolSurfV,DG.wZ[1],M,Val(NV),Val(NAUX);ndrange=ndrange) 


  NQ = N * N
  NQG = min(div(NumberThreadGPU,Nz*M),NQ)
  group = (Nz,M,NQG,1)
  ndrange = (Nz,M,NQ,NF)
  KVCart2VSp3Kernel! = VCart2VSp3Kernel!(backend,group)
  KVCart2VSp3Kernel!(F,FV,Metric.Rotate,Metric.J,DG.Glob;ndrange=ndrange)

  if Model.Coriolis
    NFG = min(div(NumberThreadGPU,N*N),NF)
    group = (N*N,NFG)
    ndrange = (N*N,NF)
    KSourceKernel! = SourceKernel!(backend,group)
    KSourceKernel!(F,U,Metric.X,Phys;ndrange=ndrange)
  end  
  if Model.Buoyancy
    NDG = min(div(NumberThreadGPU,Nz*M),size(F,3))
    group = (Nz,M,NDG)
    ndrange = (Nz,M,size(F,3))
    KBuoyancyKernel! = BuoyancyKernel!(backend,group)
    KBuoyancyKernel!(F,U,Phys;ndrange=ndrange)
  end    

end
