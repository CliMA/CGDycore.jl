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
  group = (N,N,NF)
  ndrange = (N,N,NF)
  KFluxSplitVolumeNonLinKernel! = FluxSplitVolumeNonLinKernel!(backend,group)
  KFluxSplitVolumeNonLinKernel!(Model.FluxAverage,FV,V,Aux,Metric.dXdxI,Metric.J,DG.DVT,DG.Glob,
    Val(NV),Val(NAUX);ndrange=ndrange)

  NEG = min(div(NumberThreadGPU,N),NE)
  group = (N,NEG)
  ndrange = (N,NE)
  KRiemanNonLinKernel! = RiemanNonLinKernel!(backend,group)
  KRiemanNonLinKernel!(Model.RiemannSolver,FV,V,Aux,DG.Glob,DG.IndE,Grid.EF,Grid.FE,Metric.NH,Metric.T1H,
    Metric.T2H,Metric.VolSurfH,Metric.J,DG.w[1],DG.VZ,Val(NV),Val(NAUX);ndrange=ndrange) 

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

# Source!(F,U,Metric,DG,Grid,Phys)
end
