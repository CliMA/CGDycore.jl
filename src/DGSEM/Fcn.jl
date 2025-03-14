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


function FcnGPUSplit!(F,U,DG,Model,Metric,Grid,Cache,Phys)
  backend = get_backend(F)
  FT = eltype(F)
  @views V = Cache[:,:,:,1:4]
  @views FV = Cache[:,:,:,5:8]
  @views Aux = Cache[:,:,:,9]
  @. FV = 0

  hPos = Model.RhoPos
  @views @. Aux[:,:,:,1] = FT(0.5) * Phys.Grav * U[:,:,:,hPos]^2

  NV = 4
  NAUX = 1

  group = ((DG.OrdPoly+1)*(DG.OrdPoly+1),10)
  ndrange = ((DG.OrdPoly+1)*(DG.OrdPoly+1),Grid.NumFaces)
  KVSp2VCartKernel! = VSp2VCartKernel!(backend,group)
  KVSp2VCartKernel!(V,U,Metric.Rotate;ndrange=ndrange)

  group = (DG.OrdPoly+1,DG.OrdPoly+1,10)
  ndrange = (DG.OrdPoly+1,DG.OrdPoly+1,Grid.NumFaces)
  KFluxSplitVolumeNonLinKernel! = FluxSplitVolumeNonLinKernel!(backend,group)
  KFluxSplitVolumeNonLinKernel!(Model.FluxAverage,FV,V,Aux,Metric.dXdxI,Metric.J,DG.DVT,DG.Glob,
    Val(NV),Val(NAUX);ndrange=ndrange)

  group = (DG.OrdPoly+1,20)
  ndrange = (DG.OrdPoly+1,Grid.NumEdges)
  KRiemanNonLinKernel! = RiemanNonLinKernel!(backend,group)
  KRiemanNonLinKernel!(Model.RiemannSolver,FV,V,Aux,DG.Glob,DG.IndE,Grid.EF,Grid.FE,Metric.NH,Metric.T1H,
    Metric.T2H,Metric.VolSurfH,Metric.J,DG.w[1],DG.VZ,Grid,Val(NV),Val(NAUX);ndrange=ndrange) 

  group = ((DG.OrdPoly+1)*(DG.OrdPoly+1),10)
  ndrange = ((DG.OrdPoly+1)*(DG.OrdPoly+1),Grid.NumFaces)
  KVCart2VSpKernel! = VCart2VSpKernel!(backend,group)
  KVCart2VSpKernel!(F,FV,Metric.Rotate;ndrange=ndrange)

  group = ((DG.OrdPoly+1)*(DG.OrdPoly+1),10)
  ndrange = ((DG.OrdPoly+1)*(DG.OrdPoly+1),Grid.NumFaces)
  KSourceKernel! = SourceKernel!(backend,group)
  KSourceKernel!(F,U,Metric.X,Phys;ndrange=ndrange)

# Source!(F,U,Metric,DG,Grid,Phys)
end
