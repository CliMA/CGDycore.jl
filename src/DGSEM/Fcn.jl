function FcnSplit!(F,U,DG,Metric,Phys,CacheAux,Exchange,Global)
  backend = get_backend(F)
  FT = eltype(F)
  Model = Global.Model
  Grid = Global.Grid
  Buo = Model.BuoyancyFun
  Damp = Model.Damp
  Cor = Model.CoriolisFun
  Force = Model.Force
  FT = eltype(F)
  DoF = DG.DoF
  DoFE = DG.DoFE
  M = DG.V.OrdPolyZ + 1
  NF = Grid.NumFaces
  NE = Grid.NumEdges
  Nz = Grid.nz
  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU
  Proc = Global.ParallelCom.Proc
  NV = Model.NumV
  NAUX = Model.NumAux
  @views UI = U[:,:,1:DG.NumI,:]
  @views Aux = CacheAux.Aux

  @. F = 0

  @views @. Aux[:,:,1:DG.NumI,1] = Model.Pressure(U[:,:,1:DG.NumI,5])

  @views Parallels.ExchangeData3DSendGPU(U[:,:,:,1:NV],Exchange)

  FluxSplitVolumeNonLinH(Model.FluxAverageH,F,U,Aux,DG,Metric.dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX)
  FluxSplitVolumeNonLinV(Model.FluxAverageV,F,U,Aux,DG,Metric.dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX)
  RiemannNonLinV(Model.RiemannSolver,F,U,Aux,DG,Metric,Grid,NumberThreadGPU,NV,NAUX)

  @views Parallels.ExchangeData3DRecvSetGPU!(U[:,:,:,1:NV],Exchange)
  @views @. Aux[:,:,DG.NumI+1:DG.NumG,1] = Model.Pressure(U[:,:,DG.NumI+1:DG.NumG,5])

  RiemannNonLinH(Model.RiemannSolver,F,U,Aux,DG,Metric,Grid,NumberThreadGPU,NV,NAUX)

  ScaleMassMatrix!(F,DG,Metric,Grid,NumberThreadGPU,NV)
  if Model.Coriolis
    Sources.Coriolis!(Cor,F,U,Aux,DG,Metric,NumberThreadGPU)
  end

  if Model.Damping
    Sources.Damping!(Damp,F,U,Aux,DG,Metric,NumberThreadGPU)  
  end

  if Model.Forcing
    Sources.Forcing!(Force,F,U,Aux,DG,Metric,NumberThreadGPU)  
  end

  if Model.Dimension == 2
    @. F[:,:,:,3] = 0   
  end  


  if Model.Buoyancy
    Sources.Buoyancy!(Buo,F,U,DG.Glob,Metric.X,NumberThreadGPU)
  end

end

function Fcn!(F,U,DG,Metric,Phys,CacheAux,Exchange,Global)
  backend = get_backend(F)
  Model = Global.Model
  Grid = Global.Grid
  Damp = Model.Damp
  Cor = Model.CoriolisFun
  Force = Model.Force
  FT = eltype(F)
  DoF = DG.DoF
  DoFE = DG.DoFE
  M = DG.V.OrdPolyZ + 1
  NF = Grid.NumFaces
  NE = Grid.NumEdges
  Nz = Grid.nz
  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU
  Proc = Global.ParallelCom.Proc
  NV = Model.NumV
  NAUX = Model.NumAux
  @views UI = U[:,:,1:DG.NumI,:]
  @views Aux = CacheAux.Aux


  @. F = 0

  @views @. Aux[:,:,1:DG.NumI,1] = Model.Pressure(U[:,:,1:DG.NumI,5])

  @views Parallels.ExchangeData3DSendGPU(U[:,:,:,1:NV],Exchange)

  FluxVolumeNonLinH(Model.FluxAverage,F,U,Aux,DG,Metric.dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX)
  FluxVolumeNonLinV(Model.FluxAverage,F,U,Aux,DG,Metric.dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX)
  @views Parallels.ExchangeData3DRecvSetGPU!(U[:,:,:,1:NV],Exchange)
  @views @. Aux[:,:,DG.NumI+1:DG.NumG,1] = Model.Pressure(U[:,:,DG.NumI+1:DG.NumG,5])

  RiemannNonLinH(Model.RiemannSolver,F,U,Aux,DG,Metric,Grid,NumberThreadGPU,NV,NAUX)
  RiemannNonLinV(Model.RiemannSolver,F,U,Aux,DG,Metric,Grid,NumberThreadGPU,NV,NAUX)

  ScaleMassMatrix!(F,DG,Metric,Grid,NumberThreadGPU,NV)

  if Model.Coriolis
    Sources.Coriolis!(Cor,F,U,DG.Glob,Metric.X,NumberThreadGPU)
  end

  if Model.Buoyancy
    Sources.Buoyancy!(Buo,F,U,DG.Glob,Metric.X,NumberThreadGPU)
  end  

  if Model.Damping
    Sources.Damping!(Damp,F,U,DG.Glob,Metric.X,NumberThreadGPU)  
  end

  if Model.Forcing
    Sources.Forcing!(Force,F,U,Aux,DG.Glob,Metric.X,NumberThreadGPU)  
  end


end
