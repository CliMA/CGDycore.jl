function FcnSplit!(F,U,DG,Metric,Phys,CacheAux,Exchange,Global,VelForm)
  backend = get_backend(F)
  Model = Global.Model
  Grid = Global.Grid
  GridType = Grid.Type
  Damp = Model.Damp
  Cor = Model.CoriolisFun
  Force = Model.Force
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
  @views Aux = CacheAux.Aux


  @. F = 0

  @views @. Aux[:,:,1:DG.NumI,1] = Model.Pressure(U[:,:,1:DG.NumI,5])

  @views StateVSp2VCart!(UI[:,:,:,2:4],DG,Metric,NumberThreadGPU,VelForm)  
  @views Parallels.ExchangeData3DSendGPU(U[:,:,:,1:NV],Exchange)

  FluxSplitVolumeNonLinH(Model.FluxAverage,F,U,Aux,DG,Metric.dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX,GridType)
  FluxSplitVolumeNonLinV(Model.FluxAverage,F,U,Aux,DG,Metric.dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX)
  @views Parallels.ExchangeData3DRecvSetGPU!(U[:,:,:,1:NV],Exchange)
  @views @. Aux[:,:,DG.NumI+1:DG.NumG,1] = Model.Pressure(U[:,:,DG.NumI+1:DG.NumG,5])

  RiemannNonLinH(Model.RiemannSolver,F,U,Aux,DG,Metric,Grid,NumberThreadGPU,NV,NAUX)
  RiemannNonLinV(Model.RiemannSolver,F,U,Aux,DG,Metric,Grid,NumberThreadGPU,NV,NAUX)

  ScaleMassMatrix!(F,DG,Metric,Grid,NumberThreadGPU,NV)

  if Model.Coriolis
    Sources.Coriolis!(Cor,F,U,DG.Glob,Metric.X,NumberThreadGPU)
  end

  if Model.Damping
    Sources.Damping!(Damp,F,U,DG.Glob,Metric.X,NumberThreadGPU)  
  end

  if Model.Forcing
    Sources.Forcing!(Force,F,U,Aux,DG.Glob,Metric.X,NumberThreadGPU)  
  end

  @views StateVCart2VSp!(F[:,:,:,2:4],DG,Metric,NumberThreadGPU,VelForm)  

end

function FcnSplitEx!(F,U,DG,Metric,Phys,CacheAux,Exchange,Global,VelForm)
  backend = get_backend(F)
  Model = Global.Model
  Grid = Global.Grid
  GridType = Grid.Type
  Damp = Model.Damp
  Cor = Model.CoriolisFun
  Force = Model.Force
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
  @views Aux = CacheAux.Aux


  @. F = 0

  @views @. Aux[:,:,1:DG.NumI,1] = Model.Pressure(U[:,:,1:DG.NumI,5])

  @views StateVSp2VCart!(UI[:,:,:,2:4],DG,Metric,NumberThreadGPU,VelForm)  
  @views Parallels.ExchangeData3DSendGPU(U[:,:,:,1:NV],Exchange)

  FluxSplitVolumeNonLinH(Model.FluxAverage,F,U,Aux,DG,Metric.dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX,GridType)
  @views Parallels.ExchangeData3DRecvSetGPU!(U[:,:,:,1:NV],Exchange)
  @views @. Aux[:,:,DG.NumI+1:DG.NumG,1] = Model.Pressure(U[:,:,DG.NumI+1:DG.NumG,5])

  RiemannNonLinH(Model.RiemannSolver,F,U,Aux,DG,Metric,Grid,NumberThreadGPU,NV,NAUX)

  ScaleMassMatrix!(F,DG,Metric,Grid,NumberThreadGPU,NV)

  if Model.Coriolis
    Sources.Coriolis!(Cor,F,U,DG.Glob,Metric.X,NumberThreadGPU)
  end

  if Model.Damping
    Sources.Damping!(Damp,F,U,DG.Glob,Metric.X,NumberThreadGPU)  
  end

  if Model.Forcing
    Sources.Forcing!(Force,F,U,Aux,DG.Glob,Metric.X,NumberThreadGPU)  
  end

  @views StateVCart2VSp!(F[:,:,:,2:4],DG,Metric,NumberThreadGPU,VelForm)  

end

function FcnSplitIm!(F,U,DG,Metric,Phys,CacheAux,Exchange,Global,VelForm)
  backend = get_backend(F)
  Model = Global.Model
  Grid = Global.Grid
  GridType = Grid.Type
  Damp = Model.Damp
  Cor = Model.CoriolisFun
  Force = Model.Force
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
  @views Aux = CacheAux.Aux


  @. F = 0

  @views @. Aux[:,:,1:DG.NumI,1] = Model.Pressure(U[:,:,1:DG.NumI,5])

  @views StateVSp2VCart!(UI[:,:,:,2:4],DG,Metric,NumberThreadGPU,VelForm)  

  FluxSplitVolumeNonLinV(Model.FluxAverage,F,U,Aux,DG,Metric.dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX)

  RiemannNonLinV(Model.RiemannSolver,F,U,Aux,DG,Metric,Grid,NumberThreadGPU,NV,NAUX)

  ScaleMassMatrix!(F,DG,Metric,Grid,NumberThreadGPU,NV)

  @views StateVCart2VSp!(F[:,:,:,2:4],DG,Metric,NumberThreadGPU,VelForm)  

end

function FcnSplit1!(F,U,DG,Metric,Phys,CacheAux,Exchange,Global,VelForm)
  @show "FcnSplit1"
  backend = get_backend(F)
  Model = Global.Model
  Grid = Global.Grid
  GridType = Grid.Type
  Damp = Model.Damp
  Cor = Model.CoriolisFun
  Force = Model.Force
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
  @views Aux = CacheAux.Aux


  @. F = 0

  @views @. Aux[:,:,1:DG.NumI,1] = Model.Pressure(U[:,:,1:DG.NumI,5])

  @views StateVSp2VCart!(UI[:,:,:,2:4],DG,Metric,NumberThreadGPU,VelForm)  
  @views Parallels.ExchangeData3DSendGPU(U[:,:,:,1:NV],Exchange)

  FluxSplitVolumeNonLinH1(Model.FluxAverage,F,U,Aux,DG,Metric.dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX,GridType)
  FluxSplitVolumeNonLinV(Model.FluxAverage,F,U,Aux,DG,Metric.dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX)
  @views Parallels.ExchangeData3DRecvSetGPU!(U[:,:,:,1:NV],Exchange)
  @views @. Aux[:,:,DG.NumI+1:DG.NumG,1] = Model.Pressure(U[:,:,DG.NumI+1:DG.NumG,5])

  RiemannNonLinH(Model.RiemannSolver,F,U,Aux,DG,Metric,Grid,NumberThreadGPU,NV,NAUX)
  RiemannNonLinV(Model.RiemannSolver,F,U,Aux,DG,Metric,Grid,NumberThreadGPU,NV,NAUX)

  ScaleMassMatrix!(F,DG,Metric,Grid,NumberThreadGPU,NV)

  if Model.Coriolis
    Sources.Coriolis!(Cor,F,U,DG.Glob,Metric.X,NumberThreadGPU)
  end

  if Model.Damping
    Sources.Damping!(Damp,F,U,DG.Glob,Metric.X,NumberThreadGPU)  
  end

  if Model.Forcing
    Sources.Forcing!(Force,F,U,Aux,DG.Glob,Metric.X,NumberThreadGPU)  
  end

  @views StateVCart2VSp!(F[:,:,:,2:4],DG,Metric,NumberThreadGPU,VelForm)  

end

function Fcn!(F,U,DG,Metric,Phys,CacheAux,Exchange,Global,VelForm)
  backend = get_backend(F)
  Model = Global.Model
  Grid = Global.Grid
  GridType = Grid.Type
  Damp = Model.Damp
  Cor = Model.CoriolisFun
  Force = Model.Force
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
  @views Aux = CacheAux.Aux


  @. F = 0

  @views @. Aux[:,:,1:DG.NumI,1] = Model.Pressure(U[:,:,1:DG.NumI,5])

  @views StateVSp2VCart!(UI[:,:,:,2:4],DG,Metric,NumberThreadGPU,VelForm)  
  @views Parallels.ExchangeData3DSendGPU(U[:,:,:,1:NV],Exchange)

  FluxVolumeNonLinH(Model.FluxAverage,F,U,Aux,DG,Metric.dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX,GridType)
  FluxVolumeNonLinV(Model.FluxAverage,F,U,Aux,DG,Metric.dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX)
  @views Parallels.ExchangeData3DRecvSetGPU!(U[:,:,:,1:NV],Exchange)
  @views @. Aux[:,:,DG.NumI+1:DG.NumG,1] = Model.Pressure(U[:,:,DG.NumI+1:DG.NumG,5])

  RiemannNonLinH(Model.RiemannSolver,F,U,Aux,DG,Metric,Grid,NumberThreadGPU,NV,NAUX)
  RiemannNonLinV(Model.RiemannSolver,F,U,Aux,DG,Metric,Grid,NumberThreadGPU,NV,NAUX)

  ScaleMassMatrix!(F,DG,Metric,Grid,NumberThreadGPU,NV)

  if Model.Coriolis
    Sources.Coriolis!(Cor,F,U,DG.Glob,Metric.X,NumberThreadGPU)
  end

  if Model.Damping
    Sources.Damping!(Damp,F,U,DG.Glob,Metric.X,NumberThreadGPU)  
  end

  if Model.Forcing
    Sources.Forcing!(Force,F,U,Aux,DG.Glob,Metric.X,NumberThreadGPU)  
  end

  @views StateVCart2VSp!(F[:,:,:,2:4],DG,Metric,NumberThreadGPU,VelForm)  

end

function FcnSplitSlow!(F,U,DG,Metric,Phys,CacheAux,Exchange,Global,VelForm)
  backend = get_backend(F)
  Model = Global.Model
  Grid = Global.Grid
  GridType = Grid.Type
  Damp = Model.Damp
  Cor = Model.CoriolisFun
  Force = Model.Force
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
  Aux = CacheAux.Aux

  @. F = 0

  @views @. Aux[:,:,1:DG.NumI,1] = Model.Pressure(U[:,:,1:DG.NumI,5])

  @views StateVSp2VCart!(UI[:,:,:,2:4],DG,Metric,NumberThreadGPU,VelForm)  
  @views Parallels.ExchangeData3DSendGPU(U[:,:,:,1:NV],Exchange)

  FluxSplitVolumeNonLinH(Model.FluxAverageSlow,F,U,Aux,DG,Metric.dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX,GridType)
  FluxSplitVolumeNonLinV(Model.FluxAverageSlow,F,U,Aux,DG,Metric.dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX)
  @views Parallels.ExchangeData3DRecvSetGPU!(U[:,:,:,1:NV],Exchange)
  @views @. Aux[:,:,DG.NumI+1:DG.NumG,1] = Model.Pressure(U[:,:,DG.NumI+1:DG.NumG,5])

  RiemannNonLinH(Model.RiemannSolverSlow,F,U,Aux,DG,Metric,Grid,NumberThreadGPU,NV,NAUX)
  RiemannNonLinV(Model.RiemannSolverSlow,F,U,Aux,DG,Metric,Grid,NumberThreadGPU,NV,NAUX)

  ScaleMassMatrix!(F,DG,Metric,Grid,NumberThreadGPU,NV)

  if Model.Coriolis
    Sources.Coriolis!(Cor,F,U,DG.Glob,Metric.X,NumberThreadGPU)
  end

  if Model.Damping
    Sources.Damping!(Damp,F,U,DG.Glob,Metric.X,NumberThreadGPU)  
  end

  if Model.Forcing
    Sources.Forcing!(Force,F,U,Aux,DG.Glob,Metric.X,NumberThreadGPU)  
  end

  @views StateVCart2VSp!(F[:,:,:,2:4],DG,Metric,NumberThreadGPU,VelForm)  

end

function FcnSplitFast!(F,U,DG,Metric,Phys,CacheAux,Exchange,Global,VelForm)
  backend = get_backend(F)
  Model = Global.Model
  Grid = Global.Grid
  GridType = Grid.Type
  Damp = Model.Damp
  Cor = Model.CoriolisFun
  Force = Model.Force
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
  Aux = CacheAux.Aux
  @views UI = U[:,:,1:DG.NumI,:]


  @. F = 0

  @views @. Aux[:,:,1:DG.NumI,1] = Model.Pressure(U[:,:,1:DG.NumI,5])

  @views StateVSp2VCart!(UI[:,:,:,2:4],DG,Metric,NumberThreadGPU,VelForm)  
  @views Parallels.ExchangeData3DSendGPU(U[:,:,:,1:NV],Exchange)

  FluxSplitVolumeNonLinH(Model.FluxAverageFast,F,U,Aux,DG,Metric.dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX,GridType)
  FluxSplitVolumeNonLinV(Model.FluxAverageFast,F,U,Aux,DG,Metric.dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX)
  @views Parallels.ExchangeData3DRecvSetGPU!(U[:,:,:,1:NV],Exchange)

  @views @. Aux[:,:,DG.NumI+1:DG.NumG,1] = Model.Pressure(U[:,:,DG.NumI+1:DG.NumG,5])

  RiemannNonLinH(Model.RiemannSolverFast,F,U,Aux,DG,Metric,Grid,NumberThreadGPU,NV,NAUX)
  RiemannNonLinV(Model.RiemannSolverFast,F,U,Aux,DG,Metric,Grid,NumberThreadGPU,NV,NAUX)

  ScaleMassMatrix!(F,DG,Metric,Grid,NumberThreadGPU,NV)

  @views StateVCart2VSp!(F[:,:,:,2:4],DG,Metric,NumberThreadGPU,VelForm)  

end

function FcnSplitFastSemi!(F,U,DG,Metric,Phys,CacheAux,Exchange,Global,VelForm)
  backend = get_backend(F)
  Model = Global.Model
  Grid = Global.Grid
  GridType = Grid.Type
  Damp = Model.Damp
  Cor = Model.CoriolisFun
  Force = Model.Force
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
  Aux = CacheAux.Aux
  @views UI = U[:,:,1:DG.NumI,:]


  @. F = 0

  @views StateVSp2VCart!(UI[:,:,:,2:4],DG,Metric,NumberThreadGPU,VelForm)  
  @views Parallels.ExchangeData3DSendGPU(U[:,:,:,1:NV],Exchange)

  FluxSplitVolumeNonLinH(Model.FluxAverageFast,F,U,Aux,DG,Metric.dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX,GridType)
  FluxSplitVolumeNonLinV(Model.FluxAverageFast,F,U,Aux,DG,Metric.dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX)
  @views Parallels.ExchangeData3DRecvSetGPU!(U[:,:,:,1:NV],Exchange)


  RiemannNonLinH(Model.RiemannSolverFast,F,U,Aux,DG,Metric,Grid,NumberThreadGPU,NV,NAUX)
  RiemannNonLinV(Model.RiemannSolverFast,F,U,Aux,DG,Metric,Grid,NumberThreadGPU,NV,NAUX)

  ScaleMassMatrix!(F,DG,Metric,Grid,NumberThreadGPU,NV)

  @views StateVCart2VSp!(F[:,:,:,2:4],DG,Metric,NumberThreadGPU,VelForm)  

end

function FcnFastLin!(F,U,DG,Metric,Phys,CacheAux,Exchange,Global,VelForm)
  backend = get_backend(F)
  Model = Global.Model
  Grid = Global.Grid
  GridType = Grid.Type
  Damp = Model.Damp
  Buo = Model.BuoyancyFun
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
  Aux = CacheAux.Aux
  @views UI = U[:,:,1:DG.NumI,:]


  @. F = 0

  @views StateVSp2VCart!(UI[:,:,:,2:4],DG,Metric,NumberThreadGPU,VelForm)  
  @views Parallels.ExchangeData3DSendGPU(U[:,:,:,1:NV],Exchange)

  FluxVolumeNonLinH(Model.FluxAverageFast,F,U,Aux,DG,Metric.dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX,GridType)
  FluxVolumeNonLinV(Model.FluxAverageFast,F,U,Aux,DG,Metric.dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX)
  @views Parallels.ExchangeData3DRecvSetGPU!(U[:,:,:,1:NV],Exchange)

  RiemannNonLinH(Model.RiemannSolverFast,F,U,Aux,DG,Metric,Grid,NumberThreadGPU,NV,NAUX)
  RiemannNonLinV(Model.RiemannSolverFast,F,U,Aux,DG,Metric,Grid,NumberThreadGPU,NV,NAUX)

  ScaleMassMatrix!(F,DG,Metric,Grid,NumberThreadGPU,NV)

  if Model.Buoyancy
    Sources.Buoyancy!(Buo,F,U,DG.Glob,Metric.X,NumberThreadGPU)
  end

  @views StateVCart2VSp!(F[:,:,:,2:4],DG,Metric,NumberThreadGPU,VelForm)  

end

