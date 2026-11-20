# Exner pressure Pi = (Rd rho theta / p0)^(Rd/Cvd) for the vector invariant form
function Exner!(Pi,RhoTh,Phys)
  FT = eltype(Pi)
  Rdp0 = FT(Phys.Rd) / FT(Phys.p0)
  ex = FT(Phys.Rd) / FT(Phys.Cvd)
  @. Pi = (Rdp0 * RhoTh)^ex
end

# Pressure in Aux[..,1] and, for the vector-invariant form, the Exner pressure
# in Aux[..,ExnerAuxPos] of the nodes NodeS:NodeE in one pass over rho theta,
# replaces the two broadcasts over views (same formulas as Pressure and Exner!)
function PressureExner!(Aux,U,NodeS,NodeE,Model,Phys,Form,NumberThreadGPU)
  NodeE < NodeS && return
  backend = get_backend(Aux)
  FT = eltype(Aux)
  Rdp0 = FT(Phys.Rd) / FT(Phys.p0)
  ex = FT(Phys.Rd) / FT(Phys.Cvd)
  MNz = size(U,1) * size(U,2)
  UV = reshape(U,:,size(U,4))
  AuxV = reshape(Aux,:,size(Aux,4))
  ExnerPos = Form == VectorInvariantForm() ? ExnerAuxPos : 0
  group = (NumberThreadGPU)
  ndrange = (MNz * (NodeE - NodeS + 1))
  KPressureExnerKernel! = PressureExnerKernel!(backend,group)
  KPressureExnerKernel!(AuxV,UV,Model.Pressure,MNz*(NodeS-1),ExnerPos,Rdp0,ex,ndrange=ndrange)
end

@kernel inbounds = true function PressureExnerKernel!(Aux,@Const(U),Pressure,Off,ExnerPos,Rdp0,ex)
  I, = @index(Global, NTuple)

  N = @uniform @ndrange()[1]

  if I <= N
    RhoTh = U[Off+I,5]
    Aux[Off+I,1] = Pressure(RhoTh)
    if ExnerPos > 0
      Aux[Off+I,ExnerPos] = (Rdp0 * RhoTh)^ex
    end
  end
end

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

  fill!(F,eltype(F)(0))

  Form = Formulation(Model.ModelType)

  PressureExner!(Aux,U,1,DG.NumI,Model,Phys,Form,NumberThreadGPU)

  @views Parallels.ExchangeData3DSendGPU(U[:,:,:,1:NV],Exchange)

  FluxSplitVolumeNonLin(Model.FluxAverageH,Model.FluxAverageV,F,U,Aux,DG,Metric.dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX,Form)
  RiemannNonLinV(Model.RiemannSolver,F,U,Aux,DG,Metric,Grid,NumberThreadGPU,NV,NAUX)

  @views Parallels.ExchangeData3DRecvSetGPU!(U[:,:,:,1:NV],Exchange)
  PressureExner!(Aux,U,DG.NumI+1,DG.NumG,Model,Phys,Form,NumberThreadGPU)

  RiemannNonLinH(Model.RiemannSolver,F,U,Aux,CacheAux,DG,Metric,Grid,NumberThreadGPU,NV,NAUX)

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


  fill!(F,eltype(F)(0))

  @views @. Aux[:,:,1:DG.NumI,1] = Model.Pressure(U[:,:,1:DG.NumI,5])

  @views Parallels.ExchangeData3DSendGPU(U[:,:,:,1:NV],Exchange)

  FluxVolumeNonLinH(Model.FluxAverage,F,U,Aux,DG,Metric.dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX)
  FluxVolumeNonLinV(Model.FluxAverage,F,U,Aux,DG,Metric.dXdxI,Nz,NF,NumberThreadGPU,NV,NAUX)
  @views Parallels.ExchangeData3DRecvSetGPU!(U[:,:,:,1:NV],Exchange)
  @views @. Aux[:,:,DG.NumI+1:DG.NumG,1] = Model.Pressure(U[:,:,DG.NumI+1:DG.NumG,5])

  RiemannNonLinH(Model.RiemannSolver,F,U,Aux,CacheAux,DG,Metric,Grid,NumberThreadGPU,NV,NAUX)
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
