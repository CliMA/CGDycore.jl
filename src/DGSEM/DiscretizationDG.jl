function DiscretizationDG(backend,FT,DG,Exchange,Global,zS,GridType::Grids.Quad)
# Discretization
  Grid = Global.Grid
  OP = DG.DoFE
  OPZ=DG.OrdPolyZ+1
  DoF = DG.DoF
  NumG = DG.NumG
  nz = Grid.nz;
  NF = Grid.NumFaces
  NE = Grid.NumEdges

  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU

  EFCPU = zeros(Int,2,NE)
  FECPU = zeros(Int,2,NE)
  for iE = 1 : NE
    EFCPU[:,iE] = Grid.Edges[iE].F
    FECPU[:,iE] = Grid.Edges[iE].FE
  end
  Grid.EF = KernelAbstractions.zeros(backend,Int,2,NE)
  Grid.FE = KernelAbstractions.zeros(backend,Int,2,NE)
  copyto!(Grid.EF,EFCPU)
  copyto!(Grid.FE,FECPU)


  nQuad = DoF
  FT64 = Float64
  Metric = FiniteElements.Metric!(backend,FT64,DG,Grid,NumberThreadGPU,zS,GridType)
  return Metric
end

function DiscretizationDG(backend,FT,DG,Exchange,Global,zs,GridType::Grids.Tri)
# Discretization
  Grid = Global.Grid
  OP = DG.DoFE
  OPZ=DG.OrdPolyZ+1
  DoF = DG.DoF
  NumG = DG.NumG
  nz = Grid.nz;
  NF = Grid.NumFaces
  NE = Grid.NumEdges

  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU


  nQuad = DoF
  Metric = Grids.MetricStruct{FT}(backend,DoF,OPZ,Global.Grid.NumFaces,nz,NumG)
  F = zeros(3,3,NF)
  FGPU = KernelAbstractions.zeros(backend,FT,3,3,NF)
  for iF = 1 : NF
    F[1,1,iF] = Grid.Faces[iF].P[1].x
    F[1,2,iF] = Grid.Faces[iF].P[1].y
    F[1,3,iF] = Grid.Faces[iF].P[1].z
    F[2,1,iF] = Grid.Faces[iF].P[2].x
    F[2,2,iF] = Grid.Faces[iF].P[2].y
    F[2,3,iF] = Grid.Faces[iF].P[2].z
    F[3,1,iF] = Grid.Faces[iF].P[3].x
    F[3,2,iF] = Grid.Faces[iF].P[3].y
    F[3,3,iF] = Grid.Faces[iF].P[3].z
  end  
  copyto!(FGPU,F)
  Grids.JacobiDG3GPU!(Grid.AdaptGrid,Metric.X,Metric.dXdxI,Metric.J,
    Metric.Rotate,DG,FGPU,Grid.z,zs,Grid.Rad,GridType,Grid.Form)

  Metric.zP = KernelAbstractions.zeros(backend,FT,nz,NumG)
  Metric.dz = KernelAbstractions.zeros(backend,FT,nz,NumG)
  NzG = min(div(NumberThreadGPU,DoF),nz)
  group = (DoF,NzG,1)
  ndrange = (DoF,nz,NF)
  if Global.Grid.Form == Grids.SphericalGrid()
    KGridSizeSphereKernel! = GridSizeSphereDGKernel!(backend,group)
    Rad = Global.Grid.Rad
    KGridSizeSphereKernel!(Metric.dz,Metric.X,DG.Glob,
      Rad,ndrange=ndrange)
  else
    KGridSizeCartKernel! = GridSizeCartDGKernel!(backend,group)
    KGridSizeCartKernel!(Metric.dz,Metric.X,DG.Glob,ndrange=ndrange)
  end

  NFG = min(div(NumberThreadGPU,DoF),NF)
  group = (DoF, NFG)
  ndrange = (DoF, NF)
  if Grid.Form == Grids.SphericalGrid()
    KMetricLowerBoundaryKernel! = MetricLowerBoundaryDGKernel!(backend,group)
    KMetricLowerBoundaryKernel!(Metric.xS,Metric.X,DG.Glob,ndrange=ndrange)
  end

  EFCPU = zeros(Int,2,NE)
  FECPU = zeros(Int,2,NE)
  for iE = 1 : NE
    EFCPU[:,iE] = Grid.Edges[iE].F  
    FECPU[:,iE] = Grid.Edges[iE].FE  
  end  
  Grid.EF = KernelAbstractions.zeros(backend,Int,2,NE)
  Grid.FE = KernelAbstractions.zeros(backend,Int,2,NE)
  copyto!(Grid.EF,EFCPU)
  copyto!(Grid.FE,FECPU)

  NzG = min(div(NumberThreadGPU,OP*OPZ),nz)
  group = (OPZ,OP,NzG)
  ndrange = (OPZ,OP,nz,NE)
  KNormalHTriKernel! = NormalHTriKernel!(backend,group)
  Metric.VolSurfH = KernelAbstractions.zeros(backend,FT,OPZ,OP,nz,NE)
  Metric.NH = KernelAbstractions.zeros(backend,FT,3,OPZ,OP,nz,NE)
  KNormalHTriKernel!(Metric.VolSurfH,Metric.NH,
    Metric.dXdxI,Grid.EF,Grid.FE,DG.PosDoFE,ndrange=ndrange)

  NzG = min(div(NumberThreadGPU,DoF),nz+1)
  group = (OP,OP,NzG)
  ndrange = (DoF,nz+1,NF)
  KNormalVKernel! = NormalVKernel!(backend,group)
  Metric.VolSurfV = KernelAbstractions.zeros(backend,FT,DoF,nz+1,NF)
  Metric.NV = KernelAbstractions.zeros(backend,FT,3,DoF,nz+1,NF)
  KNormalVKernel!(Metric.VolSurfV,Metric.NV,OPZ,Metric.dXdxI,ndrange=ndrange)
  return Metric
end

@kernel inbounds = true function NormalHTriKernel!(VolSurfH,NH,
  @Const(dXdxI),@Const(EF),@Const(FE),@Const(PosDoFE))

# Normal NH(3,I,K,iz,IE)

  K,I,Iz,IE = @index(Global, NTuple)

  NE = @uniform @ndrange()[4]
  NZ = @uniform @ndrange()[3]
  N = @uniform @ndrange()[2]

  if IE <= NE && Iz <= NZ
    IF1 = EF[1,IE]
    IF2 = EF[2,IE]
    if IF1 < IF2
      IF = IF1  
      IS = FE[1,IE]  
    else
      IF = IF2
      IS = FE[2,IE]  
    end  
#   IF = EF[1,IE]
#   IS = FE[1,IE]  
    ID = PosDoFE[I,IS]  
    if IS == 1
      nSLoc1 = -dXdxI[2,1,K,ID,Iz,IF]
      nSLoc2 = -dXdxI[2,2,K,ID,Iz,IF]
      nSLoc3 = -dXdxI[2,3,K,ID,Iz,IF]
    elseif IS == 2
      nSLoc1 = +dXdxI[1,1,K,ID,Iz,IF] + dXdxI[2,1,K,ID,Iz,IF]
      nSLoc2 = +dXdxI[1,2,K,ID,Iz,IF] + dXdxI[2,2,K,ID,Iz,IF]
      nSLoc3 = +dXdxI[1,3,K,ID,Iz,IF] + dXdxI[2,3,K,ID,Iz,IF]
    elseif IS == 3
      nSLoc1 = -dXdxI[1,1,K,ID,Iz,IF]
      nSLoc2 = -dXdxI[1,2,K,ID,Iz,IF]
      nSLoc3 = -dXdxI[1,3,K,ID,Iz,IF]
    end  
    n1Norm = sqrt(nSLoc1 * nSLoc1 + nSLoc2 * nSLoc2 + nSLoc3 * nSLoc3)
    nSLoc1 = nSLoc1 / n1Norm
    nSLoc2 = nSLoc2 / n1Norm
    nSLoc3 = nSLoc3 / n1Norm
    VolSurfH[K,I,Iz,IE] = n1Norm
    if IF1 < IF2
      NH[1,K,I,Iz,IE] = nSLoc1
      NH[2,K,I,Iz,IE] = nSLoc2
      NH[3,K,I,Iz,IE] = nSLoc3
    else  
      NH[1,K,I,Iz,IE] = -nSLoc1
      NH[2,K,I,Iz,IE] = -nSLoc2
      NH[3,K,I,Iz,IE] = -nSLoc3
    end  
  end
end

