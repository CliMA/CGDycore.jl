function DiscretizationDG(backend,FT,Jacobi,DG,Exchange,Global,zs,GridType::Grids.Quad)
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
  Metric = MetricStruct{FT}(backend,DoF,OPZ,Global.Grid.NumFaces,nz,NumG)
  F = zeros(4,3,NF)
  FGPU = KernelAbstractions.zeros(backend,FT,4,3,NF)
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
    F[4,1,iF] = Grid.Faces[iF].P[4].x
    F[4,2,iF] = Grid.Faces[iF].P[4].y
    F[4,3,iF] = Grid.Faces[iF].P[4].z
  end
  copyto!(FGPU,F)
  if Grid.Form == "Sphere"
    Grids.JacobiSphereDG3GPU!(Grid.AdaptGrid,Metric.X,Metric.dXdxI,Metric.J,
      Metric.Rotate,DG,FGPU,Grid.z,zs,Grid.Rad,GridType)
  else
    Grids.JacobiCartDG3GPU!(Grid.AdaptGrid,Metric.X,Metric.dXdxI,Metric.J,
      Metric.Rotate,DG,FGPU,Grid.z,zs,Grid.Rad,GridType)
  end
  Metric.dz = KernelAbstractions.zeros(backend,FT,nz,NumG)
  NzG = min(div(NumberThreadGPU,DoF),nz)
  group = (DoF,NzG,1)
  ndrange = (DoF,nz,NF)
  if Global.Grid.Form == "Sphere"
    KGridSizeSphereDGKernel! = GridSizeSphereDGKernel!(backend,group)
    Rad = Global.Grid.Rad
    KGridSizeSphereDGKernel!(Metric.dz,Metric.X,DG.Glob,
      Rad,ndrange=ndrange)
  else
    KGridSizeCartKernel! = GridSizeCartDGKernel!(backend,group)
    KGridSizeCartKernel!(Metric.dz,Metric.X,DG.Glob,ndrange=ndrange)
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
  KNormalHQuadKernel! = NormalHQuadKernel!(backend,group)
  Metric.VolSurfH = KernelAbstractions.zeros(backend,FT,OPZ,OP,nz,NE)
  Metric.NH = KernelAbstractions.zeros(backend,FT,3,OPZ,OP,nz,NE)
  KNormalHQuadKernel!(Metric.VolSurfH,Metric.NH,
    Metric.dXdxI,Grid.EF,Grid.FE,ndrange=ndrange)

  NzG = min(div(NumberThreadGPU,DoF),nz+1)
  group = (OP,OP,NzG)
  ndrange = (DoF,nz+1,NF)
  KNormalVKernel! = NormalVKernel!(backend,group)
  Metric.VolSurfV = KernelAbstractions.zeros(backend,FT,DoF,nz+1,NF)
  Metric.NV = KernelAbstractions.zeros(backend,FT,3,DoF,nz+1,NF)
  KNormalVKernel!(Metric.VolSurfV,Metric.NV,OPZ,Metric.dXdxI,ndrange=ndrange)
  return Metric
end

@kernel inbounds = true function NormalHQuadKernel!(VolSurfH,NH,
  @Const(dXdxI),@Const(EF),@Const(FE))

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
    if IS == 1
      ID = I  
      nSLoc1 = -dXdxI[2,1,K,ID,Iz,IF]
      nSLoc2 = -dXdxI[2,2,K,ID,Iz,IF]
      nSLoc3 = -dXdxI[2,3,K,ID,Iz,IF]
    elseif IS == 2
      ID = N + (I - 1) * N  
      nSLoc1 = dXdxI[1,1,K,ID,Iz,IF]
      nSLoc2 = dXdxI[1,2,K,ID,Iz,IF]
      nSLoc3 = dXdxI[1,3,K,ID,Iz,IF]
    elseif IS == 3
      ID = I + (N - 1) * N  
      nSLoc1 = -dXdxI[2,1,K,ID,Iz,IF]
      nSLoc2 = -dXdxI[2,2,K,ID,Iz,IF]
      nSLoc3 = -dXdxI[2,3,K,ID,Iz,IF]
    elseif IS == 4
      ID = 1 + (I - 1) * N  
      nSLoc1 = dXdxI[1,1,K,ID,Iz,IF]
      nSLoc2 = dXdxI[1,2,K,ID,Iz,IF]
      nSLoc3 = dXdxI[1,3,K,ID,Iz,IF]
    end  
    n1Norm = sqrt(nSLoc1 * nSLoc1 + nSLoc2 * nSLoc2 + nSLoc3 * nSLoc3)
    nSLoc1 = nSLoc1 / n1Norm
    nSLoc2 = nSLoc2 / n1Norm
    nSLoc3 = nSLoc3 / n1Norm
    VolSurfH[K,I,Iz,IE] = n1Norm
    NH[1,K,I,Iz,IE] = nSLoc1
    NH[2,K,I,Iz,IE] = nSLoc2
    NH[3,K,I,Iz,IE] = nSLoc3
  end
end

@kernel inbounds = true function NormalVKernel!(VolSurfV,NV,M,@Const(dXdxI))

  # Normal NV(3,I,J,2,iz,IF)

  ID,Iz,IF = @index(Global, NTuple)
  NF = @uniform @ndrange()[3]
  NZ = @uniform @ndrange()[2]

  if IF <= NF && Iz <= NZ
    if Iz < NZ  
      nSLoc1 = dXdxI[3,1,1,ID,Iz,IF]
      nSLoc2 = dXdxI[3,2,1,ID,Iz,IF]
      nSLoc3 = dXdxI[3,3,1,ID,Iz,IF]
    else
      nSLoc1 = dXdxI[3,1,M,ID,Iz-1,IF]
      nSLoc2 = dXdxI[3,2,M,ID,Iz-1,IF]
      nSLoc3 = dXdxI[3,3,M,ID,Iz-1,IF]
    end  
    n1Norm = sqrt(nSLoc1 * nSLoc1 + nSLoc2 * nSLoc2 + nSLoc3 * nSLoc3)
    nSLoc1 = nSLoc1 / n1Norm
    nSLoc2 = nSLoc2 / n1Norm
    nSLoc3 = nSLoc3 / n1Norm
    VolSurfV[ID,Iz,IF] = n1Norm
    NV[1,ID,Iz,IF] = nSLoc1
    NV[2,ID,Iz,IF] = nSLoc2
    NV[3,ID,Iz,IF] = nSLoc3
  end
end

function DiscretizationDG(backend,FT,Jacobi,DG,Exchange,Global,zs,GridType::Grids.Tri)
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
  Metric = MetricStruct{FT}(backend,DoF,OPZ,Global.Grid.NumFaces,nz,NumG)
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
  if Grid.Form == "Sphere"
    Grids.JacobiSphereDG3GPU!(Grid.AdaptGrid,Metric.X,Metric.dXdxI,Metric.J,
      Metric.Rotate,DG,FGPU,Grid.z,zs,Grid.Rad,GridType)
  else
    Grids.JacobiCartDG3GPU!(Grid.AdaptGrid,Metric.X,Metric.dXdxI,Metric.J,
      Metric.Rotate,DG,FGPU,Grid.z,zs,Grid.Rad,GridType)
  end

  Metric.dz = KernelAbstractions.zeros(backend,FT,nz,NumG)
  NzG = min(div(NumberThreadGPU,DoF),nz)
  group = (DoF,NzG,1)
  ndrange = (DoF,nz,NF)
  if Global.Grid.Form == "Sphere"
    KGridSizeSphereKernel! = GridSizeSphereKernel!(backend,group)
    Rad = Global.Grid.Rad
    KGridSizeSphereKernel!(Metric.dz,Metric.X,DG.Glob,
      Rad,ndrange=ndrange)
  else
    KGridSizeCartKernel! = GridSizeCartDGKernel!(backend,group)
    KGridSizeCartKernel!(Metric.dz,Metric.X,DG.Glob,ndrange=ndrange)
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

@kernel inbounds = true function GridSizeCartDGKernel!(dz,@Const(X),@Const(Glob))

  ID,Iz,IF = @index(Global, NTuple)

  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  if Iz <= Nz && IF <= NF
    ind = Glob[ID,IF]
    dz[Iz,ind] = X[ID,end,3,Iz,IF] - X[ID,1,3,Iz,IF]
  end
end

@kernel inbounds = true function GridSizeSphereDGKernel!(dz,@Const(X),@Const(Glob),Rad)

  ID,Iz,IF = @index(Global, NTuple)

  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  if Iz <= Nz && IF <= NF
    ind = Glob[ID,IF]
    rS = sqrt(X[ID,1,1,Iz,IF]^2 + X[ID,1,2,Iz,IF]^2 + X[ID,1,3,Iz,IF]^2)
    rE = sqrt(X[ID,end,1,Iz,IF]^2 + X[ID,end,2,Iz,IF]^2 + X[ID,end,3,Iz,IF]^2)
    dz[Iz,ind] = rE - rS
  end
end
