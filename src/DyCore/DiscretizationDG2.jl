function DiscretizationDG2(backend,FT,Jacobi,DG,Exchange,Global)
# Discretization
  Grid = Global.Grid
  OP=DG.OrdPoly+1
  OPZ=DG.OrdPolyZ+1
  DoF = DG.DoF
  NumG = DG.NumG
  nz = Grid.nz;
  NF = Grid.NumFaces
  NE = Grid.NumEdges

  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU


  nQuad = OP * OP
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
  if Global.Grid.Form == "Sphere"
    Grids.JacobiSphereDG2GPU!(Metric.X,Metric.dXdx,Metric.dXdxI,Metric.J,Metric.Rotate,DG,FGPU,Grid.Rad)
  else
    Grids.JacobiCartDG2GPU!(Metric.X,Metric.dXdx,Metric.dXdxI,Metric.J,Metric.Rotate,DG,FGPU,Grid.Rad)
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
  KNormalTangentHKernel! = NormalTangentHKernel!(backend,group)
  Metric.VolSurfH = KernelAbstractions.zeros(backend,FT,OPZ,OP,nz,NE)
  Metric.NH = KernelAbstractions.zeros(backend,FT,3,OPZ,OP,nz,NE)
  Metric.T1H = KernelAbstractions.zeros(backend,FT,3,OPZ,OP,nz,NE)
  Metric.T2H = KernelAbstractions.zeros(backend,FT,3,OPZ,OP,nz,NE)
  KNormalTangentHKernel!(Metric.VolSurfH,Metric.NH,Metric.T1H,Metric.T2H,
    Metric.dXdxI,Grid.EF,Grid.FE,ndrange=ndrange)

  NzG = min(div(NumberThreadGPU,DoF),nz+1)
  group = (OP,OP,NzG)
  ndrange = (DoF,nz+1,NF)
  KNormalTangentVKernel! = NormalTangentVKernel!(backend,group)
  Metric.VolSurfV = KernelAbstractions.zeros(backend,FT,DoF,nz+1,NF)
  Metric.NV = KernelAbstractions.zeros(backend,FT,3,DoF,nz+1,NF)
  Metric.T1V = KernelAbstractions.zeros(backend,FT,3,DoF,nz+1,NF)
  Metric.T2V = KernelAbstractions.zeros(backend,FT,3,DoF,nz+1,NF)
  KNormalTangentVKernel!(Metric.VolSurfV,Metric.NV,Metric.T1V,Metric.T2V,OPZ,Metric.dXdxI,ndrange=ndrange)
  return DG,Metric
end

@kernel inbounds = true function NormalTangentHKernel1!(VolSurfH,NH,T1H,T2H,
  @Const(dXdx),@Const(X),@Const(EF),@Const(FE))

# Normal NH(3,I,K,iz,IE)
# Normal TH1(3,I,K,iz,IE)
# Normal TH2(3,I,K,iz,IE)

  K,I,Iz,IE = @index(Global, NTuple)

  NE = @uniform @ndrange()[4]
  NZ = @uniform @ndrange()[3]
  N = @uniform @ndrange()[2]

  if IE <= NE && Iz <= NZ
    if EF[1,IE] < EF[2,IE]  
      IS = FE[1,IE]  
      IF = EF[1,IE]
    else  
      IS = FE[2,IE]  
      IF = EF[2,IE]
    end  
    if IS == 1
      ID = I  
      tSLoc1 = dXdx[1,1,K,ID,Iz,IF]
      tSLoc2 = dXdx[2,1,K,ID,Iz,IF]
      tSLoc3 = dXdx[3,1,K,ID,Iz,IF]
    elseif IS == 2
      ID = N + (I - 1) * N  
      tSLoc1 = dXdx[1,2,K,ID,Iz,IF]
      tSLoc2 = dXdx[2,2,K,ID,Iz,IF]
      tSLoc3 = dXdx[3,2,K,ID,Iz,IF]
    elseif IS == 3
      ID = I + (N - 1) * N  
      tSLoc1 = dXdx[1,1,K,ID,Iz,IF]
      tSLoc2 = dXdx[2,1,K,ID,Iz,IF]
      tSLoc3 = dXdx[3,1,K,ID,Iz,IF]
    elseif IS == 4
      ID = 1 + (I - 1) * N  
      tSLoc1 = dXdx[1,2,K,ID,Iz,IF]
      tSLoc2 = dXdx[2,2,K,ID,Iz,IF]
      tSLoc3 = dXdx[3,2,K,ID,Iz,IF]
    end  
    XLoc1 = X[ID,1,1,1,IF] 
    XLoc2 = X[ID,1,2,1,IF] 
    XLoc3 = X[ID,1,3,1,IF] 
    nSLoc1 = tSLoc2 * XLoc3 - tSLoc3 * XLoc2 
    nSLoc2 = -(tSLoc1 * XLoc3 - tSLoc3 * XLoc1) 
    nSLoc3 = tSLoc1 * XLoc2 - tSLoc2 * XLoc1 

    n1Norm = sqrt(nSLoc1 * nSLoc1 + nSLoc2 * nSLoc2 + nSLoc3 * nSLoc3)
    nSLoc1 = nSLoc1 / n1Norm
    nSLoc2 = nSLoc2 / n1Norm
    nSLoc3 = nSLoc3 / n1Norm

    t1Norm = sqrt(tSLoc1 * tSLoc1 + tSLoc2 * tSLoc2 + tSLoc3 * tSLoc3)
    tSLoc1 = tSLoc1 / t1Norm
    tSLoc2 = tSLoc2 / t1Norm
    tSLoc3 = tSLoc3 / t1Norm

    VolSurfH[K,I,Iz,IE] = t1Norm

    NH[1,K,I,Iz,IE] = nSLoc1
    NH[2,K,I,Iz,IE] = nSLoc2
    NH[3,K,I,Iz,IE] = nSLoc3
    T1H[1,K,I,Iz,IE] = tSLoc1
    T1H[2,K,I,Iz,IE] = tSLoc2
    T1H[3,K,I,Iz,IE] = tSLoc3
    T2H[1,K,I,Iz,IE] = 0.0
    T2H[2,K,I,Iz,IE] = 0.0
    T2H[3,K,I,Iz,IE] = 0.0
  end
end
