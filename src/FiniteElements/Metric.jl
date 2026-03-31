function Metric!(backend,FT,FE::DGElement,Grid,NumberThreadGPU,zS,::Grids.Quad)
  DoF = FE.DoF
  DoFE = FE.DoFE
  N = FE.OrdPoly + 1
  M = FE.OrdPolyZ + 1
  NF = Grid.NumFaces
  NE = Grid.NumEdges
  Nz = Grid.nz
  DW,DS,DV,DVT,B = DG.DerivativeMatrixSingle(FE.OrdPoly)
  DWZ,DSZ,DVZ,DVTZ,BZ = DG.DerivativeMatrixSingle(FE.OrdPolyZ)
  Metric = Grids.MetricStruct{FT}(backend,FE.DoF,FE.OrdPolyZ+1,NF,Nz,FE.NumG)
  F = PointsFromGrid(backend,FT,Grid)
  group = (DoF,M,1,1)
  ndrange = (DoF,M,Nz,NF)
  KFillXKernel! = FillXKernel!(backend,group)
  KFillXKernel!(Metric.X,Grid.Rad,Grid.AdaptGrid,FE.ksi,FE.xwZ,F,Grid.z,zS,Grid.Type,Grid.Form;ndrange=ndrange)
  KRotateKernel! = RotateKernel!(backend,group)
  KRotateKernel!(Metric.Rotate,Metric.X,Grid.Form;ndrange=ndrange)
  group = (N,N,M,1,1)
  ndrange = (N,N,M,Nz,NF)
  CurlMetric = false
  if CurlMetric
  KMetricKernel1! = MetricCurlQuadKernel1!(backend,group)
  KMetricKernel1!(Metric.dXdxI,Metric.X,FE.DS,FE.DSZ,Val(N),Val(M);ndrange=ndrange)
  KMetricKernel2! = MetricCurlQuadKernel2!(backend,group)
  KMetricKernel2!(Metric.dXdxI,Metric.X,FE.DS,FE.DSZ,Val(N),Val(M);ndrange=ndrange)
  KMetricKernel3! = MetricCurlQuadKernel3!(backend,group)
  KMetricKernel3!(Metric.dXdxI,Metric.X,FE.DS,Val(N),Val(M);ndrange=ndrange)
  else
    @show "NonCurlMetric"  
    KMetricKernel! = MetricQuadKernel!(backend,group)
    KMetricKernel!(Metric.dXdxI,Metric.X,DS,DSZ,Val(N),Val(M);ndrange=ndrange)
  end  
  group = (DoF,M,1,1)
  ndrange = (DoF,M,Nz,NF)
  KDetKernel! = DetKernel!(backend,group)
  KDetKernel!(Metric.J,Metric.dXdxI;ndrange=ndrange)
  Metric.zP = KernelAbstractions.zeros(backend,FT,Nz,FE.NumG)
  Metric.dz = KernelAbstractions.zeros(backend,FT,Nz,FE.NumG)

  GridSizeDGKernel!(FE,Metric,Grid.Rad,NumberThreadGPU,Grid.Form)

  NFG = min(div(NumberThreadGPU,DoF),NF)
  group = (DoF, NFG)
  ndrange = (DoF, NF)
  if Grid.Form == Grids.SphericalGrid()
    KMetricLowerBoundaryKernel! = MetricLowerBoundaryDGKernel!(backend,group)
    KMetricLowerBoundaryKernel!(Metric.xS,Metric.X,FE.Glob,ndrange=ndrange)
  end

  NzG = min(div(NumberThreadGPU,N*M),Nz)
  group = (M,N,NzG)
  ndrange = (M,N,Nz,NE)
  KNormalHQuadKernel! = NormalHQuadKernel!(backend,group)
  Metric.VolSurfH = KernelAbstractions.zeros(backend,FT,M,N,Nz,NE)
  Metric.NH = KernelAbstractions.zeros(backend,FT,3,M,N,Nz,NE)
  KNormalHQuadKernel!(Metric.VolSurfH,Metric.NH,
    Metric.dXdxI,Grid.EF,Grid.FE,ndrange=ndrange)

  NzG = min(div(NumberThreadGPU,DoF),Nz+1)
  group = (N,N,NzG)
  ndrange = (DoF,Nz+1,NF)
  KNormalVKernel! = NormalVKernel!(backend,group)
  Metric.VolSurfV = KernelAbstractions.zeros(backend,FT,DoF,Nz+1,NF)
  Metric.NV = KernelAbstractions.zeros(backend,FT,3,DoF,Nz+1,NF)
  KNormalVKernel!(Metric.VolSurfV,Metric.NV,M,Metric.dXdxI,ndrange=ndrange)

  return Metric
end

function Metric!(backend,FT,FE::DGElement,Grid,NumberThreadGPU,zS,::Grids.Tri)
  DoF = FE.DoF
  DoFE = FE.DoFE
  M = FE.OrdPolyZ + 1
  NF = Grid.NumFaces
  Nz = Grid.nz
  Metric = Grids.MetricStruct{FT}(backend,FE.DoF,FE.OrdPolyZ+1,NF,Nz,FE.NumG)
  F = PointsFromGrid(backend,FT,Grid)
  group = (DoF,M,1,1)
  ndrange = (DoF,M,Nz,NF)
  KFillXKernel! = FillXKernel!(backend,group)
  KFillXKernel!(Metric.X,Grid.Rad,Grid.AdaptGrid,FE.ksi,FE.xwZ,F,Grid.z,zS,Grid.Type,Grid.Form;ndrange=ndrange)
  KRotateKernel! = RotateKernel!(backend,group)
  KRotateKernel!(Metric.Rotate,Metric.X,Grid.Form;ndrange=ndrange)
  group = (DoF,M,1,1)
  ndrange = (DoF,M,Nz,NF)
  @show "KMetricKernel1"
  KMetricKernel1! = MetricTriKernel1!(backend,group)
  KMetricKernel1!(Metric.dXdxI,Metric.J,Metric.X,FE.Dx2,FE.DSZ,Val(DoF),Val(M);ndrange=ndrange)
  @show "KMetricKernel2"
  KMetricKernel2! = MetricTriKernel2!(backend,group)
  KMetricKernel2!(Metric.dXdxI,Metric.J,Metric.X,FE.Dx1,FE.DSZ,Val(DoF),Val(M);ndrange=ndrange)
  @show "KMetricKernel3"
  KMetricKernel3! = MetricTriKernel3!(backend,group)
  KMetricKernel3!(Metric.dXdxI,Metric.J,Metric.X,FE.Dx1,FE.Dx2,Val(DoF),Val(M);ndrange=ndrange)
  group = (DoF,M,1,1)
  ndrange = (DoF,M,Nz,NF)
  KDetKernel! = DetKernel!(backend,group)
  KDetKernel!(Metric.J,Metric.dXdxI;ndrange=ndrange)
  return Metric
end


function Metric!(backend,FT,FE::CGElement,Grid,NumberThreadGPU,zS)
  DoF = FE.DoF
  DoFE = FE.DoFE
  N = FE.OrdPoly + 1
  M = FE.OrdPolyZ + 1
  NF = Grid.NumFaces
  Nz = Grid.nz
  Metric = Grids.MetricStruct{FT}(backend,FE.DoF,FE.OrdPolyZ+1,NF,Nz,FE.NumG)
  F = PointsFromGrid(backend,FT,Grid)
  group = (DoF,M,1,1)
  ndrange = (DoF,M,Nz,NF)
  KFillXKernel! = FillXKernel!(backend,group)
  KFillXKernel!(Metric.X,Grid.AdaptGrid,FE.ksi,FE.xwZ,F,Grid.z,zS,Grid.Type,Grid.Form;ndrange=ndrange)
  group = (N,N,M,1,1)
  ndrange = (N,N,M,Nz,NF)
  KMetricKernel1! = MetricKernel1!(backend,group)
  KMetricKernel1!(Metric.dXdxI,Metric.J,Metric.X,FE.DS,FE.DSZ,Val(N),Val(M);ndrange=ndrange)
  KMetricKernel2! = MetricKernel2!(backend,group)
  KMetricKernel2!(Metric.dXdxI,Metric.J,Metric.X,FE.DS,FE.DSZ,Val(N),Val(M);ndrange=ndrange)
  KMetricKernel3! = MetricKernel3!(backend,group)
  KMetricKernel3!(Metric.dXdxI,Metric.J,Metric.X,FE.DS,FE.DSZ,Val(N),Val(M);ndrange=ndrange)
  return Metric
end

@kernel inbounds = true function MetricCurlQuadKernel3!(dXdxI,@Const(X),@Const(DH), ::Val{N}, ::Val{M}) where {N,M}

  I, J, K,   = @index(Local, NTuple)
  _,_,_,IZ,IF = @index(Global, NTuple)

  XLoc = @localmem eltype(dXdxI) (N,N,M,3)
  DXdx = @localmem eltype(dXdxI) (N,N,M)
  DXdy = @localmem eltype(dXdxI) (N,N,M)
  DYdx = @localmem eltype(dXdxI) (N,N,M)
  DYdy = @localmem eltype(dXdxI) (N,N,M)
  DZdx = @localmem eltype(dXdxI) (N,N,M)
  DZdy = @localmem eltype(dXdxI) (N,N,M)
  NZ = @uniform @ndrange()[4]
  NF = @uniform @ndrange()[5]


  if IZ <= NZ && IF <= NF
    ID = I + (J - 1) * M
    XLoc[I,J,K,1] = X[ID,K,1,IZ,IF]
    XLoc[I,J,K,2] = X[ID,K,2,IZ,IF]
    XLoc[I,J,K,3] = X[ID,K,3,IZ,IF]
  end  
  @synchronize

  if IZ <= NZ && IF <= NF
    DXdx[I,J,K] = DH[I,1] * XLoc[1,J,K,1]
    DYdx[I,J,K] = DH[I,1] * XLoc[1,J,K,2]
    DZdx[I,J,K] = DH[I,1] * XLoc[1,J,K,3]
    DXdy[I,J,K] = DH[J,1] * XLoc[I,1,K,1]
    DYdy[I,J,K] = DH[J,1] * XLoc[I,1,K,2]
    DZdy[I,J,K] = DH[J,1] * XLoc[I,1,K,3]
    @unroll for l = 2 : N
      DXdx[I,J,K] += DH[I,l] * XLoc[l,J,K,1]
      DYdx[I,J,K] += DH[I,l] * XLoc[l,J,K,2]
      DZdx[I,J,K] += DH[I,l] * XLoc[l,J,K,3]
      DXdy[I,J,K] += DH[J,l] * XLoc[I,l,K,1]
      DYdy[I,J,K] += DH[J,l] * XLoc[I,l,K,2]
      DZdy[I,J,K] += DH[J,l] * XLoc[I,l,K,3]
    end
  end  
  @synchronize

  if IZ <= NZ && IF <= NF
    ID = I + (J - 1) * M
    t1 = DH[I,1] * (XLoc[1,J,K,2] * DZdy[1,J,K] - XLoc[1,J,K,3] * DYdy[1,J,K]) -
      DH[J,1] * (XLoc[I,1,K,2] * DZdx[I,1,K] - XLoc[I,1,K,3] * DYdx[I,1,K])

    t2 = DH[I,1] * (XLoc[1,J,K,3] * DXdy[1,J,K] - XLoc[1,J,K,1] * DZdy[1,J,K]) -
      DH[J,1] * (XLoc[I,1,K,3] * DXdx[I,1,K] - XLoc[I,1,K,1] * DZdx[I,1,K])

    t3 = DH[I,1] * (XLoc[1,J,K,1] * DYdy[1,J,K] - XLoc[1,J,K,2] * DXdy[1,J,K]) -
      DH[J,1] * (XLoc[I,1,K,1] * DYdx[I,1,K] - XLoc[I,1,K,2] * DXdx[I,1,K])
    @unroll for l = 2 : N
      t1 += DH[I,l] * (XLoc[l,J,K,2] * DZdy[l,J,K] - XLoc[l,J,K,3] * DYdy[l,J,K]) -
        DH[J,l] * (XLoc[I,l,K,2] * DZdx[I,l,K] - XLoc[I,l,K,3] * DYdx[I,l,K])

      t2 += DH[I,l] * (XLoc[l,J,K,3] * DXdy[l,J,K] - XLoc[l,J,K,1] * DZdy[l,J,K]) -
        DH[J,l] * (XLoc[I,l,K,3] * DXdx[I,l,K] - XLoc[I,l,K,1] * DZdx[I,l,K])
        
      t3 += DH[I,l] * (XLoc[l,J,K,1] * DYdy[l,J,K] - XLoc[l,J,K,2] * DXdy[l,J,K]) -
        DH[J,l] * (XLoc[I,l,K,1] * DYdx[I,l,K] - XLoc[I,l,K,2] * DXdx[I,l,K])
    end  
    dXdxI[3,1,K,ID,IZ,IF] = 0.5 * t1
    dXdxI[3,2,K,ID,IZ,IF] = 0.5 * t2
    dXdxI[3,3,K,ID,IZ,IF] = 0.5 * t3
  end
end

@kernel inbounds = true function MetricQuadKernel!(dXdxI,@Const(X),@Const(DH),@Const(DV), ::Val{N}, ::Val{M}) where {N,M}

  I, J, K,   = @index(Local, NTuple)
  _,_,_,IZ,IF = @index(Global, NTuple)

  XLoc = @localmem eltype(dXdxI) (N,N,M,3)
  NZ = @uniform @ndrange()[4]
  NF = @uniform @ndrange()[5]


  if IZ <= NZ && IF <= NF
    ID = I + (J - 1) * M
    XLoc[I,J,K,1] = X[ID,K,1,IZ,IF]
    XLoc[I,J,K,2] = X[ID,K,2,IZ,IF]
    XLoc[I,J,K,3] = X[ID,K,3,IZ,IF]
  end  
  @synchronize

  if IZ <= NZ && IF <= NF
    DXdx = DH[I,1] * XLoc[1,J,K,1]
    DYdx = DH[I,1] * XLoc[1,J,K,2]
    DZdx = DH[I,1] * XLoc[1,J,K,3]
    DXdy = DH[J,1] * XLoc[I,1,K,1]
    DYdy = DH[J,1] * XLoc[I,1,K,2]
    DZdy = DH[J,1] * XLoc[I,1,K,3]
    @unroll for l = 2 : N
      DXdx += DH[I,l] * XLoc[l,J,K,1]
      DYdx += DH[I,l] * XLoc[l,J,K,2]
      DZdx += DH[I,l] * XLoc[l,J,K,3]
      DXdy += DH[J,l] * XLoc[I,l,K,1]
      DYdy += DH[J,l] * XLoc[I,l,K,2]
      DZdy += DH[J,l] * XLoc[I,l,K,3]
    end
    DXdz = DV[K,1] * XLoc[I,J,1,1]
    DYdz = DV[K,1] * XLoc[I,J,1,2]
    DZdz = DV[K,1] * XLoc[I,J,1,3]
    @unroll for l = 2 : M
      DXdz += DV[K,l] * XLoc[I,J,l,1]
      DYdz += DV[K,l] * XLoc[I,J,l,2]
      DZdz += DV[K,l] * XLoc[I,J,l,3]
    end

    ID = I + (J - 1) * M

    t11 = DH[J,1] * (XLoc[I,1,K,2] * DZdz[I,1,K] - XLoc[I,1,K,3] * DYdz[I,1,K])
    t21 = DH[J,1] * (XLoc[I,1,K,3] * DXdz[I,1,K] - XLoc[I,1,K,1] * DZdz[I,1,K])
    t31 = DH[J,1] * (XLoc[I,1,K,1] * DYdz[I,1,K] - XLoc[I,1,K,2] * DXdz[I,1,K])
    dXdxI[1,1,K,ID,IZ,IF] = DYdy * DZdz - DZdy * DYdz
    dXdxI[1,2,K,ID,IZ,IF] = DZdy * DXdz - DXdy * DZdz
    dXdxI[1,3,K,ID,IZ,IF] = DXdy * DYdz - DYdy * DXdz 

    dXdxI[2,1,K,ID,IZ,IF] = DYdz * DZdx - DZdz * DYdx
    dXdxI[2,2,K,ID,IZ,IF] = DZdz * DXdx - DXdz * DZdx
    dXdxI[2,3,K,ID,IZ,IF] = DXdz * DYdx - DYdz * DXdx 

    dXdxI[3,1,K,ID,IZ,IF] = DYdx * DZdy - DZdx * DYdy
    dXdxI[3,2,K,ID,IZ,IF] = DZdx * DXdy - DXdx * DZdy
    dXdxI[3,3,K,ID,IZ,IF] = DXdx * DYdy - DYdx * DXdy
  end

end

@kernel inbounds = true function MetricTriKernel3!(dXdxI,J,@Const(X),@Const(Dx1),
  @Const(Dx2), ::Val{DoF}, ::Val{M}) where {DoF,M}

  ID, K,   = @index(Local, NTuple)
  _,_,IZ,IF = @index(Global, NTuple)

  XLoc = @localmem eltype(dXdxI) (DoF,M,3)
  DXdx = @localmem eltype(dXdxI) (DoF,M)
  DXdy = @localmem eltype(dXdxI) (DoF,M)
  DYdx = @localmem eltype(dXdxI) (DoF,M)
  DYdy = @localmem eltype(dXdxI) (DoF,M)
  DZdx = @localmem eltype(dXdxI) (DoF,M)
  DZdy = @localmem eltype(dXdxI) (DoF,M)
  NZ = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]


  if IZ <= NZ && IF <= NF
    XLoc[ID,K,1] = X[ID,K,1,IZ,IF]
    XLoc[ID,K,2] = X[ID,K,2,IZ,IF]
    XLoc[ID,K,3] = X[ID,K,3,IZ,IF]
  end  
  @synchronize

  if IZ <= NZ && IF <= NF
    DXdx[ID,K] = Dx1[ID,1] * XLoc[1,K,1]
    DXdy[ID,K] = Dx2[ID,1] * XLoc[1,K,1]
    DYdx[ID,K] = Dx1[ID,1] * XLoc[1,K,2]
    DYdy[ID,K] = Dx2[ID,1] * XLoc[1,K,2]
    DZdx[ID,K] = Dx1[ID,1] * XLoc[1,K,3]
    DZdy[ID,K] = Dx2[ID,1] * XLoc[1,K,3]
    @unroll for l = 2 : DoF
      DXdx[ID,K] += Dx1[ID,l] * XLoc[l,K,1]
      DXdy[ID,K] += Dx2[ID,l] * XLoc[l,K,1]
      DYdx[ID,K] += Dx1[ID,l] * XLoc[l,K,2]
      DYdy[ID,K] += Dx2[ID,l] * XLoc[l,K,2]
      DZdx[ID,K] += Dx1[ID,l] * XLoc[l,K,3]
      DZdy[ID,K] += Dx2[ID,l] * XLoc[l,K,3]
    end
  end  
  @synchronize

  if IZ <= NZ && IF <= NF
    t1 = Dx1[ID,1] * (XLoc[1,K,2] * DZdy[1,K] - XLoc[1,K,3] * DYdy[1,K]) -
      Dx2[ID,1] * (XLoc[1,K,2] * DZdx[1,K] - XLoc[1,K,3] * DYdx[1,K])
    t2 = Dx1[ID,1] * (XLoc[1,K,3] * DXdy[1,K] - XLoc[1,K,1] * DZdy[1,K]) -
      Dx2[ID,1] * (XLoc[1,K,3] * DXdx[1,K] - XLoc[1,K,1] * DZdx[1,K])
    t3 = Dx1[ID,1] * (XLoc[1,K,1] * DYdy[1,K] - XLoc[1,K,2] * DXdy[1,K]) -
      Dx2[ID,1] * (XLoc[1,K,1] * DYdx[1,K] - XLoc[1,K,2] * DXdx[1,K])
    @unroll for l = 2 : DoF
      t1 += Dx1[ID,l] * (XLoc[l,K,2] * DZdy[l,K] - XLoc[l,K,3] * DYdy[l,K]) -
        Dx2[ID,l] * (XLoc[l,K,2] * DZdx[l,K] - XLoc[l,K,3] * DYdx[l,K])
      t2 += Dx1[ID,l] * (XLoc[l,K,3] * DXdy[l,K] - XLoc[l,K,1] * DZdy[l,K]) -
        Dx2[ID,l] * (XLoc[l,K,3] * DXdx[l,K] - XLoc[l,K,1] * DZdx[l,K])
      t3 += Dx1[ID,l] * (XLoc[l,K,1] * DYdy[l,K] - XLoc[l,K,2] * DXdy[l,K]) -
        Dx2[ID,l] * (XLoc[l,K,1] * DYdx[l,K] - XLoc[l,K,2] * DXdx[l,K])
    end  
    dXdxI[3,1,K,ID,IZ,IF] = 0.5 * t1
    dXdxI[3,2,K,ID,IZ,IF] = 0.5 * t2
    dXdxI[3,3,K,ID,IZ,IF] = 0.5 * t3
  end
end

@kernel inbounds = true function MetricCurlQuadKernel2!(dXdxI,@Const(X),@Const(DH),@Const(DV), ::Val{N}, ::Val{M}) where {N,M}

  I, J, K,   = @index(Local, NTuple)
  _,_,_,IZ,IF = @index(Global, NTuple)

  XLoc = @localmem eltype(dXdxI) (N,N,M,3)
  DXdx = @localmem eltype(dXdxI) (N,N,M)
  DXdz = @localmem eltype(dXdxI) (N,N,M)
  DYdx = @localmem eltype(dXdxI) (N,N,M)
  DYdz = @localmem eltype(dXdxI) (N,N,M)
  DZdx = @localmem eltype(dXdxI) (N,N,M)
  DZdz = @localmem eltype(dXdxI) (N,N,M)
  NZ = @uniform @ndrange()[4]
  NF = @uniform @ndrange()[5]


  if IZ <= NZ && IF <= NF
    ID = I + (J - 1) * M
    XLoc[I,J,K,1] = X[ID,K,1,IZ,IF]
    XLoc[I,J,K,2] = X[ID,K,2,IZ,IF]
    XLoc[I,J,K,3] = X[ID,K,3,IZ,IF]
  end  
  @synchronize

  if IZ <= NZ && IF <= NF
    DXdx[I,J,K] = DH[I,1] * XLoc[1,J,K,1]
    DYdx[I,J,K] = DH[I,1] * XLoc[1,J,K,2]
    DZdx[I,J,K] = DH[I,1] * XLoc[1,J,K,3]
    @unroll for l = 2 : N
      DXdx[I,J,K] += DH[I,l] * XLoc[l,J,K,1]
      DYdx[I,J,K] += DH[I,l] * XLoc[l,J,K,2]
      DZdx[I,J,K] += DH[I,l] * XLoc[l,J,K,3]
    end
    DXdz[I,J,K] = DV[K,1] * XLoc[I,J,1,1]
    DYdz[I,J,K] = DV[K,1] * XLoc[I,J,1,2]
    DZdz[I,J,K] = DV[K,1] * XLoc[I,J,1,3]
    @unroll for l = 2 : M
      DXdz[I,J,K] += DV[K,l] * XLoc[I,J,l,1]
      DYdz[I,J,K] += DV[K,l] * XLoc[I,J,l,2]
      DZdz[I,J,K] += DV[K,l] * XLoc[I,J,l,3]
    end  
  end  
  @synchronize

  if IZ <= NZ && IF <= NF
    ID = I + (J - 1) * M
    t11 = DV[K,1] * (XLoc[I,J,1,2] * DZdx[I,J,1] - XLoc[I,J,1,3] * DYdx[I,J,1]) 
    t21 = DV[K,1] * (XLoc[I,J,1,3] * DXdx[I,J,1] - XLoc[I,J,1,1] * DZdx[I,J,1]) 
    t31 = DV[K,1] * (XLoc[I,J,1,1] * DYdx[I,J,1] - XLoc[I,J,1,2] * DXdx[I,J,1]) 
    @unroll for l = 2 : M
      t11 += DV[K,l] * (XLoc[I,J,l,2] * DZdx[I,J,l] - XLoc[I,J,l,3] * DYdx[I,J,l]) 
      t21 += DV[K,l] * (XLoc[I,J,l,3] * DXdx[I,J,l] - XLoc[I,J,l,1] * DZdx[I,J,l]) 
      t31 += DV[K,l] * (XLoc[I,J,l,1] * DYdx[I,J,l] - XLoc[I,J,l,2] * DXdx[I,J,l]) 
    end  

    t12 = DH[I,1] * (XLoc[1,J,K,2] * DZdz[1,J,K] - XLoc[1,J,K,3] * DYdz[1,J,K])
    t22 = DH[I,1] * (XLoc[1,J,K,3] * DXdz[1,J,K] - XLoc[1,J,K,1] * DZdz[1,J,K])
    t32 = DH[I,1] * (XLoc[1,J,K,1] * DYdz[1,J,K] - XLoc[1,J,K,2] * DXdz[1,J,K])
    @unroll for l = 2 : N
      t12 += DH[I,l] * (XLoc[l,J,K,2] * DZdz[l,J,K] - XLoc[l,J,K,3] * DYdz[l,J,K])
      t22 += DH[I,l] * (XLoc[l,J,K,3] * DXdz[l,J,K] - XLoc[l,J,K,1] * DZdz[l,J,K])
      t32 += DH[I,l] * (XLoc[l,J,K,1] * DYdz[l,J,K] - XLoc[l,J,K,2] * DXdz[l,J,K])
    end  

    dXdxI[2,1,K,ID,IZ,IF] = 0.5 * (t11 - t12)
    dXdxI[2,2,K,ID,IZ,IF] = 0.5 * (t21 - t22)
    dXdxI[2,3,K,ID,IZ,IF] = 0.5 * (t31 - t32)
  end
end

@kernel inbounds = true function MetricTriKernel2!(dXdxI,J,@Const(X),@Const(Dx1),
  @Const(DV), ::Val{DoF}, ::Val{M}) where {DoF,M}

  ID, K,   = @index(Local, NTuple)
  _,_,IZ,IF = @index(Global, NTuple)

  XLoc = @localmem eltype(dXdxI) (DoF,M,3)
  DXdx = @localmem eltype(dXdxI) (DoF,M)
  DXdz = @localmem eltype(dXdxI) (DoF,M)
  DYdx = @localmem eltype(dXdxI) (DoF,M)
  DYdz = @localmem eltype(dXdxI) (DoF,M)
  DZdx = @localmem eltype(dXdxI) (DoF,M)
  DZdz = @localmem eltype(dXdxI) (DoF,M)
  NZ = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]


  if IZ <= NZ && IF <= NF
    XLoc[ID,K,1] = X[ID,K,1,IZ,IF]
    XLoc[ID,K,2] = X[ID,K,2,IZ,IF]
    XLoc[ID,K,3] = X[ID,K,3,IZ,IF]
  end  
  @synchronize

  if IZ <= NZ && IF <= NF
    DXdx[ID,K] = Dx1[ID,1] * XLoc[1,K,1]
    DYdx[ID,K] = Dx1[ID,1] * XLoc[1,K,2]
    DZdx[ID,K] = Dx1[ID,1] * XLoc[1,K,3]
    @unroll for l = 2 : DoF 
      DXdx[ID,K] += Dx1[ID,l] * XLoc[l,K,1]
      DYdx[ID,K] += Dx1[ID,l] * XLoc[l,K,2]
      DZdx[ID,K] += Dx1[ID,l] * XLoc[l,K,3]
    end 
    DXdz[ID,K] = DV[K,1] * XLoc[ID,1,1]
    DYdz[ID,K] = DV[K,1] * XLoc[ID,1,2]
    DZdz[ID,K] = DV[K,1] * XLoc[ID,1,3]
    @unroll for l = 2 : M
      DXdz[ID,K] += DV[K,l] * XLoc[ID,l,1]
      DYdz[ID,K] += DV[K,l] * XLoc[ID,l,2]
      DZdz[ID,K] += DV[K,l] * XLoc[ID,l,3]
    end  
  end  
  @synchronize

  if IZ <= NZ && IF <= NF
    t11 = DV[K,1] * (XLoc[ID,1,2] * DZdx[ID,1] - XLoc[ID,1,3] * DYdx[ID,1]) 
    t21 = DV[K,1] * (XLoc[ID,1,3] * DXdx[ID,1] - XLoc[ID,1,1] * DZdx[ID,1]) 
    t31 = DV[K,1] * (XLoc[ID,1,1] * DYdx[ID,1] - XLoc[ID,1,2] * DXdx[ID,1]) 

    @unroll for l = 2 : M
      t11 += DV[K,l] * (XLoc[ID,l,2] * DZdx[ID,l] - XLoc[ID,l,3] * DYdx[ID,l]) 
      t21 += DV[K,l] * (XLoc[ID,l,3] * DXdx[ID,l] - XLoc[ID,l,1] * DZdx[ID,l]) 
      t31 += DV[K,l] * (XLoc[ID,l,1] * DYdx[ID,l] - XLoc[ID,l,2] * DXdx[ID,l]) 
    end  

    t12 = Dx1[ID,1] * (XLoc[1,K,2] * DZdz[1,K] - XLoc[1,K,3] * DYdz[1,K])
    t22 = Dx1[ID,1] * (XLoc[1,K,3] * DXdz[1,K] - XLoc[1,K,1] * DZdz[1,K])
    t32 = Dx1[ID,1] * (XLoc[1,K,1] * DYdz[1,K] - XLoc[1,K,2] * DXdz[1,K])

    @unroll for l = 2 : DoF
      t12 += Dx1[ID,l] * (XLoc[l,K,2] * DZdz[l,K] - XLoc[l,K,3] * DYdz[l,K])
      t22 += Dx1[ID,l] * (XLoc[l,K,3] * DXdz[l,K] - XLoc[l,K,1] * DZdz[l,K])
      t32 += Dx1[ID,l] * (XLoc[l,K,1] * DYdz[l,K] - XLoc[l,K,2] * DXdz[l,K])
    end  

    dXdxI[2,1,K,ID,IZ,IF] = 0.5 * (t11 - t12)
    dXdxI[2,2,K,ID,IZ,IF] = 0.5 * (t21 - t22)
    dXdxI[2,3,K,ID,IZ,IF] = 0.5 * (t31 - t32)
  end
end

@kernel inbounds = true function MetricCurlQuadKernel1!(dXdxI,@Const(X),@Const(DH),@Const(DV), ::Val{N}, ::Val{M}) where {N,M}

  I, J, K,   = @index(Local, NTuple)
  _,_,_,IZ,IF = @index(Global, NTuple)

  XLoc = @localmem eltype(dXdxI) (N,N,M,3)
  DXdy = @localmem eltype(dXdxI) (N,N,M)
  DXdz = @localmem eltype(dXdxI) (N,N,M)
  DYdy = @localmem eltype(dXdxI) (N,N,M)
  DYdz = @localmem eltype(dXdxI) (N,N,M)
  DZdy = @localmem eltype(dXdxI) (N,N,M)
  DZdz = @localmem eltype(dXdxI) (N,N,M)
  NZ = @uniform @ndrange()[4]
  NF = @uniform @ndrange()[5]


  if IZ <= NZ && IF <= NF
    ID = I + (J - 1) * M
    XLoc[I,J,K,1] = X[ID,K,1,IZ,IF]
    XLoc[I,J,K,2] = X[ID,K,2,IZ,IF]
    XLoc[I,J,K,3] = X[ID,K,3,IZ,IF]
  end  
  @synchronize

  if IZ <= NZ && IF <= NF
    DXdy[I,J,K] = DH[J,1] * XLoc[I,1,K,1]
    DYdy[I,J,K] = DH[J,1] * XLoc[I,1,K,2]
    DZdy[I,J,K] = DH[J,1] * XLoc[I,1,K,3]
    @unroll for l = 2 : N
      DXdy[I,J,K] += DH[J,l] * XLoc[I,l,K,1]
      DYdy[I,J,K] += DH[J,l] * XLoc[I,l,K,2]
      DZdy[I,J,K] += DH[J,l] * XLoc[I,l,K,3]
    end
    DXdz[I,J,K] = DV[K,1] * XLoc[I,J,1,1]
    DYdz[I,J,K] = DV[K,1] * XLoc[I,J,1,2]
    DZdz[I,J,K] = DV[K,1] * XLoc[I,J,1,3]
    @unroll for l = 2 : M
      DXdz[I,J,K] += DV[K,l] * XLoc[I,J,l,1]
      DYdz[I,J,K] += DV[K,l] * XLoc[I,J,l,2]
      DZdz[I,J,K] += DV[K,l] * XLoc[I,J,l,3]
    end  
  end  
  @synchronize

  if IZ <= NZ && IF <= NF
    ID = I + (J - 1) * M
    t11 = DH[J,1] * (XLoc[I,1,K,2] * DZdz[I,1,K] - XLoc[I,1,K,3] * DYdz[I,1,K])
    t21 = DH[J,1] * (XLoc[I,1,K,3] * DXdz[I,1,K] - XLoc[I,1,K,1] * DZdz[I,1,K])
    t31 = DH[J,1] * (XLoc[I,1,K,1] * DYdz[I,1,K] - XLoc[I,1,K,2] * DXdz[I,1,K])
    t12 = DV[K,1] * (XLoc[I,J,1,2] * DZdy[I,J,1] - XLoc[I,J,1,3] * DYdy[I,J,1]) 
    t22 = DV[K,1] * (XLoc[I,J,1,3] * DXdy[I,J,1] - XLoc[I,J,1,1] * DZdy[I,J,1]) 
    t32 = DV[K,1] * (XLoc[I,J,1,1] * DYdy[I,J,1] - XLoc[I,J,1,2] * DXdy[I,J,1]) 

    @unroll for l = 2 : M
      t12 += DV[K,l] * (XLoc[I,J,l,2] * DZdy[I,J,l] - XLoc[I,J,l,3] * DYdy[I,J,l]) 
      t22 += DV[K,l] * (XLoc[I,J,l,3] * DXdy[I,J,l] - XLoc[I,J,l,1] * DZdy[I,J,l]) 
      t32 += DV[K,l] * (XLoc[I,J,l,1] * DYdy[I,J,l] - XLoc[I,J,l,2] * DXdy[I,J,l]) 
    end  
    @unroll for l = 2 : N
      t11 += DH[J,l] * (XLoc[I,l,K,2] * DZdz[I,l,K] - XLoc[I,l,K,3] * DYdz[I,l,K])
      t21 += DH[J,l] * (XLoc[I,l,K,3] * DXdz[I,l,K] - XLoc[I,l,K,1] * DZdz[I,l,K])
      t31 += DH[J,1] * (XLoc[I,l,K,1] * DYdz[I,l,K] - XLoc[I,l,K,2] * DXdz[I,l,K])
    end  

    dXdxI[1,1,K,ID,IZ,IF] = 0.5 * (t11 - t12)
    dXdxI[1,2,K,ID,IZ,IF] = 0.5 * (t21 - t22)
    dXdxI[1,3,K,ID,IZ,IF] = 0.5 * (t31 - t32)
  end
end

@kernel inbounds = true function MetricTriKernel1!(dXdxI,J,@Const(X),@Const(Dx2),
  @Const(DV), ::Val{DoF}, ::Val{M}) where {DoF,M}

  ID, K,   = @index(Local, NTuple)
  _,_,IZ,IF = @index(Global, NTuple)

  XLoc = @localmem eltype(dXdxI) (DoF,M,3)
  DXdy = @localmem eltype(dXdxI) (DoF,M)
  DXdz = @localmem eltype(dXdxI) (DoF,M)
  DYdy = @localmem eltype(dXdxI) (DoF,M)
  DYdz = @localmem eltype(dXdxI) (DoF,M)
  DZdy = @localmem eltype(dXdxI) (DoF,M)
  DZdz = @localmem eltype(dXdxI) (DoF,M)
  NZ = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]


  if IZ <= NZ && IF <= NF
    XLoc[ID,K,1] = X[ID,K,1,IZ,IF]
    XLoc[ID,K,2] = X[ID,K,2,IZ,IF]
    XLoc[ID,K,3] = X[ID,K,3,IZ,IF]
  end  
  @synchronize

  if IZ <= NZ && IF <= NF
    DXdy[ID,K] = Dx2[ID,1] * XLoc[1,K,1]
    DYdy[ID,K] = Dx2[ID,1] * XLoc[1,K,2]
    DZdy[ID,K] = Dx2[ID,1] * XLoc[1,K,3]
    @unroll for l = 2 : DoF
      DXdy[ID,K] += Dx2[ID,l] * XLoc[l,K,1]
      DYdy[ID,K] += Dx2[ID,l] * XLoc[l,K,2]
      DZdy[ID,K] += Dx2[ID,l] * XLoc[l,K,3]
    end
    DXdz[ID,K] = DV[K,1] * XLoc[ID,1,1]
    DYdz[ID,K] = DV[K,1] * XLoc[ID,1,2]
    DZdz[ID,K] = DV[K,1] * XLoc[ID,1,3]
    @unroll for l = 2 : M
      DXdz[ID,K] += DV[K,l] * XLoc[ID,l,1]
      DYdz[ID,K] += DV[K,l] * XLoc[ID,l,2]
      DZdz[ID,K] += DV[K,l] * XLoc[ID,l,3]
    end  
  end  

  @synchronize

  if IZ <= NZ && IF <= NF
    t11 = Dx2[ID,1] * (XLoc[1,K,2] * DZdz[1,K] - XLoc[1,K,3] * DYdz[1,K])
    t21 = Dx2[ID,1] * (XLoc[1,K,3] * DXdz[1,K] - XLoc[1,K,1] * DZdz[1,K])
    t31 = Dx2[ID,1] * (XLoc[1,K,1] * DYdz[1,K] - XLoc[1,K,2] * DXdz[1,K])
    @unroll for l = 2 : DoF
      t11 += Dx2[ID,l] * (XLoc[l,K,2] * DZdz[l,K] - XLoc[l,K,3] * DYdz[l,K])
      t21 += Dx2[ID,l] * (XLoc[l,K,3] * DXdz[l,K] - XLoc[l,K,1] * DZdz[l,K])
      t31 += Dx2[ID,l] * (XLoc[l,K,1] * DYdz[l,K] - XLoc[l,K,2] * DXdz[l,K])
    end  

    t12 = DV[K,1] * (XLoc[ID,1,2] * DZdy[ID,1] - XLoc[ID,1,3] * DYdy[ID,1]) 
    t22 = DV[K,1] * (XLoc[ID,1,3] * DXdy[ID,1] - XLoc[ID,1,1] * DZdy[ID,1]) 
    t32 = DV[K,1] * (XLoc[ID,1,1] * DYdy[ID,1] - XLoc[ID,1,2] * DXdy[ID,1]) 
    @unroll for l = 2 : M
      t12 += DV[K,l] * (XLoc[ID,l,2] * DZdy[ID,l] - XLoc[ID,l,3] * DYdy[ID,l]) 
      t22 += DV[K,l] * (XLoc[ID,l,3] * DXdy[ID,l] - XLoc[ID,l,1] * DZdy[ID,l]) 
      t32 += DV[K,l] * (XLoc[ID,l,1] * DYdy[ID,l] - XLoc[ID,l,2] * DXdy[ID,l]) 
    end  

    dXdxI[1,1,K,ID,IZ,IF] = 0.5 * (t11 - t12)
    dXdxI[1,2,K,ID,IZ,IF] = 0.5 * (t21 - t22)
    dXdxI[1,3,K,ID,IZ,IF] = 0.5 * (t31 - t32)
  end
end


@kernel inbounds = true function RotateKernel!(Rotate,@Const(X),::Grids.SphericalGrid)

  ID, K, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  if Iz <= Nz
    lon,lat,_ = Grids.cart2sphere(X[ID,K,1,Iz,IF],X[ID,K,2,Iz,IF],X[ID,K,3,Iz,IF])
    MR = Grids.MCart2Sphere(lon,lat)
    @. Rotate[:,:,K,ID,Iz,IF] =  MR
  end
end  

@kernel inbounds = true function RotateKernel!(Rotate,@Const(X),::Grids.CartesianGrid)

  ID, K, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  if Iz <= Nz
    Rotate[1,1,K,ID,Iz,IF] =  eltype(X)(1)
    Rotate[2,2,K,ID,Iz,IF] =  eltype(X)(1)
    Rotate[3,3,K,ID,Iz,IF] =  eltype(X)(1)
  end
end  

@kernel inbounds = true function FillXKernel!(X,Rad,AdaptGrid,@Const(ksi),@Const(zeta),@Const(F),@Const(z),@Const(zs),
  ElemType,GridForm)

  ID, K, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]
  H = @uniform z[Nz+1]

  if Iz <= Nz
    z1 = z[Iz]
    z2 = z[Iz+1]
    @views XPoint!(AdaptGrid,X[ID,K,:,Iz,IF],Rad,
    ksi[1,ID],ksi[2,ID],zeta[K],F[:,:,IF],z1,z2,H,zs[ID,IF],ElemType,GridForm)
  end
end  

@kernel inbounds = true function DetKernel!(J,@Const(dXdxI))

  ID, K, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  if Iz <= Nz && IF <= NF
    @views JLoc = Det3(dXdxI[:,:,K,ID,Iz,IF])
    J[ID,K,Iz,IF] = sqrt(abs(JLoc))
  end
end  

@inline function XPoint!(AdaptGrid,X,Rad,ksi1,ksi2,ksi3,F,z1,z2,H,zs,::Grids.Quad,::Grids.SphericalGrid)
  zero = eltype(X)(0)
  one = eltype(X)(1)
  half = eltype(X)(1/2)
  quarter = eltype(X)(1/4)
  X1 = quarter * (F[1,1] * (one-ksi1)*(one-ksi2) +
   F[2,1] * (one+ksi1)*(one-ksi2) +
   F[3,1] * (one+ksi1)*(one+ksi2) +
   F[4,1] * (one-ksi1)*(one+ksi2))
  X2 = quarter * (F[1,2] * (one-ksi1)*(one-ksi2) +
   F[2,2] * (one+ksi1)*(one-ksi2) +
   F[3,2] * (one+ksi1)*(one+ksi2) +
   F[4,2] * (one-ksi1)*(one+ksi2))
  X3 = quarter * (F[1,3] * (one-ksi1)*(one-ksi2) +
   F[2,3] * (one+ksi1)*(one-ksi2) +
   F[3,3] * (one+ksi1)*(one+ksi2) +
   F[4,3] * (one-ksi1)*(one+ksi2))
  zLoc = half * ((one-ksi3) * z1 + (one+ksi3) * z2)
  hR, = AdaptGrid(zLoc,zs)
  r = sqrt(X1 * X1 + X2 * X2 + X3 * X3)
  X[1] = X1 / r * (Rad + hR)
  X[2] = X2 / r * (Rad + hR)
  X[3] = X3 / r * (Rad + hR)
end

@inline function XPoint!(AdaptGrid,X,Rad,ksi1,ksi2,ksi3,F,z1,z2,H,zs,::Grids.Quad,::Grids.CartesianGrid)
  zero = eltype(X)(0)
  one = eltype(X)(1)
  half = eltype(X)(1/2)
  quarter = eltype(X)(1/4)
  X1 = quarter * (F[1,1] * (one-ksi1)*(one-ksi2) +
   F[2,1] * (one+ksi1)*(one-ksi2) +
   F[3,1] * (one+ksi1)*(one+ksi2) +
   F[4,1] * (one-ksi1)*(one+ksi2))
  X2 = quarter * (F[1,2] * (one-ksi1)*(one-ksi2) +
   F[2,2] * (one+ksi1)*(one-ksi2) +
   F[3,2] * (one+ksi1)*(one+ksi2) +
   F[4,2] * (one-ksi1)*(one+ksi2))
  X3 = quarter * (F[1,3] * (one-ksi1)*(one-ksi2) +
   F[2,3] * (one+ksi1)*(one-ksi2) +
   F[3,3] * (one+ksi1)*(one+ksi2) +
   F[4,3] * (one-ksi1)*(one+ksi2))
  zLoc = half * ((one-ksi3) * z1 + (one+ksi3) * z2)
  hR, = AdaptGrid(zLoc,zs)
  X[1] = X1
  X[2] = X2
  X[3] = X3 + hR
end

@inline function XPoint!(AdaptGrid,X,Rad,ksi1,ksi2,ksi3,F,z1,z2,H,zs,::Grids.Tri,::Grids.SphericalGrid)
  zero = eltype(X)(0)
  one = eltype(X)(1)
  half = eltype(X)(1/2)
  quarter = eltype(X)(1/4)
  X1 = half * (F[1,1] * (-ksi1 - ksi2) +
   F[2,1] * (eltype(X)(1)+ksi1) +
   F[3,1] * (eltype(X)(1)+ksi2))
  X2 = half * (F[1,2] * (-ksi1 - ksi2) +
   F[2,2] * (eltype(X)(1)+ksi1) +
   F[3,2] * (eltype(X)(1)+ksi2))
  X3 = half * (F[1,3] * (-ksi1 - ksi2) +
   F[2,3] * (eltype(X)(1)+ksi1) +
   F[3,3] * (eltype(X)(1)+ksi2))
  zLoc = half * ((one-ksi3) * z1 + (one+ksi3) * z2)
  hR, = AdaptGrid(zLoc,zs)
  r = sqrt(X1 * X1 + X2 * X2 + X3 * X3)
  X[1] = X1 / r * (Rad + hR)
  X[2] = X2 / r * (Rad + hR)
  X[3] = X3 / r * (Rad + hR)
end


function PointsFromGrid(backend,FT,Grid)
  if Grid.Type == Grids.Quad()
    nP = 4
  elseif Grid.Type == Grids.Tri()  
    nP = 3
  end  
  NF = Grid.NumFaces
  F = zeros(nP,3,NF)
  FGPU = KernelAbstractions.zeros(backend,FT,nP,3,NF)
  for iF = 1 : NF
    for iP = 1 : nP  
      F[iP,1,iF] = Grid.Faces[iF].P[iP].x
      F[iP,2,iF] = Grid.Faces[iF].P[iP].y
      F[iP,3,iF] = Grid.Faces[iF].P[iP].z
    end  
  end  
  copyto!(FGPU,F)
  return FGPU
end  

@inline function Det3(A)
  A[1,1] * (A[2,2] * A[3,3] - A[2,3] * A[3,2]) -
  A[1,2] * (A[2,1] * A[3,3] - A[2,3] * A[3,1]) +
  A[1,3] * (A[2,1] * A[3,2] - A[2,2] * A[3,1])
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

function GridSizeDGKernel!(DG,Metric,Rad,NumberThreadGPU,::Grids.CartesianGrid)
  backend = get_backend(Metric.dz)
  nz = size(Metric.dz,1) 
  DoF = size(Metric.X,1)
  NF = size(Metric.X,5)
  NzG = min(div(NumberThreadGPU,DoF),nz)
  group = (DoF,NzG,1)
  ndrange = (DoF,nz,NF)
  KGridSizeDGKernel! = GridSizeCartDGKernel!(backend,group)
  KGridSizeDGKernel!(Metric.dz,Metric.X,DG.Glob,ndrange=ndrange)
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

function GridSizeDGKernel!(DG,Metric,Rad,NumberThreadGPU,::Grids.SphericalGrid)
  backend = get_backend(Metric.dz)
  nz = size(Metric.dz,1) 
  DoF = size(Metric.X,1)
  NF = size(Metric.X,5)
  NzG = min(div(NumberThreadGPU,DoF),nz)
  group = (DoF,NzG,1)
  ndrange = (DoF,nz,NF)
  KGridSizeDGKernel! = GridSizeSphereDGKernel!(backend,group)
  KGridSizeDGKernel!(Metric.dz,Metric.X,DG.Glob,Rad,ndrange=ndrange)
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

@kernel inbounds = true function MetricLowerBoundaryDGKernel!(xS,@Const(X),@Const(Glob))

  ID,IF = @index(Global, NTuple)

  NF = @uniform @ndrange()[2]

  if IF <= NF
    ind = Glob[ID,IF]
    x = X[ID,1,1,1,IF]
    y = X[ID,1,2,1,IF]
    z = X[ID,1,3,1,IF]
    lon,lat,_ = Grids.cart2sphere(x,y,z)
    xS[1,ind] = lon
    xS[2,ind] = lat
  end
end
