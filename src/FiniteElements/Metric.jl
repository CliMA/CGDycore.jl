mutable struct MetricDGStruct{FT<:AbstractFloat,
                            AT2<:AbstractArray,
                            AT3<:AbstractArray,
                            AT4<:AbstractArray,
                            AT5<:AbstractArray,
                            AT6<:AbstractArray}
  J::AT4
  X::AT5
  dXdxI::AT6
  Rotate::AT6
  dz::AT2
  zP::AT2
  xS::AT2
  VolSurfH::AT4
  NH::AT5
  VolSurfV::AT3
  NV::AT4
end

function MetricCreate(backend,FT,nQuad,OPZ,NF,nz,NumG,::DGElement)
    J      = KernelAbstractions.zeros(backend,FT,nQuad,OPZ,nz,NF)
    X      = KernelAbstractions.zeros(backend,FT,nQuad,OPZ,3,nz,NF)
    dXdxI  = KernelAbstractions.zeros(backend,FT,3,3,OPZ,nQuad,nz,NF)
    Rotate  = KernelAbstractions.zeros(backend,FT,3,3,OPZ,nQuad,nz,NF)
    dz = KernelAbstractions.zeros(backend,FT,0,0)
    zP = KernelAbstractions.zeros(backend,FT,0,0)
    xS    = KernelAbstractions.zeros(backend,FT,2,NumG)
    VolSurfH = KernelAbstractions.zeros(backend,FT,0,0,0,0)
    NH = KernelAbstractions.zeros(backend,FT,0,0,0,0,0)
    VolSurfV = KernelAbstractions.zeros(backend,FT,0,0,0)
    NV = KernelAbstractions.zeros(backend,FT,0,0,0,0)
    return MetricDGStruct{FT,
                        typeof(zP),
                        typeof(VolSurfV),
                        typeof(J),
                        typeof(X),
                        typeof(dXdxI)}(
        J,
        X,
        dXdxI,
        Rotate,
        dz,
        zP,
        xS,
        VolSurfH,
        NH,
        VolSurfV,
        NV,
    )
end

mutable struct MetricCGStruct{FT<:AbstractFloat,
                            AT2<:AbstractArray,
                            AT3<:AbstractArray,
                            AT4<:AbstractArray,
                            AT5<:AbstractArray,
                            AT6<:AbstractArray}
  J::AT4
  X::AT5
  dXdxI::AT6
  dXdx::AT6
  nSS::AT2
  nS::AT3
  FS::AT2
  dz::AT2
  zP::AT2
  JC::AT3
  JCW::AT3
  xS::AT2
  M::AT3
  MMass::AT2
end

function MetricCreate(backend,FT,nQuad,OPZ,NF,nz,NumG,::CGElement)
    J      = KernelAbstractions.zeros(backend,FT,nQuad,OPZ,nz,NF)
    X      = KernelAbstractions.zeros(backend,FT,nQuad,OPZ,3,nz,NF)
    dXdxI  = KernelAbstractions.zeros(backend,FT,3,3,OPZ,nQuad,nz,NF)
    dXdx   = KernelAbstractions.zeros(backend,FT,3,3,OPZ,nQuad,nz,NF)
    nSS  = KernelAbstractions.zeros(backend,FT,3,NumG)
    nS = KernelAbstractions.zeros(backend,FT,nQuad,3,NF)
    FS = KernelAbstractions.zeros(backend,FT,nQuad,NF)
    dz = KernelAbstractions.zeros(backend,FT,0,0)
    zP = KernelAbstractions.zeros(backend,FT,0,0)
    JC     = KernelAbstractions.zeros(backend,FT,0,0,0)
    JCW    = KernelAbstractions.zeros(backend,FT,0,0,0)
    xS    = KernelAbstractions.zeros(backend,FT,2,NumG)
    VolSurfH = KernelAbstractions.zeros(backend,FT,0,0,0,0)
    NH = KernelAbstractions.zeros(backend,FT,0,0,0,0,0)
    VolSurfV = KernelAbstractions.zeros(backend,FT,0,0,0)
    NV = KernelAbstractions.zeros(backend,FT,0,0,0,0)
    M = KernelAbstractions.zeros(backend,FT,0,0,0)
    MMass = KernelAbstractions.zeros(backend,FT,0,0)
    return MetricCGStruct{FT,
                        typeof(zP),
                        typeof(nS),
                        typeof(J),
                        typeof(X),
                        typeof(dXdxI)}(
        J,
        X,
        dXdxI,
        dXdx,
        nSS,
        nS, 
        FS, 
        dz,
        zP,
        JC,
        JCW,
        xS,
        M,
        MMass,
    )
end

function MetricCompute(backend,FT,FE::DGElement,Model,Exchange,Grid,NumberThreadGPU,zS)
  DoF = FE.DoF
  DoFE = FE.DoFE
  M = FE.OrdPolyZ + 1
  NF = Grid.NumFaces
  NE = Grid.NumEdges
  Nz = Grid.nz

  EdgeFace!(backend,Grid)
  Metric = MetricCreate(backend,FT,DoF,M,NF,Nz,FE.NumG,FE)
  Metric.zP = KernelAbstractions.zeros(backend,FT,Nz,FE.NumG)
  Metric.dz = KernelAbstractions.zeros(backend,FT,Nz,FE.NumG)
  F = PointsFromGrid(backend,FT,Grid)
  FillX!(backend,Metric,FE,F,Grid,zS)
  FillRotate!(backend,Metric,FE,Grid)
  FillContravariant!(backend,Metric,FE,Grid,Grid.Type,Model.MetricType)
  FillDet!(backend,Metric,FE,Grid)
  GridSizeDGKernel!(FE,Metric,Grid.Rad,NumberThreadGPU,Grid.Form)
  MetricLowerBoundary!(backend,Metric,FE,Grid,NumberThreadGPU,Grid.Form)
  NormalH!(backend,Metric,FE,Grid,NumberThreadGPU,Grid.Type)
  NormalV!(backend,Metric,FE,Grid,NumberThreadGPU)  

  return Metric
end

function MetricCompute(backend,FT,FE::CGElement,Model,Exchange,Grid,NumberThreadGPU,zS)
  DoF = FE.DoF
  DoFE = FE.DoFE
  M = FE.OrdPolyZ + 1
  NumG = FE.NumG
  NF = Grid.NumFaces
  NE = Grid.NumEdges
  Nz = Grid.nz

  EdgeFace!(backend,Grid)
  Metric = MetricCreate(backend,FT,DoF,M,NF,Nz,FE.NumG,FE)
  Metric.zP = KernelAbstractions.zeros(backend,FT,Nz,FE.NumG)
  Metric.dz = KernelAbstractions.zeros(backend,FT,Nz,FE.NumG)
  F = PointsFromGrid(backend,FT,Grid)

  if occursin("DGMetric",Model.MetricType)
    FillX!(backend,Metric,FE,F,Grid,zS)
    FillContravariant!(backend,Metric,FE,Grid,Grid.Type,Metric.Type)
    FillDet!(backend,Metric,FE,Grid)
  else  
    Grids.JacobiDG3GPU!(Grid.AdaptGrid,Metric.X,Metric.dXdxI,Metric.J,FE,F,Grid.z,zS,
      Grid.Rad,Grid.Type,Grid.Form)
    Grids.JacobiSphere3GPU!(Grid.AdaptGrid,Metric.X,Metric.dXdxI,Metric.J,FE,F,Grid.z,zS,
    Grid.Rad,Models.CompressibleShallow())
  end  

  MassMatrix!(backend,FT,Metric,FE,Exchange,NumberThreadGPU)

  NzG = min(div(NumberThreadGPU,DoF),Nz)
  group = (DoF,NzG,1)
  ndrange = (DoF,Nz,NF)
  if Grid.Form == Grids.SphericalGrid()
    KGridSizeSphereKernel! = GridSizeSphereKernel!(backend,group)
    Rad = Grid.Rad
    KGridSizeSphereKernel!(Metric.zP,Metric.dz,Metric.X,FE.Glob,
      Rad,ndrange=ndrange)
  else
    KGridSizeCartKernel! = GridSizeCartKernel!(backend,group)
    KGridSizeCartKernel!(Metric.zP,Metric.dz,Metric.X,FE.Glob,ndrange=ndrange)
  end
  NFG = min(div(NumberThreadGPU,DoF),NF)
  group = (DoF, NFG)
  ndrange = (DoF, NF)
  KSurfaceNormalKernel! = SurfaceNormalKernel!(backend,group)
  KSurfaceNormalKernel!(Metric.FS,Metric.nS,Metric.dXdxI,ndrange=ndrange)
  KMetricLowerBoundaryKernel! = MetricLowerBoundaryKernel!(backend,group)
  KMetricLowerBoundaryKernel!(Metric.nSS,Metric.xS,Metric.dXdxI,Metric.X,Metric.M,FE.Glob,ndrange=ndrange)
  Parallels.ExchangeData!(Metric.nSS,Exchange)
  groupS = (NumberThreadGPU)
  ndrangeS = (NumG)
  KMetricLowerBoundaryScaleKernel! = MetricLowerBoundaryScaleKernel!(backend,groupS)
  KMetricLowerBoundaryScaleKernel!(Metric.nSS,ndrange=ndrangeS)

  HorLimit = false
  if HorLimit
    Metric.JC = KernelAbstractions.zeros(backend,FT,size(Metric.J,1),size(Metric.J,3),size(Metric.J,4))
    Metric.JCW = KernelAbstractions.zeros(backend,FT,size(Metric.J,1),size(Metric.J,3),size(Metric.J,4))
#   NFG = min(div(NumberThreadGPU,),NF)
    group = (Nz, 1)
    ndrange = (Nz, NF)
    KCenterJacobiansKernel! = CenterJacobiansKernel!(backend,group)
    KCenterJacobiansKernel!(FE.OrdPoly+1,Metric.JC,Metric.JCW,Metric.J,FE.w,ndrange=ndrange)
  end

  return Metric
end

function FillContravariant!(backend,Metric,FE::DGElement,Grid,::Grids.Quad,MetricType)
  FT = eltype(Metric.X)
  DoF = FE.DoF
  DoFE = FE.DoFE
  N = FE.OrdPoly + 1
  M = FE.OrdPolyZ + 1
  NF = Grid.NumFaces
  Nz = Grid.nz
  group = (N,N,M,1,1)
  ndrange = (N,N,M,Nz,NF)
  CurlMetric = false
  _,DS,_,_,_ = DG.DerivativeMatrixSingle(FE.OrdPoly)
  DSGPU = KernelAbstractions.zeros(backend,FT,size(DS))
  copyto!(DSGPU,DS)
  _,DSZ,_,_,_ = DG.DerivativeMatrixSingle(FE.OrdPolyZ)
  DSZGPU = KernelAbstractions.zeros(backend,FT,size(DSZ))
  copyto!(DSZGPU,DSZ)
  if occursin("Curl",MetricType) 
    KFillContraKernel1! = FillContraCurlQuadKernel1!(backend,group)
    KFillContraKernel1!(Metric.dXdxI,Metric.X,FE.DS,FE.DSZ,Val(N),Val(M);ndrange=ndrange)
    KFillContraKernel2! = FillContraCurlQuadKernel2!(backend,group)
    KFillContraKernel2!(Metric.dXdxI,Metric.X,FE.DS,FE.DSZ,Val(N),Val(M);ndrange=ndrange)
    KFillContraKernel3! = FillContraCurlQuadKernel3!(backend,group)
    KFillContraKernel3!(Metric.dXdxI,Metric.X,FE.DS,Val(N),Val(M);ndrange=ndrange)
  else
    KFillContraKernel! = MetricQuadKernel!(backend,group)
    KFillContraKernel!(Metric.dXdxI,Metric.X,DS,DSZ,Val(N),Val(M);ndrange=ndrange)
  end
end

function FillContravariant!(backend,Metric,FE::CGElement,Grid,::Grids.Quad,MetricType)
  FT = eltype(Metric.X)
  DoF = FE.DoF
  DoFE = FE.DoFE
  N = FE.OrdPoly + 1
  M = FE.OrdPolyZ + 1
  NF = Grid.NumFaces
  Nz = Grid.nz
  group = (N,N,M,1,1)
  ndrange = (N,N,M,Nz,NF)
  _,DS,_,_,_ = DG.DerivativeMatrixSingle(FE.OrdPoly)
  DSGPU = KernelAbstractions.zeros(backend,FT,size(DS))
  copyto!(DSGPU,DS)
  _,DSZ,_,_,_ = DG.DerivativeMatrixSingle(FE.OrdPolyZ)
  DSZGPU = KernelAbstractions.zeros(backend,FT,size(DSZ))
  copyto!(DSZGPU,DSZ)
  if occursin("Curl",MetricType) 
    KFillContraKernel1! = FillContraCurlQuadKernel1!(backend,group)
    KFillContraKernel1!(Metric.dXdxI,Metric.X,FE.DS,FE.DSZ,Val(N),Val(M);ndrange=ndrange)
    KFillContraKernel2! = FillContraCurlQuadKernel2!(backend,group)
    KFillContraKernel2!(Metric.dXdxI,Metric.X,FE.DS,FE.DSZ,Val(N),Val(M);ndrange=ndrange)
    KFillContraKernel3! = FillContraCurlQuadKernel3!(backend,group)
    KFillContraKernel3!(Metric.dXdxI,Metric.X,FE.DS,Val(N),Val(M);ndrange=ndrange)
  else
    KFillContraKernel! = MetricQuadKernel!(backend,group)
    KFillContraKernel!(Metric.dXdxI,Metric.X,DS,DSZ,Val(N),Val(M);ndrange=ndrange)
  end  
  group = (DoF,M,1,1)
  ndrange = (DoF,M,Nz,NF)
  KApplyRotateKernel! = ApplyRotateKernel!(backend,group)
  KApplyRotateKernel!(Metric.dXdxI,Metric.X,Grid.Form;ndrange=ndrange)
end

function FillContravariant!(backend,Metric,FE,Grid,::Grids.Tri)
  DoF = FE.DoF
  M = FE.OrdPolyZ + 1
  NF = Grid.NumFaces
  Nz = Grid.nz
  group = (DoF,M,1,1)
  ndrange = (DoF,M,Nz,NF)
  CurlMetric = true
  if CurlMetric
    KFillContraKernel1! = MetricTriKernel1!(backend,group)
    KFillContraKernel1!(Metric.dXdxI,Metric.J,Metric.X,FE.Dx2,FE.DSZ,Val(DoF),Val(M);ndrange=ndrange)
    KFillContraKernel2! = MetricTriKernel2!(backend,group)
    KFillContraKernel2!(Metric.dXdxI,Metric.J,Metric.X,FE.Dx1,FE.DSZ,Val(DoF),Val(M);ndrange=ndrange)
    KFillContraKernel3! = MetricTriKernel3!(backend,group)
    KFillContraKernel3!(Metric.dXdxI,Metric.J,Metric.X,FE.Dx1,FE.Dx2,Val(DoF),Val(M);ndrange=ndrange)
  end
end

@kernel inbounds = true function FillContraCurlQuadKernel3!(dXdxI,@Const(X),@Const(DH), ::Val{N}, ::Val{M}) where {N,M}

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
    ID = I + (J - 1) * N
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

    ID = I + (J - 1) * N

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

@kernel inbounds = true function FillContraCurlQuadKernel2!(dXdxI,@Const(X),@Const(DH),@Const(DV), ::Val{N}, ::Val{M}) where {N,M}

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
    ID = I + (J - 1) * N
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
    ID = I + (J - 1) * N
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

@kernel inbounds = true function FillContraCurlQuadKernel1!(dXdxI,@Const(X),@Const(DH),@Const(DV), ::Val{N}, ::Val{M}) where {N,M}

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
    ID = I + (J - 1) * N
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
    ID = I + (J - 1) * N
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

@kernel inbounds = true function ApplyRotateKernel!(dXdxI,@Const(X),::Grids.SphericalGrid)

  ID, K, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  if Iz <= Nz && IF <= NF
    lon,lat,_ = Grids.cart2sphere(X[ID,K,1,Iz,IF],X[ID,K,2,Iz,IF],X[ID,K,3,Iz,IF])
    MR = Grids.MCart2Sphere(lon,lat)
    dXdxI[:,:,K,ID,Iz,IF] = dXdxI[:,:,K,ID,Iz,IF] * MR'
  end
end  

function FillRotate!(backend,Metric,FE,Grid)
  DoF = FE.DoF
  M = FE.OrdPolyZ + 1
  NF = Grid.NumFaces
  Nz = Grid.nz
  group = (DoF,M,1,1)
  ndrange = (DoF,M,Nz,NF)
  KRotateKernel! = RotateKernel!(backend,group)
  KRotateKernel!(Metric.Rotate,Metric.X,Grid.Form;ndrange=ndrange)
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

function FillX!(backend,Metric,FE,F,Grid,zS)
  DoF = FE.DoF
  M = FE.OrdPolyZ + 1
  NF = Grid.NumFaces
  Nz = Grid.nz
  group = (DoF,M,1,1)
  ndrange = (DoF,M,Nz,NF)
  KFillXKernel! = FillXKernel!(backend,group)
  KFillXKernel!(Metric.X,Grid.Rad,Grid.AdaptGrid,FE.ksi,FE.xwZ,F,Grid.z,zS,Grid.Type,Grid.Form;ndrange=ndrange)
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

function FillDet!(backend,Metric,FE,Grid)
  DoF = FE.DoF
  M = FE.OrdPolyZ + 1
  NF = Grid.NumFaces
  Nz = Grid.nz
  group = (DoF,M,1,1)
  ndrange = (DoF,M,Nz,NF)
  KDetKernel! = DetKernel!(backend,group)
  KDetKernel!(Metric.J,Metric.dXdxI;ndrange=ndrange)
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

function NormalH!(backend,Metric,FE::DGElement,Grid,NumberThreadGPU,::Grids.Quad)
  N = FE.OrdPoly + 1
  M = FE.OrdPolyZ + 1
  NE = Grid.NumEdges
  Nz = Grid.nz
  NzG = min(div(NumberThreadGPU,N*M),Nz)
  group = (M,N,NzG)
  ndrange = (M,N,Nz,NE)
  FT = eltype(Metric.dXdxI)
  KNormalHQuadKernel! = NormalHQuadKernel!(backend,group)
  Metric.VolSurfH = KernelAbstractions.zeros(backend,FT,M,N,Nz,NE)
  Metric.NH = KernelAbstractions.zeros(backend,FT,3,M,N,Nz,NE)
  KNormalHQuadKernel!(Metric.VolSurfH,Metric.NH,
    Metric.dXdxI,Grid.EF,Grid.FE,ndrange=ndrange)
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
      if IF1 > 0  
        IF = IF1  
        IS = FE[1,IE]  
      else  
        IF = IF2
        IS = FE[2,IE]  
      end  
    else
      if IF2 > 0  
        IF = IF2
        IS = FE[2,IE]  
      else
        IF = IF1  
        IS = FE[1,IE]  
      end  
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

function NormalH!(backend,Metric,FE::DGElement,Grid,NumberThreadGPU,::Grids.Tri)
  FT = eltype(Metric.dXdxI)
  N = FE.DoFE
  M = FE.OrdPolyZ + 1
  NE = Grid.NumEdges
  Nz = Grid.nz
  NzG = min(div(NumberThreadGPU,N*M),Nz)
  group = (M,N,NzG)
  ndrange = (M,N,Nz,NE)
  KNormalHTriKernel! = NormalHTriKernel!(backend,group)
  Metric.VolSurfH = KernelAbstractions.zeros(backend,FT,M,N,Nz,NE)
  Metric.NH = KernelAbstractions.zeros(backend,FT,3,M,N,Nz,NE)
  KNormalHTriKernel!(Metric.VolSurfH,Metric.NH,
    Metric.dXdxI,Grid.EF,Grid.FE,FE.PosDoFE,ndrange=ndrange)
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

function NormalV!(backend,Metric,FE::DGElement,Grid,NumberThreadGPU)  
  DoF = FE.DoF
  M = FE.OrdPolyZ + 1
  NF = Grid.NumFaces
  Nz = Grid.nz
  FT = eltype(Metric.dXdxI)
  NzG = min(div(NumberThreadGPU,DoF),Nz+1)
  group = (DoF,NzG,1)
  ndrange = (DoF,Nz+1,NF)
  KNormalVKernel! = NormalVKernel!(backend,group)
  Metric.VolSurfV = KernelAbstractions.zeros(backend,FT,DoF,Nz+1,NF)
  Metric.NV = KernelAbstractions.zeros(backend,FT,3,DoF,Nz+1,NF)
  KNormalVKernel!(Metric.VolSurfV,Metric.NV,M,Metric.dXdxI,ndrange=ndrange)
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

function MetricLowerBoundary!(backend,Metric,FE::DGElement,Grid,NumberThreadGPU,::Grids.CartesianGrid)
end  

function MetricLowerBoundary!(backend,Metric,FE::DGElement,Grid,NumberThreadGPU,::Grids.SphericalGrid)
  DoF = FE.DoF
  NF = Grid.NumFaces
  NFG = min(div(NumberThreadGPU,DoF),NF)
  group = (DoF, NFG)
  ndrange = (DoF, NF)
  KMetricLowerBoundaryKernel! = MetricLowerBoundaryDGKernel!(backend,group)
  KMetricLowerBoundaryKernel!(Metric.xS,Metric.X,FE.Glob,ndrange=ndrange)
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

function EdgeFace!(backend,Grid)
  NE = Grid.NumEdges
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
end  

function MassMatrix!(backend,FT,Metric,FE::CGElement,Exchange,NumberThreadGPU)
  N = FE.OrdPoly + 1
  DoF = FE.DoF
  NumG = FE.NumG
  w = FE.w
  J = Metric.J
  Glob = FE.Glob
  Nz = size(J,3)
  NF = size(Glob,2)
  Metric.M = KernelAbstractions.zeros(backend,FT,Nz,NumG,2)
  Metric.MMass = KernelAbstractions.zeros(backend,FT,Nz,NumG)
  M = Metric.M
  MMass = Metric.MMass

  M .= FT(0)
  MMass .= FT(0)


  NzG = min(div(NumberThreadGPU,N*N),Nz)
  group = (N, N, NzG, 1)
  ndrange = (N, N, Nz, NF)

  KMassCGKernel! = MassCGKernel!(backend,group)
  KMassCGKernel!(M,MMass,J,w,Glob,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)

  @views Parallels.ExchangeData!(M[:,:,1],Exchange)
  @views Parallels.ExchangeData!(M[:,:,2],Exchange)
  Parallels.ExchangeData!(MMass,Exchange)
end

@kernel inbounds = true function MassCGKernel!(M,MMass,@Const(JJ),@Const(w),@Const(Glob))
  I,J,Iz,IF = @index(Global, NTuple)

  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  if Iz <= Nz && IF <= NF
    ID = I + (J - 1) * N  
    ind = Glob[ID,IF]  
    @atomic :monotonic M[Iz,ind,1] += JJ[ID,1,Iz,IF]
    @atomic :monotonic M[Iz,ind,2] += JJ[ID,2,Iz,IF]
    @atomic :monotonic MMass[Iz,ind] += eltype(M)(0.5) * (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]) * w[I] * w[J]
  end  
end

@kernel inbounds = true function GridSizeSphereKernel!(zP,dz,@Const(X),@Const(Glob),Rad)

  ID,Iz,IF = @index(Global, NTuple)

  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]


  if Iz <= Nz && IF <= NF
    ind = Glob[ID,IF]
    x = eltype(X)(0.5) * (X[ID,1,1,Iz,IF] + X[ID,2,1,Iz,IF])
    y = eltype(X)(0.5) * (X[ID,1,2,Iz,IF] + X[ID,2,2,Iz,IF])
    z = eltype(X)(0.5) * (X[ID,1,3,Iz,IF] + X[ID,2,3,Iz,IF])
    r = sqrt(x * x + y * y + z * z)
    zP[Iz,ind] = max(r-Rad, eltype(X)(0))
    x = X[ID,1,1,Iz,IF]
    y = X[ID,1,2,Iz,IF]
    z = X[ID,1,3,Iz,IF]
    r1 = sqrt(x * x + y * y + z * z)
    x = X[ID,2,1,Iz,IF]
    y = X[ID,2,2,Iz,IF]
    z = X[ID,2,3,Iz,IF]
    r2 = sqrt(x * x + y * y + z * z)
    dz[Iz,ind] =  r2 - r1
  end
end

@kernel inbounds = true function SurfaceNormalKernel!(FS,nS,@Const(dXdxI))

  ID,IF = @index(Global, NTuple)

  NF = @uniform @ndrange()[2]

  if IF <= NF
    FS[ID,IF] = sqrt(dXdxI[3,1,1,ID,1,IF] * dXdxI[3,1,1,ID,1,IF] +
      dXdxI[3,2,1,ID,1,IF] * dXdxI[3,2,1,ID,1,IF] +
      dXdxI[3,3,1,ID,1,IF] * dXdxI[3,3,1,ID,1,IF])
    nS[ID,1,IF] = dXdxI[3,1,1,ID,1,IF] / FS[ID,IF]
    nS[ID,2,IF] = dXdxI[3,2,1,ID,1,IF] / FS[ID,IF]
    nS[ID,3,IF] = dXdxI[3,3,1,ID,1,IF] / FS[ID,IF]
  end
end

@kernel inbounds = true function CenterJacobiansKernel!(N,JC,JCW,@Const(J),@Const(w))
  Iz,IF = @index(Global, NTuple)

  FT = eltype(JC)
  Nz = @uniform @ndrange()[1]
  if Iz <= Nz
    @views sumJ = sum(J[:,:,Iz,IF])
    for j = 1 : N
      for i = 1 : N
        ID = i + (j - 1) * N 
        JC[ID,Iz,IF] = J[ID,1,Iz,IF] + J[ID,2,Iz,IF]
        JCW[ID,Iz,IF] = JC[ID,Iz,IF] * w[i] * w[j] / sumJ
      end
    end 
  end 
end

@kernel inbounds = true function MetricLowerBoundaryKernel!(nS,xS,@Const(dXdxI),@Const(X),@Const(M),@Const(Glob))

  ID,IF = @index(Global, NTuple)

  NF = @uniform @ndrange()[2]

  if IF <= NF
    ind = Glob[ID,IF]
    @atomic :monotonic nS[1,ind] += dXdxI[3,1,1,ID,1,IF] / (M[1,ind,1] + M[1,ind,2])
    @atomic :monotonic nS[2,ind] += dXdxI[3,2,1,ID,1,IF] / (M[1,ind,1] + M[1,ind,2])
    @atomic :monotonic nS[3,ind] += dXdxI[3,3,1,ID,1,IF] / (M[1,ind,1] + M[1,ind,2])
    x = X[ID,1,1,1,IF]
    y = X[ID,1,2,1,IF]
    z = X[ID,1,3,1,IF]
    lon,lat,_ = Grids.cart2sphere(x,y,z)
    xS[1,ind] = lon
    xS[2,ind] = lat
  end
end

@kernel inbounds = true function MetricLowerBoundaryScaleKernel!(nS)

  IC, = @index(Global, NTuple)

  NumG = @uniform @ndrange()[1]

  if IC <= NumG
    fac = eltype(nS)(1) / sqrt(nS[1,IC]^2 + nS[2,IC]^2 + nS[3,IC]^2)
    nS[1,IC] *= fac
    nS[2,IC] *= fac
    nS[3,IC] *= fac
  end
end
