mutable struct CGQuad{FT<:AbstractFloat,
                        AT1<:AbstractArray,
                        AT2<:AbstractArray,
                        AT3<:AbstractArray,
                        IT1<:AbstractArray,
                        IT2<:AbstractArray}
    OrdPoly::Int
    OrdPolyZ::Int
    DoF::Int
    Glob::IT2
    Stencil::IT2
    NumG::Int
    NumI::Int
    w::AT1
    xw::AT1
    xwCPU::Array{FT, 1}
    xe::Array{FT, 1}
    IntXE2F::Array{FT, 2}
    xwZ::AT1
    xwZCPU::Array{FT, 1}
    IntZE2F::Array{FT, 2}
    DW::AT2
    DWT::Array{FT, 2}
    DS::AT2
    DST::Array{FT, 2}
    DSZ::AT2
    S::Array{FT, 2}
    M::AT3
    MMass::AT2
    BoundaryDoF::Array{Int, 1}
    MasterSlave::IT1
end

function CGQuad{FT}(backend,OrdPoly,OrdPolyZ,Grid) where FT<:AbstractFloat
# Discretization
  OP=OrdPoly+1
  OPZ=OrdPolyZ+1
  nz = Grid.nz

# CG = CGStruct{FT}(backend)
  OrdPoly=OrdPoly
  OrdPolyZ=OrdPolyZ
  DoF = OP * OP

  xwCPU, wCPU = gausslobatto(OrdPoly+1)
  w = KernelAbstractions.zeros(backend,FT,size(wCPU))
  xw = KernelAbstractions.zeros(backend,FT,size(xwCPU))
  copyto!(w,wCPU)
  copyto!(xw,xwCPU)
  
  xwZCPU, wZ = gausslobatto(OrdPolyZ+1)
  xwZ = KernelAbstractions.zeros(backend,FT,size(xwZCPU))
  copyto!(xwZ,xwZCPU)
  xe = zeros(OrdPoly+1)
  xe[1] = -1.0
  for i = 2 : OrdPoly
    xe[i] = xe[i-1] + 2.0/OrdPoly
  end
  xe[OrdPoly+1] = 1.0

  IntXE2F = zeros(OrdPoly+1,OrdPoly+1)
  for j = 1 : OrdPoly + 1
    for i = 1 : OrdPoly +1
      IntXE2F[i,j] = DG.Lagrange(xwCPU[i],xe,j)
    end
  end

  IntZE2F = zeros(OrdPolyZ+1,OrdPolyZ+1)
  for j = 1 : OrdPolyZ + 1
    for i = 1 : OrdPolyZ +1
      IntZE2F[i,j] = DG.Lagrange(xwZCPU[i],xwZCPU,j)
    end
  end

  (DWCPU,DSCPU)=DG.DerivativeMatrixSingle(OrdPoly)
  DS = KernelAbstractions.zeros(backend,FT,size(DSCPU))
  copyto!(DS,DSCPU)
  DW = KernelAbstractions.zeros(backend,FT,size(DWCPU))
  copyto!(DW,DWCPU)
  DST=DS'
  DWT=DW'

  Q = diagm(wCPU) * DSCPU
  S = Q - Q'
  (DWZ,DSZCPU)=DG.DerivativeMatrixSingle(OrdPolyZ)
  DSZ = KernelAbstractions.zeros(backend,FT,size(DSZCPU))
  copyto!(DSZ,DSZCPU)
  (GlobCPU,NumG,NumI,StencilCPU,MasterSlaveCPU,BoundaryDoFCPU) =
    NumberingFemCGQuad(Grid,OrdPoly)  

  Glob = KernelAbstractions.zeros(backend,Int,size(GlobCPU))
  Stencil = KernelAbstractions.zeros(backend,Int,size(StencilCPU))
  copyto!(Stencil,StencilCPU)
  MasterSlave = KernelAbstractions.zeros(backend,Int,size(MasterSlaveCPU))
  copyto!(MasterSlave,MasterSlaveCPU)
  BoundaryDoF = KernelAbstractions.zeros(backend,Int,size(BoundaryDoFCPU))
  copyto!(BoundaryDoF,BoundaryDoFCPU)
  copyto!(Glob,GlobCPU)
  M = KernelAbstractions.zeros(backend,FT,nz,NumG,2)
  MMass = KernelAbstractions.zeros(backend,FT,nz,NumG)
  return CGQuad{FT,
                 typeof(w),
                 typeof(DW),
                 typeof(M),
                 typeof(MasterSlave),
                 typeof(Stencil)}(
    OrdPoly,
    OrdPolyZ,
    DoF,
    Glob,
    Stencil,
    NumG,
    NumI,
    w,
    xw,
    xwCPU,
    xe,
    IntXE2F,
    xwZ,
    xwZCPU,
    IntZE2F,
    DW,
    DWT,
    DS,
    DST,
    DSZ,
    S,
    M,
    MMass,
    BoundaryDoF,
    MasterSlave,
 )
end

mutable struct DGQuad{FT<:AbstractFloat,
                        AT1<:AbstractArray,
                        AT2<:AbstractArray,
                        IT1<:AbstractArray,
                        IT2<:AbstractArray,
                        IT3<:AbstractArray} 
    OrdPoly::Int
    OrdPolyZ::Int
    DoF::Int
    Glob::IT2
    GlobE::IT3
    Stencil::IT2
    NumG::Int
    NumI::Int
    w::AT1
    xw::AT1
    xwCPU::Array{FT, 1}
    xe::Array{FT, 1}
    IntXE2F::Array{FT, 2}
    wZ::AT1
    xwZ::AT1
    xwZCPU::Array{FT, 1}
    IntZE2F::Array{FT, 2}
    DW::AT2
    DWT::Array{FT, 2}
    DS::AT2
    DST::Array{FT, 2}
    DSZ::AT2
    DWZ::AT2
    DV::AT2
    DVT::AT2
    DVZ::AT2
    DVZT::AT2
    S::Array{FT, 2}
    BoundaryDoF::Array{Int, 1}
    MasterSlave::IT1
end

function DGQuad{FT}(backend,OrdPoly,OrdPolyZ,Grid,Proc) where FT<:AbstractFloat
# Discretization
  OP=OrdPoly+1
  OPZ=OrdPolyZ+1
  nz = Grid.nz

  DoF = OP * OP

  xwCPU, wCPU = gausslobatto(OrdPoly+1)
  w = KernelAbstractions.zeros(backend,FT,size(wCPU))
  xw = KernelAbstractions.zeros(backend,FT,size(xwCPU))
  copyto!(w,wCPU)
  copyto!(xw,xwCPU)
  
  if OrdPolyZ == 0
    xwZCPU = zeros(1)  
    wZCPU = 2 * ones(1)
  else    
    xwZCPU, wZCPU = gausslobatto(OrdPolyZ+1)
  end  
  xwZ = KernelAbstractions.zeros(backend,FT,size(xwZCPU))
  wZ = KernelAbstractions.zeros(backend,FT,size(wZCPU))
  copyto!(xwZ,xwZCPU)
  copyto!(wZ,wZCPU)
  xe = zeros(OrdPoly+1)
  xe[1] = -1.0
  for i = 2 : OrdPoly
    xe[i] = xe[i-1] + 2.0/OrdPoly
  end
  xe[OrdPoly+1] = 1.0

  IntXE2F = zeros(OrdPoly+1,OrdPoly+1)
  for j = 1 : OrdPoly + 1
    for i = 1 : OrdPoly +1
      IntXE2F[i,j] = DG.Lagrange(xwCPU[i],xe,j)
    end
  end

  IntZE2F = zeros(OrdPolyZ+1,OrdPolyZ+1)
  for j = 1 : OrdPolyZ + 1
    for i = 1 : OrdPolyZ +1
      IntZE2F[i,j] = DG.Lagrange(xwZCPU[i],xwZCPU,j)
    end
  end

  (DWCPU,DSCPU)=DG.DerivativeMatrixSingle(OrdPoly)
  DS = KernelAbstractions.zeros(backend,FT,size(DSCPU))
  copyto!(DS,DSCPU)
  DW = KernelAbstractions.zeros(backend,FT,size(DWCPU))
  copyto!(DW,DWCPU)
  DST=DS'
  DWT=DW'

  DV = KernelAbstractions.zeros(backend,FT,size(DSCPU))
  DVCPU = 2 * DSCPU
  DVCPU[1,1] += 1 / wCPU[1]
  DVCPU[OrdPoly+1,OrdPoly+1] += -1 / wCPU[OrdPoly+1]
  copyto!(DV,DVCPU)
  DVT=DV'

  Q = diagm(wCPU) * DSCPU
  S = Q - Q'
  (DWZCPU,DSZCPU)=DG.DerivativeMatrixSingle(OrdPolyZ)
  DSZ = KernelAbstractions.zeros(backend,FT,size(DSZCPU))
  copyto!(DSZ,DSZCPU)
  DWZ = KernelAbstractions.zeros(backend,FT,size(DWZCPU))
  copyto!(DWZ,DWZCPU)
  (GlobCPU,GlobECPU,NumG,NumI,StencilCPU,MasterSlaveCPU,BoundaryDoFCPU) =
    NumberingFemDGQuad(Grid,OrdPoly,Proc)  

  DVZ = KernelAbstractions.zeros(backend,FT,size(DSZCPU))
  DVZCPU = 2 * DSZCPU
  DVZCPU[1,1] += 1 / wZCPU[1]
  DVZCPU[OrdPolyZ+1,OrdPolyZ+1] += -1 / wZCPU[OrdPolyZ+1]
  copyto!(DVZ,DVZCPU)
  DVZT=DVZ'

  Glob = KernelAbstractions.zeros(backend,Int,size(GlobCPU))
  GlobE = KernelAbstractions.zeros(backend,Int,size(GlobECPU))
  Stencil = KernelAbstractions.zeros(backend,Int,size(StencilCPU))
  copyto!(Stencil,StencilCPU)
  MasterSlave = KernelAbstractions.zeros(backend,Int,size(MasterSlaveCPU))
  copyto!(MasterSlave,MasterSlaveCPU)
  BoundaryDoF = KernelAbstractions.zeros(backend,Int,size(BoundaryDoFCPU))
  copyto!(BoundaryDoF,BoundaryDoFCPU)
  copyto!(Glob,GlobCPU)
  copyto!(GlobE,GlobECPU)
  return DGQuad{FT,
                 typeof(w),
                 typeof(DW),
                 typeof(MasterSlave),
                 typeof(Stencil),
                 typeof(GlobE)}(
    OrdPoly,
    OrdPolyZ,
    DoF,
    Glob,
    GlobE,
    Stencil,
    NumG,
    NumI,
    w,
    xw,
    xwCPU,
    xe,
    IntXE2F,
    wZ,
    xwZ,
    xwZCPU,
    IntZE2F,
    DW,
    DWT,
    DS,
    DST,
    DSZ,
    DWZ,
    DV,
    DVT,
    DVZ,
    DVZT,
    S,
    BoundaryDoF,
    MasterSlave,
 )
end
