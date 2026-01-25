abstract type DGElement end
mutable struct CGQuad{FT<:AbstractFloat,
                        AT1<:AbstractArray,
                        AT2<:AbstractArray,
                        AT3<:AbstractArray,
                        IT1<:AbstractArray,
                        IT2<:AbstractArray}
    DoFN::Int                    
    DoFE::Int                    
    DoF::Int
    OrdPoly::Int
    OrdPolyZ::Int
    Glob::IT2
    Stencil::IT2
    NumG::Int
    NumI::Int
    ksi::AT2
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
    InterOutputH::AT2
    InterOutputV::AT2
end

function CGQuad{FT}(backend,OrdPoly,OrdPolyZ,OrdPrint,Grid) where FT<:AbstractFloat
# Discretization
  OP=OrdPoly+1
  OPZ=OrdPolyZ+1
  nz = Grid.nz

# CG = CGStruct{FT}(backend)
  DoF = OP * OP
  DoFN = 1
  DoFE = OrdPoly - 1

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
  ksi = zeros(2,(OrdPoly+1)*(OrdPoly+1))
  iDoF = 1
  for j = 1 : OrdPoly + 1
    for i = 1 : OrdPoly +1
      ksi[1,iDoF] = xwCPU[i]
      ksi[2,iDoF] = xwCPU[j]
      iDoF += 1
    end
  end

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
  DoFPrint = (OrdPrint + 1) * (OrdPrint + 1)
  InterOutputHCPU = zeros(Float64,DoFPrint,DoF)
  dd = 2 / (OrdPrint + 1)
  eta0=-1.0
  iDoFP = 1
  for jRef = 1 : (OrdPrint + 1)
    ksi0 = -1.0
    eta1 = eta0 + dd
    for iRef = 1 : (OrdPrint + 1)
      ksi1 = ksi0 + dd
      iDoF = 1
      for j = 1 : OrdPoly + 1
        for i = 1 : OrdPoly + 1
          InterOutputHCPU[iDoFP,iDoF] = DG.Lagrange(0.5*(ksi0+ksi1),xwCPU,i)*
              DG.Lagrange(0.5*(eta0+eta1),xwCPU,j)
          iDoF += 1    
        end
      end
      iDoFP += 1
      ksi0 = ksi1
    end
    eta0 = eta1
  end
  InterOutputH = KernelAbstractions.zeros(backend,FT,size(InterOutputHCPU))
  copyto!(InterOutputH,InterOutputHCPU)
  InterOutputV = KernelAbstractions.ones(backend,FT,1,1)

  return CGQuad{FT,
                 typeof(w),
                 typeof(DW),
                 typeof(M),
                 typeof(MasterSlave),
                 typeof(Stencil)}(
    DoFN,
    DoFE,
    DoF,
    OrdPoly,
    OrdPolyZ,
    Glob,
    Stencil,
    NumG,
    NumI,
    ksi,
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
    InterOutputH,
    InterOutputV,
 )
end

mutable struct DGQuad{FT<:AbstractFloat,
                        AT1<:AbstractArray,
                        AT2<:AbstractArray,
                        IT1<:AbstractArray,
                        IT2<:AbstractArray,
                        IT3<:AbstractArray} <: DGElement 
    DoFN::Int                    
    DoFE::Int                    
    DoF::Int                    
    PosDoFE::IT2
    PosDoFECPU::Array{Int, 2}
    OrdPoly::Int
    OrdPolyZ::Int
    Glob::IT2
    GlobE::IT3
    Stencil::IT2
    NumG::Int
    NumI::Int
    ksi::Array{FT, 2}
    w::AT1
    wF::AT1
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
    InterOutputH::AT2
    InterOutputV::AT2
end

function DGQuad{FT}(backend,OrdPoly,OrdPolyZ,OrdPrint,OrdPrintZ,Grid,Proc) where FT<:AbstractFloat
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
  wFCPU = zeros(OP)
  @. wFCPU = wCPU[1]
  wF = KernelAbstractions.zeros(backend,FT,size(wFCPU))
  copyto!(wF,wFCPU)
  
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
  ksi = zeros(2,(OrdPoly+1)*(OrdPoly+1))
  iDoF = 1
  for j = 1 : OrdPoly + 1
    for i = 1 : OrdPoly +1
      ksi[1,iDoF] = xwCPU[i]
      ksi[2,iDoF] = xwCPU[j]
      iDoF += 1
    end
  end
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

  (DWCPU,DSCPU,DVCPU)=DG.DerivativeMatrixSingle(OrdPoly)
  DS = KernelAbstractions.zeros(backend,FT,size(DSCPU))
  copyto!(DS,DSCPU)
  DW = KernelAbstractions.zeros(backend,FT,size(DWCPU))
  copyto!(DW,DWCPU)
  DST=DS'
  DWT=DW'
  DV = KernelAbstractions.zeros(backend,FT,size(DVCPU))
  copyto!(DV,DVCPU)
  DVT=DV'

  Q = diagm(wCPU) * DSCPU
  S = Q - Q'

  (DWZCPU,DSZCPU,DVZCPU)=DG.DerivativeMatrixSingle(OrdPolyZ)
  DSZ = KernelAbstractions.zeros(backend,FT,size(DSZCPU))
  copyto!(DSZ,DSZCPU)
  DWZ = KernelAbstractions.zeros(backend,FT,size(DWZCPU))
  copyto!(DWZ,DWZCPU)
  DVZ = KernelAbstractions.zeros(backend,FT,size(DVZCPU))
  copyto!(DVZ,DVZCPU)
  DVZT=DVZ'

  (GlobCPU,GlobECPU,NumG,NumI,StencilCPU,MasterSlaveCPU,BoundaryDoFCPU) =
    NumberingFemDGQuad(Grid,OrdPoly,Proc)  

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
  DoFPrint = (OrdPrint + 1) * (OrdPrint + 1)
  InterOutputHCPU = zeros(Float64,DoFPrint,DoF)
  dd = 2 / (OrdPrint + 1)
  eta0=-1.0
  iDoFP = 1
  for jRef=1:OrdPrint + 1
    ksi0=-1.0
    eta1=eta0+dd
    for iRef=1:OrdPrint + 1
      ksi1=ksi0+dd
      iDoF = 1
      for j=1:OrdPoly+1
        for i=1:OrdPoly+1
          InterOutputHCPU[iDoFP,iDoF] = DG.Lagrange(0.5*(ksi0+ksi1),xwCPU,i)*
              DG.Lagrange(0.5*(eta0+eta1),xwCPU,j)
          iDoF += 1
        end
      end
      ksi0 = ksi1
      iDoFP += 1
    end
    eta0 = eta1
  end
  InterOutputH = KernelAbstractions.zeros(backend,FT,size(InterOutputHCPU))
  copyto!(InterOutputH,InterOutputHCPU)

  InterOutputVCPU = zeros(FT,OrdPrintZ+1,OrdPolyZ+1)
  dd = 2 / (OrdPrintZ + 1)
  zeta0 = -1.0
  for kRef = 1 : OrdPrintZ + 1
    zeta1 = zeta0 + dd
    for k = 1 : OrdPolyZ+1
      InterOutputVCPU[kRef,k] = DG.Lagrange(0.5*(zeta0+zeta1),xwZCPU,k)
    end
    zeta0 = zeta1
  end
  InterOutputV = KernelAbstractions.zeros(backend,FT,size(InterOutputVCPU))
  copyto!(InterOutputV,InterOutputVCPU)
  DoFN = 0
  DoFE = OrdPoly + 1
  DoF = DoFE * DoFE 
  PosDoFECPU = zeros(Int,DoFE,4)
  for i = 1 : DoFE
    PosDoFECPU[i,1] = i  
    PosDoFECPU[i,2] = i * DoFE
    PosDoFECPU[i,3] = i + (DoFE - 1) *  DoFE
    PosDoFECPU[i,4] = 1 + (i - 1) * DoFE
  end  
  PosDoFE = KernelAbstractions.zeros(backend,Int,size(PosDoFECPU))
  copyto!(PosDoFE,PosDoFECPU)

  return DGQuad{FT,
                 typeof(w),
                 typeof(DW),
                 typeof(MasterSlave),
                 typeof(PosDoFE),
                 typeof(GlobE)}(
    DoFN,
    DoFE,
    DoF,
    PosDoFE,
    PosDoFECPU,
    OrdPoly,
    OrdPolyZ,
    Glob,
    GlobE,
    Stencil,
    NumG,
    NumI,
    ksi,
    w,
    wF,
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
    InterOutputH,
    InterOutputV,
 )
end

mutable struct DGTri{FT<:AbstractFloat,
                        AT1<:AbstractArray,
                        AT2<:AbstractArray,
                        IT2<:AbstractArray,
                        IT3<:AbstractArray} <: DGElement 
    k::Int                   
    DoFN::Int                   
    DoFE::Int                   
    DoF::Int
    OrdPolyZ::Int
    NumG::Int
    NumI::Int
    ksiCPU::Array{FT, 2}
    ksi::AT2
    Glob::IT2
    GlobE::IT3
    Dx1::AT2
    Dx2::AT2
    DSx1::AT2
    DSx2::AT2
    DVZ::AT2
    DVZT::AT2
    w::AT1
    wF::AT1
    PosDoFE::IT2
    PosDoFECPU::Array{Int, 2}
    wZ::AT1
    xwZ::AT1
    xwZCPU::Array{FT, 1}
    DSZ::AT2
    DWZ::AT2
    PL2::AT2
    PL2CPU::Array{FT, 2}
    phi
    InterOutputH::AT2
    InterOutputV::AT2
end    

function DGTri{FT}(backend,Method,OrdPolyZ,OrdPrint,OrdPrintZ,Grid,Proc) where FT<:AbstractFloat

  N1 = [0.0 -1.0]
  N2 = [1.0 1.0]
  N3 = [-1.0 0.0]
  if Method  == "Kubatko0"
    k = 0  
    n = 3
    ksi = zeros(2,n)
    N = zeros(2,n)
    wFR = zeros(n)
    ksi[:,1] = [0.0  -1]
    ksi[:,2] = [0.0 0.0]
    ksi[:,3] = [-1.0 0.0]
    N[:,1] = N1
    N[:,2] = N2
    N[:,3] = N3
    w = ones(n) * 2/3
    _, wF = gausslegendre(k+1)
    wFR[1:k+1] = wF
    wFR[k+1+1:2*(k+1)] = wF
    wFR[2*(k+1)+1:3*(k+1)] = wF
    DoFE = 1
    PosDoFE = zeros(Int,DoFE,3)
    PosDoFE[1,1] = 1
    PosDoFE[1,2] = 2 
    PosDoFE[1,3] = 3 
    @. wF = w[1:k+1] / wF
  elseif Method  == "Kubatko1"
    k = 1  
    n = 6
    ksi = zeros(2,n)
    N = zeros(2,n)
    wFR = zeros(n)
    ksi[:,1] = [-0.577350269189626  -1]
    ksi[:,2] = [0.577350269189626 -1]
    ksi[:,3] = [-1 -0.577350269189626]
    ksi[:,4] = [-1 0.577350269189626]
    ksi[:,5] = [-0.577350269189626 0.577350269189626]
    ksi[:,6] = [ 0.577350269189626 -0.577350269189626]

    N[:,1] = N1
    N[:,2] = N1
    N[:,3] = N3
    N[:,4] = N3
    N[:,5] = N2
    N[:,6] = N2

    w = ones(n) * 1/3
    _, wF = gausslegendre(k+1)
    wFR[1:k+1] = wF
    wFR[k+1+1:2*(k+1)] = wF
    wFR[2*(k+1)+1:3*(k+1)] = wF
    DoFE = 2
    PosDoFE = zeros(Int,DoFE,3)
    PosDoFE[1,1] = 1
    PosDoFE[2,1] = 2 
    PosDoFE[1,2] = 6 
    PosDoFE[2,2] = 5 
    PosDoFE[1,3] = 3 
    PosDoFE[2,3] = 4 
    @. wF = w[1:k+1] / wF
  elseif Method  == "Kubatko2LGL"
    k = 2
    n = 10
    ksi = zeros(2,n)
    N = zeros(2,n)
    w = zeros(n)
    wFR = zeros(n)

    ksi[:,1] = [-1 -1] 
    ksi[:,2] = [-0.447213595499958 -1] 
    ksi[:,3] = [0.447213595499958 -1] 
    ksi[:,4] = [1 -1] 
    ksi[:,5] = [0.447213595499958 -0.447213595499958]
    ksi[:,6] = [-0.447213595499958 0.447213595499958]
    ksi[:,7] = [-1 1] 
    ksi[:,8] = [-1 0.447213595499958] 
    ksi[:,9] = [-1 -0.447213595499958] 
    ksi[:,10] = [-0.333333333333333 -0.333333333333333]

    N[:,1] = N1 + N3
    N[:,2] = N1
    N[:,3] = N1
    N[:,4] = N1 + N2
    N[:,5] = N2
    N[:,6] = N2
    N[:,7] = N2 + N3
    N[:,8] = N3
    N[:,9] = N3

    w[1] = 0.033333333333333 
    w[2] = 0.166666666666667
    w[3] = 0.166666666666667
    w[4] = 0.033333333333333 
    w[5] = 0.166666666666667
    w[6] = 0.166666666666667
    w[7] = 0.033333333333333 
    w[8] = 0.166666666666667
    w[9] = 0.166666666666667
    w[10] = 0.9

    _, wF = gausslobatto(k+2)
    @. wFR[1:4] = wF
    @. wFR[4:7] = wF
    @. wFR[7:9] = wF[1:3]
    wFR[1] = wF[4]

    DoFE = 4
    PosDoFE = zeros(Int,DoFE,3)

    PosDoFE[1,1] = 1
    PosDoFE[2,1] = 2
    PosDoFE[3,1] = 3
    PosDoFE[4,1] = 4

    PosDoFE[1,2] = 4
    PosDoFE[2,2] = 5
    PosDoFE[3,2] = 6
    PosDoFE[4,2] = 7

    PosDoFE[1,3] = 7
    PosDoFE[2,3] = 8
    PosDoFE[3,3] = 9
    PosDoFE[4,3] = 1

    wE = [w[1], w[2], w[3], w[4]] 
    @. wF = wE / wF

  elseif Method  == "Kubatko2"
    k = 2
    n = 10
    ksi = zeros(2,n)
    N = zeros(2,n)
    w = zeros(n)
    wFR = zeros(n)

    ksi[:,1] = [-0.774596669241483 -1] 
    ksi[:,2] = [0 -1] 
    ksi[:,3] = [0.774596669241483 -1] 

    ksi[:,4] = [0.774596669241483 -0.774596669241483] 
    ksi[:,5] = [0 0] 
    ksi[:,6] = [-0.774596669241483 0.774596669241483] 

    ksi[:,7] = [-1 0.774596669241483]
    ksi[:,8] = [-1 0] 
    ksi[:,9] = [-1 -0.774596669241483] 

    ksi[:,10] = [-0.333333333333333 -0.333333333333333]
    N[:,1] = N1
    N[:,2] = N1
    N[:,3] = N1
    N[:,4] = N2
    N[:,5] = N2
    N[:,6] = N2
    N[:,7] = N3
    N[:,8] = N3
    N[:,9] = N3

    w[1] = 0.083333333333333 
    w[2] = 0.2 
    w[3] = 0.083333333333333 
    w[4] = 0.083333333333333 
    w[5] = 0.2
    w[6] = 0.083333333333333 
    w[7] = 0.083333333333333 
    w[8] = 0.2 
    w[9] = 0.083333333333333 
    w[10] = 0.9


    _, wF = gausslegendre(k+1)
    wFR[1:k+1] = wF
    wFR[k+1+1:2*(k+1)] = wF
    wFR[2*(k+1)+1:3*(k+1)] = wF
    DoFE = k + 1

    PosDoFE = zeros(Int,DoFE,3)
    PosDoFE[1,1] = 1
    PosDoFE[2,1] = 2
    PosDoFE[3,1] = 3

    PosDoFE[1,2] = 4
    PosDoFE[2,2] = 5
    PosDoFE[3,2] = 6

    PosDoFE[1,3] = 7
    PosDoFE[2,3] = 8
    PosDoFE[3,3] = 9

    @. wF = w[1:k+1] / wF
  elseif Method  == "Kubatko5"
    k = 5
    n = 34
    ksi = zeros(2,n)
    N = zeros(2,n)
    wFR = zeros(n)
    for i = 1 : k + 1
      N[:,i] = N1
      N[:,i+k+1] = N2
      N[:,i+2*(k+1)] = N3
    end  

    ksiT = [
     -0.932469514203153  -1.000000000000000
     -0.661209386466264  -1.000000000000000
     -0.238619186083197  -1.000000000000000
      0.238619186083197  -1.000000000000000
      0.661209386466263  -1.000000000000000
      0.932469514203152  -1.000000000000000

      0.932469514203153  -0.932469514203153
      0.661209386466264  -0.661209386466264
      0.238619186083197  -0.238619186083197
     -0.238619186083197   0.238619186083197
     -0.661209386466263   0.661209386466263
     -0.932469514203152   0.932469514203152

     -1.000000000000000   0.932469514203153
     -1.000000000000000   0.661209386466264
     -1.000000000000000   0.238619186083197
     -1.000000000000000  -0.238619186083197
     -1.000000000000000  -0.661209386466263
     -1.000000000000000  -0.932469514203152

     -0.333333333333333  -0.333333333333333
     -0.550932095205921  -0.550932095205921
     -0.550932095205921   0.101864190411842
      0.101864190411842  -0.550932095205921
      0.611382635695146  -0.744884752682278
     -0.744884752682278   0.611382635695146
      0.611382635695146  -0.866497883012868
     -0.866497883012868   0.611382635695146
     -0.744884752682278  -0.866497883012868
     -0.866497883012868  -0.744884752682278
      0.154088161163644  -0.813213547833420
     -0.813213547833420   0.154088161163644
      0.154088161163644  -0.340874613330224
     -0.340874613330224   0.154088161163644
     -0.813213547833420  -0.340874613330224
     -0.340874613330224  -0.813213547833420]

    for i = 1 : n
      ksi[1,i] = ksiT[i,1]
      ksi[2,i] = ksiT[i,2]
    end  


    w = [
      0.006658503841347
      0.016313739546177
      0.025165807590991
      0.025165807590991
      0.016313739546177
      0.006658503841347
      0.006658503841347
      0.016313739546177
      0.025165807590991
      0.025165807590991
      0.016313739546177
      0.006658503841347
      0.006658503841347
      0.016313739546177
      0.025165807590991
      0.025165807590991
      0.016313739546177
      0.006658503841347
      0.209054393645172
      0.095225751929022
      0.095225751929022
      0.095225751929022
      0.060381677758004
      0.060381677758004
      0.060381677758004
      0.060381677758004
      0.060381677758004
      0.060381677758004
      0.142358329691442
      0.142358329691442
      0.142358329691442
      0.142358329691442
      0.142358329691442
      0.142358329691442]

    _, wF = gausslegendre(k+1)
    wFR[1:k+1] = wF
    wFR[k+1+1:2*(k+1)] = wF
    wFR[2*(k+1)+1:3*(k+1)] = wF
    DoFE = k + 1
    PosDoFE = zeros(Int,DoFE,3)
    for i = 1 : k + 1
      PosDoFE[i,1] = i
      PosDoFE[i,2] = i + k + 1 
      PosDoFE[i,3] = i + 2 * (k + 1) 
    end  
    @. wF = w[1:k+1] / wF
    
  elseif Method == "Hicken1"
    k = 1
    n = 7
    ksi = zeros(2,n)
    N = zeros(2,n)
    w = zeros(n)
    wF = zeros(n)
    wFx1 = zeros(n)
    wFx1 = zeros(n)
    ksi[:,1] = [-1.0 -1.0]
    N[:,1] = [0.0 -1.0] + [-1.0 0.0]
    ksi[:,2] = [1.0 -1.0]
    N[:,2] = [0.0 -1.0] + [1.0 1.0]
    ksi[:,3] = [-1.0 1.0]
    N[:,3] = [1.0 1.0] + [-1.0 0.0]
    ksi[:,4] = [0.0 -1.0]
    N[:,4] = [0.0 -1.0]
    ksi[:,5] = [0.0 0.0]
    N[:,5] = [1.0 1.0]
    ksi[:,6] = [-1.0 0.0]
    N[:,6] = [-1.0 0.0]
    ksi[:,7] = [-1/3 -1/3]

    w[1] = 0.09999999999999999
    w[2] =0.09999999999999999
    w[3] =0.09999999999999999
    w[4] =0.26666666666666666
    w[5] =0.26666666666666666
    w[6] =0.26666666666666666
    w[7] =0.9000000000000002
    wFR = zeros(n)
    wFR[1] = 0.333333333333333 
    wFR[2] = 0.333333333333333 
    wFR[3] = 0.333333333333333 
    wFR[4] = 1.333333333333333
    wFR[5] = 1.333333333333333
    wFR[6] = 1.333333333333333
    DoFE = 3
    PosDoFE = zeros(Int,DoFE,3)
    PosDoFE[1,1] = 1
    PosDoFE[2,1] = 4 
    PosDoFE[3,1] = 2 
    PosDoFE[1,2] = 2 
    PosDoFE[2,2] = 5 
    PosDoFE[3,2] = 3 
    PosDoFE[1,3] = 1 
    PosDoFE[2,3] = 6 
    PosDoFE[3,3] = 3 
  end

  DoF = n
  wFx1 = wFR .* N[1,:]
  wFx2 = wFR .* N[2,:]
  s = @polyvar x[1:2]
  phi = DG.Polynomial_k(k,s)
  nSt = size(phi,1)
  phiDx1 = Array{Polynomial,2}(undef,nSt,1)
  phiDx2 = Array{Polynomial,2}(undef,nSt,1)
  for i = 1 : nSt
    phiDx1[i] = differentiate(phi[i],x[1])  
    phiDx2[i] = differentiate(phi[i],x[2])  
  end  
  V = zeros(n,nSt)
  VDx1 = zeros(n,nSt)
  VDx2 = zeros(n,nSt)
  for j = 1 : nSt
    for i = 1 : n
      V[i,j] = phi[j](ksi[1,i],ksi[2,i])  
      VDx1[i,j] = phiDx1[j](ksi[1,i],ksi[2,i])  
      VDx2[i,j] = phiDx2[j](ksi[1,i],ksi[2,i])  
    end
  end  

  Q,R = qr(diagm(sqrt.(w)) * V)
  PL2 = R \ (Q[:,1:nSt]' * diagm(sqrt.(w)))
  VO = V / R 
  VODx1 = VDx1 / R 
  VODx2 = VDx2 / R 

  Dx1CPU = 0.5 * (diagm(1.0 ./ w) + VO * VO') * diagm(wFx1) * (I - VO * VO' * diagm(w)) +
    VODx1 * VO' * diagm(w)
  Dx2CPU = 0.5 * (diagm(1.0 ./ w) + VO * VO') * diagm(wFx2) * (I - VO * VO' * diagm(w)) +
    VODx2 * VO' * diagm(w)

  Dx1 = KernelAbstractions.zeros(backend,FT,size(Dx1CPU))
  copyto!(Dx1,Dx1CPU)
  Dx2 = KernelAbstractions.zeros(backend,FT,size(Dx2CPU))
  copyto!(Dx2,Dx2CPU)
  for i = 1 : n
    Dx1CPU[i,i] = 0.0  
    Dx2CPU[i,i] = 0.0  
  end
  Dx1CPU .*= 2.0
  Dx2CPU .*= 2.0
  DSx1 = KernelAbstractions.zeros(backend,FT,size(Dx1CPU))
  copyto!(DSx1,Dx1CPU)
  DSx2 = KernelAbstractions.zeros(backend,FT,size(Dx2CPU))
  copyto!(DSx2,Dx2CPU)

  Stencil = zeros(Int,0,0)
  NumG = n * Grid.NumFaces
  NumI = NumG

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
  (DWZCPU,DSZCPU)=DG.DerivativeMatrixSingle(OrdPolyZ)
  DSZ = KernelAbstractions.zeros(backend,FT,size(DSZCPU))
  copyto!(DSZ,DSZCPU)
  DWZ = KernelAbstractions.zeros(backend,FT,size(DWZCPU))
  copyto!(DWZ,DWZCPU)
  DVZ = KernelAbstractions.zeros(backend,FT,size(DSZCPU))
  DVZCPU = 2 * DSZCPU
  DVZCPU[1,1] += 1 / wZCPU[1]
  DVZCPU[OrdPolyZ+1,OrdPolyZ+1] += -1 / wZCPU[OrdPolyZ+1]
  copyto!(DVZ,DVZCPU)
  DVZT=DVZ'

  GlobCPU,GlobECPU,NumG,NumI,Stencil,MasterSlave,BoundaryDoF = NumberingFemDGTri(Grid,n,DoFE,PosDoFE,Proc)
  Glob = KernelAbstractions.zeros(backend,Int,size(GlobCPU))
  GlobE = KernelAbstractions.zeros(backend,Int,size(GlobECPU))
  copyto!(Glob,GlobCPU)
  copyto!(GlobE,GlobECPU)
  _,RefineMidPoints = Grids.PrintPoints(OrdPrint,Grids.Tri())

  DoFPrint = size(RefineMidPoints,1)
  InterOutputHCPU = zeros(Float64,DoFPrint,DoF)
  for iM = 1 : DoFPrint
    ksi1 = RefineMidPoints[iM,1]  
    ksi2 = RefineMidPoints[iM,2]  
    for iDoF = 1 : DoF
      @views XR = PL2[:,iDoF] 
      for i = 1 : size(XR,1)
        InterOutputHCPU[iM,iDoF] += phi[i](ksi1,ksi2) * PL2[i,iDoF]
      end
    end
  end
  InterOutputH = KernelAbstractions.zeros(backend,FT,size(InterOutputHCPU))
  copyto!(InterOutputH,InterOutputHCPU)

  InterOutputVCPU = zeros(Float64,OrdPrintZ+1,OrdPolyZ+1)
  dd = 2 / (OrdPrintZ + 1)
  zeta0 = -1.0
  for kRef = 1 : OrdPrintZ + 1
    zeta1 = zeta0 + dd
    for k = 1 : OrdPolyZ+1
      InterOutputVCPU[kRef,k] = DG.Lagrange(0.5*(zeta0+zeta1),xwZCPU,k)
    end
    zeta0 = zeta1
  end
  InterOutputV = KernelAbstractions.zeros(backend,FT,size(InterOutputVCPU))
  copyto!(InterOutputV,InterOutputVCPU)
  DoFN = 0
  ksiGPU = KernelAbstractions.zeros(backend,FT,size(ksi))
  copyto!(ksiGPU,ksi)
  wGPU = KernelAbstractions.zeros(backend,FT,size(w))
  copyto!(wGPU,w)
  PL2GPU = KernelAbstractions.zeros(backend,FT,size(PL2))
  copyto!(PL2GPU,PL2)
  PosDoFEGPU = KernelAbstractions.zeros(backend,FT,size(PosDoFE))
  copyto!(PosDoFEGPU,PosDoFE)

  return DGTri{FT,
               typeof(wGPU),
               typeof(Dx1),
               typeof(Glob),
               typeof(GlobE)}(
    k,           
    DoFN,
    DoFE,
    DoF,
    OrdPolyZ,
    NumG,
    NumI,
    ksi,
    ksiGPU,
    Glob,
    GlobE,
    Dx1,
    Dx2,
    DSx1,
    DSx2,
    DVZ,
    DVZT,
    wGPU,
    wF,
    PosDoFEGPU,
    PosDoFE,
    wZ,
    xwZ,
    xwZCPU,
    DSZ,
    DWZ,
    PL2GPU,
    PL2,
    phi,
    InterOutputH,
    InterOutputV,
  )
end

mutable struct DG1{FT<:AbstractFloat,
                        AT1<:AbstractArray,
                        AT2<:AbstractArray} <: DGElement 
    OrdPolyZ::Int                    
    wZ::AT1
    xwZ::AT1
    xwZCPU::Array{FT, 1}
    DSZ::AT2
    DWZ::AT2
    DVZ::AT2
    DVZT::AT2
end    
function DG1{FT}(backend,OrdPolyZ,OrdPrintZ) where FT<:AbstractFloat

  DoF = OrdPolyZ

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

  (DWZCPU,DSZCPU)=DG.DerivativeMatrixSingle(OrdPolyZ)
  DSZ = KernelAbstractions.zeros(backend,FT,size(DSZCPU))
  copyto!(DSZ,DSZCPU)
  DWZ = KernelAbstractions.zeros(backend,FT,size(DWZCPU))
  copyto!(DWZ,DWZCPU)
  DVZ = KernelAbstractions.zeros(backend,FT,size(DSZCPU))
  DVZCPU = 2 * DSZCPU
  DVZCPU[1,1] += 1 / wZCPU[1]
  DVZCPU[OrdPolyZ+1,OrdPolyZ+1] += -1 / wZCPU[OrdPolyZ+1]
  copyto!(DVZ,DVZCPU)
  DVZT=DVZ'

  return DG1{FT,
             typeof(wZ),
             typeof(DVZ)}(
    OrdPolyZ,         
    wZ,
    xwZ,
    xwZCPU,
    DSZ,
    DWZ,
    DVZ,
    DVZT,
  )
end

