abstract type DGElement end
abstract type CGElement end
mutable struct CGTri{FT<:AbstractFloat,
                        AT2<:AbstractArray} <: CGElement
  DoF::Int                        
  Dx1::Array{FT, 2}
  Dx2::Array{FT, 2}
end                        

function CGTri{FT}(backend,OrdPoly,OrdPolyZ,OrdPrint,Grid) where FT<:AbstractFloat
  DoF = 0
  Dx1 = zeros(FT,0,0)
  Dx2 = zeros(FT,0,0)
  return CGTri{FT,
                 typeof(Dx1)}(
    DoF,
    Dx1,
    Dx2
  )
end

mutable struct CGQuad{FT<:AbstractFloat,
                        AT1<:AbstractArray,
                        AT2<:AbstractArray,
                        IT1<:AbstractArray,
                        IT2<:AbstractArray} <: CGElement
    DoFN::Int                    
    DoFE::Int                    
    DoF::Int
    OrdPoly::Int
    OrdPolyZ::Int
    Mz::Int
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
    DS::AT2
    DST::Array{FT, 2}
    DSZ::AT2
    S::Array{FT, 2}
    BoundaryDoF::Array{Int, 1}
    MasterSlave::IT1
    InterOutputH::AT2
    InterOutputV::AT2
end

function CGQuad{FT}(backend,OrdPoly,OrdPolyZ,OrdPrint,Grid) where FT<:AbstractFloat
# Discretization
  OP=OrdPoly+1
  OPZ=OrdPolyZ+1
  Mz = OrdPolyZ
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
  ksiGPU = KernelAbstractions.zeros(backend,FT,size(ksi))
  copyto!(ksiGPU,ksi)

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
                 typeof(MasterSlave),
                 typeof(Stencil)}(
    DoFN,
    DoFE,
    DoF,
    OrdPoly,
    OrdPolyZ,
    Mz,
    Glob,
    Stencil,
    NumG,
    NumI,
    ksiGPU,
    w,
    xw,
    xwCPU,
    xe,
    IntXE2F,
    xwZ,
    xwZCPU,
    IntZE2F,
    DW,
    DS,
    DST,
    DSZ,
    S,
    BoundaryDoF,
    MasterSlave,
    InterOutputH,
    InterOutputV,
 )
end


