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


# Boundary nodes
#=
  Boundary = zeros(Int,0)
  for iF = 1 : Grid.NumFaces
    Side = 0
    for iE in Grid.Faces[iF].E
       Side += 1 
       if Grid.Edges[iE].F[1] == 0 || Grid.Edges[iE].F[2] == 0
         @views GlobLoc = reshape(GlobCPU[:,iF],OP,OP) 
         if Side == 1
           for i in GlobLoc[1:OP-1,1]   
             push!(Boundary,i)
           end
         elseif Side == 2  
           for i in GlobLoc[OP,1:OP-1]
             push!(Boundary,i)
           end
         elseif Side == 3  
           for i in GlobLoc[2:OP,OP]
             push!(Boundary,i)
           end
         elseif Side == 4  
           for i in GlobLoc[1,2:OP]
             push!(Boundary,i)
           end
         end  
       end  
    end
  end  
=#  
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
    Glob::IT3
    IndE::IT2
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
    S::Array{FT, 2}
    BoundaryDoF::Array{Int, 1}
    MasterSlave::IT1
    VZ::AT1
end

function DGQuad{FT}(backend,OrdPoly,OrdPolyZ,Grid) where FT<:AbstractFloat
# Discretization
  OP=OrdPoly+1
  OPZ=OrdPolyZ+1
  nz = Grid.nz

# CG = CGStruct{FT}(backend)
  DoF = OP * OP

  xwCPU, wCPU = gausslobatto(OrdPoly+1)
  w = KernelAbstractions.zeros(backend,FT,size(wCPU))
  xw = KernelAbstractions.zeros(backend,FT,size(xwCPU))
  copyto!(w,wCPU)
  copyto!(xw,xwCPU)
  
  if OrdPolyZ == 0
    xwZCPU = zeros(1)  
    wZ = 2 * ones(1)
  else    
    xwZCPU, wZ = gausslobatto(OrdPolyZ+1)
  end  
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
  (GlobCPU,NumG,NumI,StencilCPU,MasterSlaveCPU,BoundaryDoFCPU) =
    NumberingFemDGQuad(Grid,OrdPoly)  

  Glob = KernelAbstractions.zeros(backend,Int,size(GlobCPU))
  Stencil = KernelAbstractions.zeros(backend,Int,size(StencilCPU))
  copyto!(Stencil,StencilCPU)
  MasterSlave = KernelAbstractions.zeros(backend,Int,size(MasterSlaveCPU))
  copyto!(MasterSlave,MasterSlaveCPU)
  BoundaryDoF = KernelAbstractions.zeros(backend,Int,size(BoundaryDoFCPU))
  copyto!(BoundaryDoF,BoundaryDoFCPU)
  IndECPU = zeros(Int,4,OrdPoly+1)
  @inbounds for i = 1 : OrdPoly + 1
    IndECPU[1,i] = i                               #  1  2  3  4  5
    IndECPU[2,i] = i * (OrdPoly + 1)               #  5 10 15 20 25
    IndECPU[3,i] = i  + OrdPoly * (OrdPoly + 1)    # 21 22 23 24 25
    IndECPU[4,i] = 1  + (i-1) * (OrdPoly + 1)      #  1  6 11 16 21
  end
  IndE = KernelAbstractions.zeros(backend,Int,4,OrdPoly+1)
  copyto!(IndE,IndECPU)
  VZCPU = zeros(4)
  VZCPU[1] = -1.0
  VZCPU[2] = -1.0
  VZCPU[3] = 1.0
  VZCPU[4] = 1.0
  VZ = KernelAbstractions.zeros(backend,FT,4)
  copyto!(VZ,VZCPU)

# Boundary nodes
#=
  Boundary = zeros(Int,0)
  for iF = 1 : Grid.NumFaces
    Side = 0
    for iE in Grid.Faces[iF].E
       Side += 1 
       if Grid.Edges[iE].F[1] == 0 || Grid.Edges[iE].F[2] == 0
         @views GlobLoc = reshape(GlobCPU[:,iF],OP,OP) 
         if Side == 1
           for i in GlobLoc[1:OP-1,1]   
             push!(Boundary,i)
           end
         elseif Side == 2  
           for i in GlobLoc[OP,1:OP-1]
             push!(Boundary,i)
           end
         elseif Side == 3  
           for i in GlobLoc[2:OP,OP]
             push!(Boundary,i)
           end
         elseif Side == 4  
           for i in GlobLoc[1,2:OP]
             push!(Boundary,i)
           end
         end  
       end  
    end
  end  
=#  
  copyto!(Glob,GlobCPU)
  return DGQuad{FT,
                 typeof(w),
                 typeof(DW),
                 typeof(MasterSlave),
                 typeof(Stencil),
                 typeof(Glob)}(
    OrdPoly,
    OrdPolyZ,
    DoF,
    Glob,
    IndE,
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
    S,
    BoundaryDoF,
    MasterSlave,
    VZ,
 )
end
