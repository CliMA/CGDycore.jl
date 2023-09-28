mutable struct RT0Struct{FT<:AbstractFloat,
                        IT2<:AbstractArray}
  Glob::IT2                        
  Stencil::IT2
end
function RT0Struct{FT}(backend) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  Stencil = KernelAbstractions.zeros(backend,Int,0,0)
 return CGStruct{FT,
                 typeof(Glob)}( 
    Glob,
    Stencil,
  )
end

mutable struct CGStruct{FT<:AbstractFloat,
                        AT2<:AbstractArray,
                        IT1<:AbstractArray,
                        IT2<:AbstractArray}
    OrdPoly::Int
    OrdPolyZ::Int
    DoF::Int
    Glob::IT2
    Stencil::IT2
    NumG::Int
    NumI::Int
    w::Array{FT, 1}
    xw::Array{FT, 1}
    xe::Array{FT, 1}
    IntXE2F::Array{FT, 2}
    xwZ::Array{FT, 1}
    IntZE2F::Array{FT, 2}
    DW::AT2
    DWT::Array{FT, 2}
    DS::AT2
    DST::Array{FT, 2}
    DSZ::Array{FT, 2}
    S::Array{FT, 2}
    M::AT2
    MMass::AT2
    MW::AT2
    Boundary::Array{Int, 1}
    MasterSlave::IT1
end
function CGStruct{FT}(backend) where FT<:AbstractFloat
  OrdPoly = 0
  OrdPolyZ = 0
  DoF = 0
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  Stencil = KernelAbstractions.zeros(backend,Int,0,0)
  NumG = 0
  NumI = 0
  w = zeros(FT,0)
  xw = zeros(FT,0)
  xe = zeros(FT,0)
  IntXE2F = zeros(FT,0,0)
  xwZ = zeros(FT,0)
  IntZE2F = zeros(FT,0,0)
  DW = KernelAbstractions.zeros(backend,FT,0,0)
  DWT = zeros(FT,0,0)
  DS = KernelAbstractions.zeros(backend,FT,0,0)
  DST = zeros(FT,0,0)
  DSZ = zeros(FT,0,0)
  S = zeros(FT,0,0)
  M = KernelAbstractions.zeros(backend,FT,0,0)
  MMass = KernelAbstractions.zeros(backend,FT,0,0)
  MW = KernelAbstractions.zeros(backend,FT,0,0)
  Boundary = zeros(Int,0)
  MasterSlave = KernelAbstractions.zeros(backend,Int,0)
 return CGStruct{FT,
                 typeof(DW),
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
    xe,
    IntXE2F,
    xwZ,
    IntZE2F,
    DW,
    DWT,
    DS,
    DST,
    DSZ,
    S,
    M,
    MMass,
    MW,
    Boundary,
    MasterSlave,
 )
end 

function CGStruct(backend,FT,OrdPoly,OrdPolyZ,Grid) 
# Discretization
  OP=OrdPoly+1
  OPZ=OrdPolyZ+1

  CG = CGStruct{FT}(backend)
  CG.OrdPoly=OrdPoly
  CG.OrdPolyZ=OrdPolyZ
  CG.DoF = OP * OP

  (CG.w,CG.xw)=GaussLobattoQuad(CG.OrdPoly)
  (wZ,CG.xwZ)=GaussLobattoQuad(CG.OrdPolyZ)
  CG.xe = zeros(OrdPoly+1)
  CG.xe[1] = -1.0
  for i = 2 : OrdPoly
    CG.xe[i] = CG.xe[i-1] + 2.0/OrdPoly
  end
  CG.xe[OrdPoly+1] = 1.0

  CG.IntXE2F = zeros(OrdPoly+1,OrdPoly+1)
  for j = 1 : OrdPoly + 1
    for i = 1 : OrdPoly +1
      CG.IntXE2F[i,j] = Lagrange(CG.xw[i],CG.xe,j)
    end
  end

  CG.IntZE2F = zeros(OrdPolyZ+1,OrdPolyZ+1)
  for j = 1 : OrdPolyZ + 1
    for i = 1 : OrdPolyZ +1
      CG.IntZE2F[i,j] = Lagrange(CG.xwZ[i],CG.xwZ,j)
    end
  end

  (DW,DS)=DerivativeMatrixSingle(CG.OrdPoly)
  CG.DS = KernelAbstractions.zeros(backend,FT,size(DS))
  copyto!(CG.DS,DS)
  CG.DW = KernelAbstractions.zeros(backend,FT,size(DW))
  copyto!(CG.DW,DW)
  CG.DST=DS'
  CG.DWT=DW'

  Q = diagm(CG.w) * DS
  CG.S = Q - Q'
  (DWZ,CG.DSZ)=DerivativeMatrixSingle(CG.OrdPolyZ)
  (Glob,CG.NumG,CG.NumI,Stencil,MasterSlave) =
    NumberingFemCG(Grid,OrdPoly)  

  CG.Glob = KernelAbstractions.zeros(backend,Int,size(Glob))
  copyto!(CG.Glob,Glob)
  CG.Stencil = KernelAbstractions.zeros(backend,Int,size(Stencil))
  copyto!(CG.Stencil,Stencil)
  CG.MasterSlave = KernelAbstractions.zeros(backend,Int,size(MasterSlave))
  copyto!(CG.MasterSlave,MasterSlave)


# Boundary nodes
  for iF = 1 : Grid.NumFaces
    Side = 0
    for iE in Grid.Faces[iF].E
       Side += 1 
       if Grid.Edges[iE].F[1] == 0 || Grid.Edges[iE].F[2] == 0
          @views GlobLoc = reshape(Glob[:,iF],OP,OP) 
         if Side == 1
           for i in GlobLoc[1:OP-1,1]   
             push!(CG.Boundary,i)
           end
         elseif Side == 2  
           for i in GlobLoc[OP,1:OP-1]
             push!(CG.Boundary,i)
           end
         elseif Side == 3  
           for i in GlobLoc[2:OP,OP]
             push!(CG.Boundary,i)
           end
         elseif Side == 4  
           for i in GlobLoc[1,2:OP]
             push!(CG.Boundary,i)
           end
         end  
       end  
    end
  end  
  return CG
end
