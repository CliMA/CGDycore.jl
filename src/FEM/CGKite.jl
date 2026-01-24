mutable struct CG0KitePrimalStruct{FT<:AbstractFloat,
                        IT2<:AbstractArray} <: ScalarKitePElement
  Glob::IT2
  DoF::Int
  Comp::Int                      
  phi::Array{Polynomial,2}                       
  NumG::Int
  NumI::Int
  Type::Grids.ElementType
  M::AbstractSparseMatrix
  LUM::SparseArrays.UMFPACK.UmfpackLU{Float64, Int64}
end

function CG0KitePrimalStruct{FT}(::Grids.Quad,backend,Grid) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  Type = Grids.QuadPrimal()
  DoF = 1
  Comp = 1
  @polyvar x1 x2
  phi = Array{Polynomial,2}(undef,DoF,Comp)
  x, w = gaussradau(1)
  phi[1,1] = 1 + 0*x1 + 0*x2
  
  NumNodes = Grid.NumNodes
  NumFaces = Grid.NumFaces
  Nodes = Grid.Nodes


  iKite = 1
  iOff = 1
  NumG = 0
  Glob = KernelAbstractions.zeros(backend,Int,DoF,NumFaces)
  GlobCPU = zeros(Int,DoF,NumFaces)
  for iN = 1 : NumNodes
    if Nodes[iN].Type == 'F'
      NumF = length(Nodes[iN].F)
      Offset = 1 
      for i = 1 : NumF
        GlobCPU[1,iKite] = iOff
        iKite +=1 
      end 
      iOff += Offset
      NumG += Offset
    end 
  end
  NumI = NumG
  copyto!(Glob,GlobCPU)
  M = sparse([1],[1],[1.0])
  LUM = lu(M)
  return CG0KitePrimalStruct{FT,
                  typeof(Glob)}( 
    Glob,
    DoF,
    Comp,
    phi,                      
    NumG,
    NumI,
    Type,
    M,
    LUM,
      )
end

mutable struct CG1KitePrimalStruct{FT<:AbstractFloat,
                        IT2<:AbstractArray} <: ScalarKitePElement
  Glob::IT2
  DoF::Int
  Comp::Int                      
  phi::Array{Polynomial,2}                       
  rotphi::Array{Polynomial,2}                       
  Gradphi::Array{Polynomial,3}                       
  NumG::Int
  NumI::Int
  Type::Grids.ElementType
  points::Array{Float64,2}
  M::AbstractSparseMatrix
  LUM::SparseArrays.UMFPACK.UmfpackLU{Float64, Int64}
end

function CG1KitePrimalStruct{FT}(::Grids.Quad,backend,Grid) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  Type = Grids.QuadPrimal()
  Type = Grids.Quad()
  DoF = 4
  points = KernelAbstractions.zeros(backend,Float64,DoF,2)
  Comp = 1
  @polyvar x y
  phi = Array{Polynomial,2}(undef,DoF,Comp)
  rotphi = Array{Polynomial,2}(undef,DoF,2)
  xP, w = gaussradau(2)
  xP .= -xP
  lx0 = (x - xP[2])/(xP[1] - xP[2])
  lx1 = (x - xP[1])/(xP[2] - xP[1])
  ly0 = (y - xP[2])/(xP[1] - xP[2])
  ly1 = (y - xP[1])/(xP[2] - xP[1])
  phi[1,1] = lx0 * ly0
  phi[2,1] = lx1 * ly0
  phi[3,1] = lx0 * ly1
  phi[4,1] = lx1 * ly1
  Gradphi = Array{Polynomial,3}(undef,DoF,1,2)
  for i = 1 : DoF
    Gradphi[i,1,1] = differentiate(phi[i,1],x)
    Gradphi[i,1,2] = differentiate(phi[i,1],y)
    rotphi[i,1] = differentiate(phi[i,1],y)
    rotphi[i,2] = -differentiate(phi[i,1],x)
  end  
  iDoF = 1
  @inbounds for j = 1 : 2
    @inbounds for i = 1 : 2
      points[iDoF,1] = xP[i]
      points[iDoF,2] = xP[j]
      iDoF += 1
    end
  end
  
  NumNodes = Grid.NumNodes
  NumEdges = Grid.NumEdges
  NumFaces = Grid.NumFaces
  Nodes = Grid.Nodes
  Edges = Grid.Edges
  Faces = Grid.Faces


  NumNodesD = 0
  NumNodesP = 0
  NodesP = zeros(Int,NumNodes)
  for iN = 1 : NumNodes
    if Nodes[iN].Type == 'F'
      NumNodesP += 1  
      NodesP[iN] = NumNodesP
    end
  end  

  NumEdgesD = 0 
  NumEdgesP = 0 
  EdgesP = zeros(Int,NumEdges)
  for iE = 1 : NumEdges
    if Edges[iE].Type == "E1" || Edges[iE].Type == "E2"  
      NumEdgesD += 1  
    else  
      NumEdgesP += 1  
      EdgesP[iE] = NumEdgesP
    end
  end  
        
  NumG = 0
  Glob = KernelAbstractions.zeros(backend,Int,DoF,NumFaces)
  GlobCPU = zeros(Int,DoF,NumFaces)
  for iF = 1 : NumFaces
    iN = NodesP[Faces[iF].N[3]]  
    GlobCPU[1,iF] = iN  
    iE = EdgesP[Faces[iF].E[3]]
    GlobCPU[2,iF] = NumNodesP + iE
    iE = EdgesP[Faces[iF].E[2]]
    GlobCPU[3,iF] = NumNodesP + iE
    GlobCPU[4,iF] = NumNodesP + NumEdgesP + iF
  end    
  NumG = NumNodesP + NumEdgesP + NumFaces
  NumI = NumG
  copyto!(Glob,GlobCPU)
  M = sparse([1],[1],[1.0])
  LUM = lu(M)
  return CG1KitePrimalStruct{FT,
                  typeof(Glob)}( 
    Glob,
    DoF,
    Comp,
    phi,                      
    rotphi,                      
    Gradphi,                      
    NumG,
    NumI,
    Type,
    points,
    M,
    LUM,
      )
end

mutable struct CG1KiteDualStruct{FT<:AbstractFloat,
                        IT2<:AbstractArray} <: ScalarKitePElement
  Glob::IT2
  DoF::Int
  Comp::Int                      
  phi::Array{Polynomial,2}                       
  rotphi::Array{Polynomial,2}                       
  Gradphi::Array{Polynomial,3}                       
  NumG::Int
  NumI::Int
  Type::Grids.ElementType
  points::Array{Float64,2}
  M::AbstractSparseMatrix
  LUM::SparseArrays.UMFPACK.UmfpackLU{Float64, Int64}
end

function CG1KiteDualStruct{FT}(::Grids.Quad,backend,Grid) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  Type = Grids.QuadPrimal()
  Type = Grids.Quad()
  DoF = 4
  points = KernelAbstractions.zeros(backend,Float64,DoF,2)
  Comp = 1
  @polyvar x y
  phi = Array{Polynomial,2}(undef,DoF,Comp)
  Gradphi = Array{Polynomial,3}(undef,DoF,1,2)
  rotphi = Array{Polynomial,2}(undef,DoF,2)
  xP, w = gaussradau(2)
  lx0 = (x - xP[2])/(xP[1] - xP[2])
  lx1 = (x - xP[1])/(xP[2] - xP[1])
  ly0 = (y - xP[2])/(xP[1] - xP[2])
  ly1 = (y - xP[1])/(xP[2] - xP[1])
  phi[1,1] = lx0 * ly0
  phi[2,1] = lx1 * ly0
  phi[3,1] = lx0 * ly1
  phi[4,1] = lx1 * ly1
  for i = 1 : DoF
    Gradphi[i,1,1] = differentiate(phi[i,1],x)
    Gradphi[i,1,2] = differentiate(phi[i,1],y)
    rotphi[i,1] = differentiate(phi[i,1],y)
    rotphi[i,2] = -differentiate(phi[i,1],x)
  end  
  iDoF = 1
  @inbounds for j = 1 : 2
    @inbounds for i = 1 : 2
      points[iDoF,1] = xP[i]
      points[iDoF,2] = xP[j]
      iDoF += 1
    end
  end
  
  NumFaces = Grid.NumFaces
  NumEdges = Grid.NumEdges
  NumNodes = Grid.NumNodes
  Faces = Grid.Faces
  Edges = Grid.Edges
  Nodes = Grid.Nodes

  NumNodesD = 0
  NumNodesP = 0
  NodesD = zeros(Int,NumNodes)
  for iN = 1 : NumNodes
    if Nodes[iN].Type == 'N'
      NumNodesD += 1
      NodesD[iN] = NumNodesD
    end
  end

  NumEdgesD = 0
  NumEdgesP = 0
  EdgesD = zeros(Int,NumEdges)
  for iE = 1 : NumEdges
    if Edges[iE].Type == "E1" || Edges[iE].Type == "E2"
      NumEdgesD += 1
      EdgesD[iE] = NumEdgesD
    else
      NumEdgesP += 1
    end
  end

  Glob = KernelAbstractions.zeros(backend,Int,DoF,NumFaces)
  GlobCPU = zeros(Int,DoF,NumFaces)
  for iF = 1 : NumFaces
    iN = NodesD[Faces[iF].N[1]]
    GlobCPU[1,iF] = iN
    iE = EdgesD[Faces[iF].E[1]]
    GlobCPU[2,iF] = NumNodesD + iE
    iE = EdgesD[Faces[iF].E[4]]
    GlobCPU[3,iF] = NumNodesD + iE
    GlobCPU[4,iF] = NumNodesD + NumEdgesD + iF
  end
  NumG = NumNodesD + NumEdgesD + NumFaces
  NumI = NumG
  copyto!(Glob,GlobCPU)
  M = sparse([1],[1],[1.0])
  LUM = lu(M)
  return CG1KiteDualStruct{FT,
                  typeof(Glob)}( 
    Glob,
    DoF,
    Comp,
    phi,                      
    rotphi,                      
    Gradphi,                      
    NumG,
    NumI,
    Type,
    points,
    M,
    LUM,
      )
end

mutable struct CGKitePrimalStruct{FT<:AbstractFloat,
                        IT2<:AbstractArray} <: ScalarKitePElement
  Order::Int                      
  Glob::IT2
  DoF::Int
  DoFE::Int
  DoFN::Int
  Comp::Int
  phi::Array{Polynomial,2}
  rotphi::Array{Polynomial,2}
  Gradphi::Array{Polynomial,3}
  NumG::Int
  NumI::Int
  Type::Grids.ElementType
  points::Array{Float64,2}
  NumNodesP::Int
  NodesP::Array{Int,1}
  EdgesP::Array{Int,1}
  M::AbstractSparseMatrix
  LUM
# LUM::SparseArrays.UMFPACK.UmfpackLU{Float64, Int64}
end

function CGKitePrimalStruct{FT}(k,::Grids.Quad,backend,Grid) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  Type = Grids.QuadPrimal()
  DoF = (k + 1) * (k + 1)
  DoFF = k * k
  DoFE = k
  DoFN = 1
  points = KernelAbstractions.zeros(backend,Float64,DoF,2)
  Comp = 1
  @polyvar x[1:2]
  phi = Array{Polynomial,2}(undef,DoF,Comp)
  Gradphi = Array{Polynomial,3}(undef,DoF,1,2)
  rotphi = Array{Polynomial,2}(undef,DoF,2)
  xP, w = gaussradau(k+1)
  xP .= -xP
  L1 = Array{Polynomial,1}(undef,k + 1)
  L2 = Array{Polynomial,1}(undef,k + 1)
  @inbounds for i = 1 : k + 1
    L1[i] = Lagrange(x[1],xP,i)
    L2[i] = Lagrange(x[2],xP,i)
  end
  iDoF = 1
  phi[iDoF,1] = L1[1] * L2[1]
  points[iDoF,1] = xP[1]
  points[iDoF,2] = xP[1]
  iDoF += 1
  for i = 2 : k + 1
    phi[iDoF,1] = L1[i] * L2[1]
    points[iDoF,1] = xP[i]
    points[iDoF,2] = xP[1]
    iDoF += 1
  end  
  for j = 2 : k + 1
    phi[iDoF,1] = L1[1] * L2[j]
    points[iDoF,1] = xP[1]
    points[iDoF,2] = xP[j]
    iDoF += 1
  end  
  for j = 2 : k + 1
    for i = 2 : k + 1
      phi[iDoF,1] = L1[i] * L2[j]
      points[iDoF,1] = xP[i]
      points[iDoF,2] = xP[j]
      iDoF += 1
    end
  end  
      
  for iDoF = 1 : DoF
    Gradphi[iDoF,1,1] = differentiate(phi[iDoF,1],x[1])
    Gradphi[iDoF,1,2] = differentiate(phi[iDoF,1],x[2])
    rotphi[iDoF,1] = differentiate(phi[iDoF,1],x[1])
    rotphi[iDoF,2] = -differentiate(phi[iDoF,1],x[2])
  end  
  
  NumFaces = Grid.NumFaces
  NumEdges = Grid.NumEdges
  NumNodes = Grid.NumNodes
  Faces = Grid.Faces
  Edges = Grid.Edges
  Nodes = Grid.Nodes

  NumNodesD = 0
  NumNodesP = 0
  NodesP = zeros(Int,NumNodes)
  for iN = 1 : NumNodes
    if Nodes[iN].Type == 'F'
      NumNodesP += 1
      NodesP[iN] = NumNodesP
    end
  end

  NumEdgesP = 0
  EdgesP = zeros(Int,NumEdges)
  for iE = 1 : NumEdges
    if Edges[iE].Type == "EM" 
      NumEdgesP += 1
      EdgesP[iE] = NumEdgesP
    end
  end

  Glob = KernelAbstractions.zeros(backend,Int,DoF,NumFaces)
  GlobCPU = zeros(Int,DoF,NumFaces)
  for iF = 1 : NumFaces
    ii = 1
    iN = NodesP[Faces[iF].N[3]]
    GlobCPU[ii,iF] = iN
    ii += 1
    iE = EdgesP[Faces[iF].E[3]]
    for i = 1 : k 
      GlobCPU[ii,iF] = NumNodesP + (iE -1) * DoFE + i
      ii += 1
    end  
    iE = EdgesP[Faces[iF].E[2]]
    for i = 1 : k 
      GlobCPU[ii,iF] = NumNodesP + (iE -1) * DoFE + i
      ii += 1
    end  
    jj = 1
    for j = 1 : k 
      for j = 1 : k 
        GlobCPU[ii,iF] = NumNodesP + NumEdgesP * DoFE + (iF-1) * DoFF + jj
        ii += 1
        jj += 1
      end  
    end  
  end
  NumG = NumNodesP + DoFE * NumEdgesP + DoFF * NumFaces
  NumI = NumG
  copyto!(Glob,GlobCPU)
  M = sparse([1],[1],[1.0])
  LUM = lu(M)
  return CGKitePrimalStruct{FT,
                  typeof(Glob)}( 
    k,              
    Glob,
    DoF,
    DoFE,
    DoFN,
    Comp,
    phi,                      
    rotphi,                      
    Gradphi,                      
    NumG,
    NumI,
    Type,
    points,
    NumNodesP,
    NodesP,
    EdgesP[1:Grid.NumEdgesB],
    M,
    LUM,
      )
end

mutable struct CGKiteDualStruct{FT<:AbstractFloat,
                        IT2<:AbstractArray} <: ScalarKitePElement
  Order::Int                      
  Glob::IT2
  DoF::Int
  Comp::Int
  phi::Array{Polynomial,2}
  rotphi::Array{Polynomial,2}
  Gradphi::Array{Polynomial,3}
  NumG::Int
  NumI::Int
  Type::Grids.ElementType
  points::Array{Float64,2}
  M::AbstractSparseMatrix
  LUM::SparseArrays.UMFPACK.UmfpackLU{Float64, Int64}
end

function CGKiteDualStruct{FT}(k,::Grids.Quad,backend,Grid) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  Type = Grids.QuadPrimal()
  Type = Grids.Quad()
  DoF = (k + 1) * (k + 1)
  points = KernelAbstractions.zeros(backend,Float64,DoF,2)
  Comp = 1
  @polyvar x[1:2]
  phi = Array{Polynomial,2}(undef,DoF,Comp)
  Gradphi = Array{Polynomial,3}(undef,DoF,1,2)
  rotphi = Array{Polynomial,2}(undef,DoF,2)
  xP, w = gaussradau(k+1)
  L1 = Array{Polynomial,1}(undef,k + 1)
  L2 = Array{Polynomial,1}(undef,k + 1)
  @inbounds for i = 1 : k + 1
    L1[i] = Lagrange(x[1],xP,i)
    L2[i] = Lagrange(x[2],xP,i)
  end
  iDoF = 1
  @inbounds for j = 1 : k + 1
    @inbounds for i = 1 : k + 1
      phi[iDoF,1] = L1[i] * L2[j]
      points[iDoF,1] = xP[i]
      points[iDoF,2] = xP[j]
      iDoF += 1
    end
  end
  for iDoF = 1 : DoF
    Gradphi[iDoF,1,1] = differentiate(phi[iDoF,1],x[1])
    Gradphi[iDoF,1,2] = differentiate(phi[iDoF,1],x[2])
    rotphi[iDoF,1] = differentiate(phi[iDoF,1],x[1])
    rotphi[iDoF,2] = -differentiate(phi[iDoF,1],x[2])
  end  
  
  NumFaces = Grid.NumFaces
  NumEdges = Grid.NumEdges
  NumNodes = Grid.NumNodes
  Faces = Grid.Faces
  Edges = Grid.Edges
  Nodes = Grid.Nodes

  NumNodesD = 0
  NumNodesP = 0
  NodesD = zeros(Int,NumNodes)
  for iN = 1 : NumNodes
    if Nodes[iN].Type == 'N'
      NumNodesD += 1
      NodesD[iN] = NumNodesD
    end
  end

  NumEdgesD = 0
  NumEdgesP = 0
  EdgesD = zeros(Int,NumEdges)
  for iE = 1 : NumEdges
    if Edges[iE].Type == "E1" || Edges[iE].Type == "E2"
      NumEdgesD += 1
      EdgesD[iE] = NumEdgesD
    else
      NumEdgesP += 1
    end
  end

  Glob = KernelAbstractions.zeros(backend,Int,DoF,NumFaces)
  GlobCPU = zeros(Int,DoF,NumFaces)
  for iF = 1 : NumFaces
    ii = 1
    iN = NodesD[Faces[iF].N[1]]
    GlobCPU[ii,iF] = iN
    ii += 1
    iE = EdgesD[Faces[iF].E[1]]
    for i = 2 : k + 1
      GlobCPU[ii,iF] = NumNodesD + k * iE - (k + 1 - i)
      ii += 1
    end  
    iE = EdgesD[Faces[iF].E[4]]
    jj = k * k
    for j = 2 : k + 1
      GlobCPU[ii,iF] = NumNodesD + k * iE - (k + 1 - j)
      ii += 1
      for i = 2 : k + 1
        jj -= 1  
        GlobCPU[ii,iF] = NumNodesD + k * NumEdgesD + k * k * iF - jj
        ii += 1
      end    
    end  
  end
  NumG = NumNodesD + k * NumEdgesD + k * k * NumFaces
  NumI = NumG
  copyto!(Glob,GlobCPU)
  M = sparse([1],[1],[1.0])
  LUM = lu(M)
  return CGKiteDualStruct{FT,
                  typeof(Glob)}( 
    k,
    Glob,
    DoF,
    Comp,
    phi,                      
    rotphi,                      
    Gradphi,                      
    NumG,
    NumI,
    Type,
    points,
    M,
    LUM,
      )
end

mutable struct CG1KiteDualHDiv{FT<:AbstractFloat,
                        IT2<:AbstractArray} <: HDivKiteDElement
  Glob::IT2
  DoF::Int
  Comp::Int                      
  phi::Array{Polynomial,2}                       
  Divphi::Array{Polynomial,2}                       
  NumG::Int
  NumI::Int
  Type::Grids.ElementType
  points::Array{Float64,2}
  M::AbstractSparseMatrix
  LUM::SparseArrays.UMFPACK.UmfpackLU{Float64, Int64}
end

function CG1KiteDualHDiv{FT}(::Grids.Quad,backend,Grid) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  Type = Grids.QuadDual()
  Type = Grids.Quad()
  DoF = 8
  Comp = 2
  points = KernelAbstractions.zeros(backend,Float64,4,2)
  @polyvar x y
  phi = Array{Polynomial,2}(undef,DoF,Comp)
  Divphi = Array{Polynomial,2}(undef,DoF,1)
  xP, w = gaussradau(2)
  iDoF = 1
  @inbounds for j = 1 : 2
    @inbounds for i = 1 : 2
      points[iDoF,1] = xP[i]
      points[iDoF,2] = xP[j]
      iDoF += 1
    end
  end
  lx0 = (x - xP[2])/(xP[1] - xP[2])
  lx1 = (x - xP[1])/(xP[2] - xP[1])
  ly0 = (y - xP[2])/(xP[1] - xP[2])
  ly1 = (y - xP[1])/(xP[2] - xP[1])
  p0 = 0.0*x + 0.0*y
# x line 1 --> 2, y0
  phi[1,1] = p0
  phi[1,2] = lx0 * ly0
  phi[2,1] = p0
  phi[2,2] = lx1 * ly0
  
# y line 1 --> 4, x0
  phi[3,1] = -lx0 * ly0
  phi[3,2] = p0
  phi[4,1] = -lx0 * ly1
  phi[4,2] = p0

# x line 3 --> 4, y1
  phi[5,1] = p0
  phi[5,2] = lx0 * ly1
  phi[6,1] = p0
  phi[6,2] = lx1 * ly1

# y line 2 --> 3, x1
  phi[7,1] = -lx1 * ly0
  phi[7,2] = p0
  phi[8,1] = -lx1 * ly1
  phi[8,2] = p0

  @. phi = -phi # Oswald

  Divphi[1,1] =  (differentiate(phi[1,1],x) + differentiate(phi[1,2],y))
  Divphi[2,1] =  (differentiate(phi[2,1],x) + differentiate(phi[2,2],y))
  Divphi[3,1] =  (differentiate(phi[3,1],x) + differentiate(phi[3,2],y))
  Divphi[4,1] =  (differentiate(phi[4,1],x) + differentiate(phi[4,2],y))
  Divphi[5,1] =  (differentiate(phi[5,1],x) + differentiate(phi[5,2],y))
  Divphi[6,1] =  (differentiate(phi[6,1],x) + differentiate(phi[6,2],y))
  Divphi[7,1] =  (differentiate(phi[7,1],x) + differentiate(phi[7,2],y))
  Divphi[8,1] =  (differentiate(phi[8,1],x) + differentiate(phi[8,2],y))

  NumFaces = Grid.NumFaces
  NumEdges = Grid.NumEdges
  NumNodes = Grid.NumNodes
  Faces = Grid.Faces
  Edges = Grid.Edges
  Nodes = Grid.Nodes

  NumEdgesD = 0
  EdgesD = zeros(Int,NumEdges)
  for iE = 1 : NumEdges
    if Edges[iE].Type == "E1" || Edges[iE].Type == "E2"
      NumEdgesD += 1
      EdgesD[iE] = NumEdgesD
    end
  end

  Glob = KernelAbstractions.zeros(backend,Int,DoF,NumFaces)
  GlobCPU = zeros(Int,DoF,NumFaces)
  for iF = 1 : NumFaces
    iE = EdgesD[Faces[iF].E[1]]
    GlobCPU[1,iF] = 2 * iE - 1
    GlobCPU[2,iF] = 2 * iE
    iE = EdgesD[Faces[iF].E[4]]
    GlobCPU[3,iF] = 2 * iE - 1
    GlobCPU[4,iF] = 2 * iE
    GlobCPU[5,iF] = 2 * NumEdgesD + 4 * iF - 3
    GlobCPU[6,iF] = 2 * NumEdgesD + 4 * iF - 2
    GlobCPU[7,iF] = 2 * NumEdgesD + 4 * iF - 1
    GlobCPU[8,iF] = 2 * NumEdgesD + 4 * iF
  end

  NumG = 2 * NumEdgesD + 4 * NumFaces
  NumI = NumG
# Add Boundary

  copyto!(Glob,GlobCPU)
  M = sparse([1],[1],[1.0])
  LUM = lu(M)
  return CG1KiteDualHDiv{FT,
                  typeof(Glob)}( 
    Glob,
    DoF,
    Comp,
    phi,                      
    Divphi,                      
    NumG,
    NumI,
    Type,
    points,
    M,
    LUM,
      )
end

mutable struct CG1KiteDualHCurl{FT<:AbstractFloat,
                        IT1<:AbstractArray,
                        IT2<:AbstractArray} <: HCurlKiteDElement
  Glob::IT2
  DoF::Int
  Comp::Int                      
  phi::Array{Polynomial,2}                       
  Curlphi::Array{Polynomial,2}                       
  NumG::Int
  NumI::Int
  ListB::IT1
  Type::Grids.ElementType
  M::AbstractSparseMatrix
  LUM::SparseArrays.UMFPACK.UmfpackLU{Float64, Int64}
end

function CG1KiteDualHCurl{FT}(::Grids.Quad,backend,Grid) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  Type = Grids.QuadDual()
  DoF = 8
  Comp = 2
  @polyvar x y
  phi = Array{Polynomial,2}(undef,DoF,Comp)
  Curlphi = Array{Polynomial,2}(undef,DoF,1)
  xP, w = gaussradau(2)
  xP .= -xP
  lx0 = (x - xP[2])/(xP[1] - xP[2])
  lx1 = (x - xP[1])/(xP[2] - xP[1])
  ly0 = (y - xP[2])/(xP[1] - xP[2])
  ly1 = (y - xP[1])/(xP[2] - xP[1])
# x line 3 --> 4, y0
  p0 = 0.0*x + 0.0*y
  phi[1,1] = lx0 * ly0
  phi[1,2] = p0
  phi[2,1] = lx1 * ly0
  phi[2,2] = p0
  
# y line 2 --> 3, x0
  p0 = 0.0*x + 0.0*y
  phi[3,1] = p0
  phi[3,2] = lx0 * ly0
  phi[4,1] = p0
  phi[4,2] = lx0 * ly1

# x line 3 --> 4, y1
  phi[5,1] = -lx0 * ly1
  phi[5,2] = p0
  phi[6,1] = -lx1 * ly1
  phi[6,2] = p0

# y line 2 --> 3, x1
  phi[7,1] = p0
  phi[7,2] = lx1 * ly0
  phi[8,1] = p0
  phi[8,2] = lx1 * ly1

  for i = 1 : DoF
    Curlphi[i,1] = -differentiate(phi[i,1],y) + differentiate(phi[i,2],x)
  end

  NumFaces = Grid.NumFaces
  NumEdges = Grid.NumEdges
  NumNodes = Grid.NumNodes
  Faces = Grid.Faces
  Edges = Grid.Edges
  Nodes = Grid.Nodes

# Kite list  
  KiteList = Dict()
  iKite = 1
  for iN = 1 : NumNodes
    if Nodes[iN].Type == 'F' 
      NumF = length(Nodes[iN].F)
      Offset = 1 + 2 * NumF
      for i = 1 : NumF
        iF = Nodes[iN].F[i]  
        N1 = Faces[iF].N[1]
        N3 = Faces[iF].N[3]  
        KiteList[(N1,N3)] = iKite
        iKite +=1 
      end 
    end 
  end 



  Glob = KernelAbstractions.zeros(backend,Int,DoF,NumFaces)
  GlobCPU = zeros(Int,DoF,NumFaces)
  iOff = 0
  ListBCPU = Int64[]
  for iN = 1 : NumNodes
    if Nodes[iN].Type == 'N'  
      NumF = length(Nodes[iN].F)  
      OffsetE = 2 * NumF
      for i = 1 : NumF
        iF = Nodes[iN].F[i]
        N1 = Faces[iF].N[1]
        N3 = Faces[iF].N[3]
        iKite = KiteList[(N1,N3)]
        #Edge e1
        GlobCPU[1,iKite] = iOff + 2 * i -1 
        GlobCPU[2,iKite] = iOff + 2 * i
        #Edge e2
        if i < NumF
          GlobCPU[3,iKite] = iOff + 2 * (i + 1) - 1
          GlobCPU[4,iKite] = iOff + 2 * (i + 1)
        else  
          GlobCPU[3,iKite] = iOff + 2 * 1 - 1
          GlobCPU[4,iKite] = iOff + 2 * 1
        end    
        # Interior e1
        GlobCPU[5,iKite] = iOff + OffsetE + 4 * i - 3
        GlobCPU[6,iKite] = iOff + OffsetE + 4 * i - 2
        # Interior e2
        GlobCPU[7,iKite] = iOff + OffsetE + 4 * i - 1
        GlobCPU[8,iKite] = iOff + OffsetE + 4 * i 
      end 
      iOff += 6 * NumF
    elseif Nodes[iN].Type == 'B' || Nodes[iN].Type == 'P'
      NumF = length(Nodes[iN].F)
      OffsetE = 2 * (NumF + 1)
      for i = 1 : NumF
        iF = Nodes[iN].F[i]
        N1 = Faces[iF].N[1]
        N3 = Faces[iF].N[3]
        iKite = KiteList[(N1,N3)]
        #Edge e1
        GlobCPU[1,iKite] = iOff + 2 * i -1
        GlobCPU[2,iKite] = iOff + 2 * i
        if i == 1 
          push!(ListBCPU,iOff + 2 * i -1)
          push!(ListBCPU,iOff + 2 * i)
        end  
        #Edge e2
        GlobCPU[3,iKite] = iOff + 2 * (i + 1) - 1
        GlobCPU[4,iKite] = iOff + 2 * (i + 1)
        if i == NumF
          push!(ListBCPU,iOff + 2 * (i + 1) - 1)
          push!(ListBCPU,iOff + 2 * (i + 1))
        end  
            
        # Interior e1
        GlobCPU[5,iKite] = iOff + OffsetE + 4 * i - 3
        GlobCPU[6,iKite] = iOff + OffsetE + 4 * i - 2
        # Interior e2
        GlobCPU[7,iKite] = iOff + OffsetE + 4 * i - 1
        GlobCPU[8,iKite] = iOff + OffsetE + 4 * i
      end
      iOff += 4 * NumF + 2 * (NumF + 1)
    end 
  end
  NumG = iOff
  NumI = NumG
  copyto!(Glob,GlobCPU)
  M = sparse([1],[1],[1.0])
  LUM = lu(M)
  ListB = KernelAbstractions.zeros(backend,Int,size(ListBCPU))
  copyto!(ListB,ListBCPU)
  return CG1KiteDualHCurl{FT,
                  typeof(ListB),
                  typeof(Glob)}( 
    Glob,
    DoF,
    Comp,
    phi,                      
    Curlphi,                      
    NumG,
    NumI,
    ListB,
    Type,
    M,
    LUM,
      )
end

