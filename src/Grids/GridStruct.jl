mutable struct BoundaryStruct
 WE::String
 SN::String
 BT::String
end
function Boundary()
  WE = ""
  SN = ""
  BT = ""
  return BoundaryStruct(
    WE,
    SN,
    BT
  )
end  

mutable struct GridStruct{FT<:AbstractFloat,
                          AT1<:AbstractArray}
  nz::Int
  zP::Array{FT, 1}
  z::AT1
  dzeta::Array{FT, 1}
  H::FT
  NumFaces::Int
  NumGhostFaces::Int
  Faces::Array{Face, 1}
  NumEdges::Int
  Edges::Array{Edge, 1}
  NumNodes::Int
  Nodes::Array{Node, 1}
  Form::String
  Type ::ElementType
  Dim::Int
  Rad::FT
  NumEdgesI::Int
  NumEdgesB::Int
  nBar3::Array{FT, 2}
  nBar::Array{FT, 2}
  colors::Array{Array{Int, 1}, 1}
  NumBoundaryFaces::Int
end
function GridStruct{FT}(backend,nz) where FT <: AbstractFloat
  zP=zeros(nz)
  z=KernelAbstractions.zeros(backend,FT,nz+1)
  dzeta=zeros(nz)
  H=0.0
  NumFaces=0
  NumGhostFaces=0
  Faces=Array{Face}(undef, 0)
  NumEdges=0
  Edges=Array{Edge}(undef, 0)
  NumNodes=0
  Nodes=Array{Node}(undef, 0)
  Form=""
  Type=nothing
  Dim=0
  Rad=0.0
  NumEdgesI=0
  NumEdgesB=0
  nBar3=zeros(0,0)
  nBar=zeros(0,0)
  colors=[[]]
  NumBoundaryFaces = 0
   return GridStruct{FT,
                     typeof(z)}(
    nz,
    zP,
    z,
    dzeta,
    H,
    NumFaces,
    NumGhostFaces,
    Faces,
    NumEdges,
    Edges, 
    NumNodes,
    Nodes, 
    Form,
    Type,
    Dim,
    Rad,
    NumEdgesI,
    NumEdgesB,
    nBar3,
    nBar,
    colors,
    NumBoundaryFaces,
    )
end   
    
