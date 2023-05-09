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

mutable struct GridStruct
    nz::Int
    zP::Array{Float64, 1}
    z::Array{Float64, 1}
    dzeta::Array{Float64, 1}
    H::Float64
    NumFaces::Int
    NumGhostFaces::Int
    Faces::Array{Face, 1}
    NumEdges::Int
    Edges::Array{Edge, 1}
    NumNodes::Int
    Nodes::Array{Node, 1}
    Form::String
    Type ::String
    Dim::Int
    Rad::Float64
    NumEdgesI::Int
    NumEdgesB::Int
    nBar3::Array{Float64, 2}
    nBar::Array{Float64, 2}
    Topography::NamedTuple
    colors::Array{Array{Int, 1}, 1}
    Spline_2d::Dierckx.Spline2D
    BoundaryFaces::Array{Int,1}
    InteriorFaces::Array{Int,1}
end
function InitGrid(nz,Topography)
  zP=zeros(nz)
  z=zeros(nz+1)
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
  Type=""
  Dim=0
  Rad=0.0
  NumEdgesI=0
  NumEdgesB=0
  nBar3=zeros(0,0)
  nBar=zeros(0,0)
  colors=[[]]
  Spline_2d = Spline2D(zeros(0),zeros(0),zeros(0),0,0,0.0)
  BoundaryFaces = zeros(Int,0)
  InteriorFaces = zeros(Int,0)
   return GridStruct(
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
    Topography,
    colors,
    Spline_2d,
    BoundaryFaces,
    InteriorFaces,
    )
end   
    
