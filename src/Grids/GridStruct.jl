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
  AdaptGrid::Any

end
    
