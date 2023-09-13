Rad = 1.0

#struct Point{FT}
#  x::FT
#  y::FT
#  z::FT
#end 

function Norm(P::Point)
  sqrt(P.x * P.x + P.y * P.y + P.z * P.z)
end  

function Div(P::Point,s)
  return Point(P.x / s, P.y / s, P.z / s)
end  

function MidPoint(P1::Point,P2::Point)
  P = Point(0.5 * (P1.x + P2.x),
    0.5 * (P1.y + P2.y), 0.5 * (P1.z + P2.z))
  M = Div(P, Norm(P) / Rad)
end  

mutable struct NodeTri_T
  P::Point
  Edge::Array{Int,1}
  Number::Int
end

function NodeTri_T(P,Edge)
  NodeTri_T(P,Edge,0)
end  

mutable struct EdgeTri_T{TN}
  Node1::ListNode{TN}
  Node2::ListNode{TN}
  Face::Array{Int,1}
  Number::Int
end

function EdgeTri_T(Node1,Node2,Face)
  EdgeTri_T(Node1,Node2,Face,0)
end  

mutable struct FaceTri_T{TE}
  Edge1::ListNode{TE}
  Edge2::ListNode{TE}
  Edge3::ListNode{TE}
  OrientE1::Int
  OrientE2::Int
  OrientE3::Int
  Number::Int
end

function FaceTri_T(Edge1,Edge2,Edge3,OrientE1,OrientE2,OrientE3)
  FaceTri_T(Edge1,Edge2,Edge3,OrientE1,OrientE2,OrientE3,0)
end  
  
mutable struct TriangularGrid_T{TN,TE,TF}
  NodeList::DoublyLinkedList{TN}
  EdgeList::DoublyLinkedList{TE}
  FaceList::DoublyLinkedList{TF}
end

function TriangularGrid_T()
  NodeList = DoublyLinkedList{NodeTri_T}()
  EdgeList = DoublyLinkedList{EdgeTri_T}()
  FaceList = DoublyLinkedList{FaceTri_T}()
  return TriangularGrid_T{NodeTri_T,EdgeTri_T,FaceTri_T}(
    NodeList,
    EdgeList,
    FaceList,
    )
end  

function CreateIcosahedronGrid()
  Rad = 1.0
  IcosahedronGrid = TriangularGrid_T()

# Nodes
  Node = NodeTri_T(Point(0.0,0.0,Rad),[0])
  push!(IcosahedronGrid.NodeList, Node)
  phi = atan(0.5)
  lam = 0.0
  for i = 1 : 5
    Node = NodeTri_T(Point(Rad*cos(lam)*cos(phi),Rad*sin(lam)*cos(phi),Rad*sin(phi)),[0])  
    push!(IcosahedronGrid.NodeList, Node)
    lam = lam + 2.0 * pi / 5.0
  end
  phi = -atan(0.5)
  lam = pi / 5.0
  for i = 1 : 5
    Node = NodeTri_T(Point(Rad*cos(lam)*cos(phi),Rad*sin(lam)*cos(phi),Rad*sin(phi)),[0])  
    push!(IcosahedronGrid.NodeList, Node)
    lam = lam + 2.0 * pi / 5.0
  end
  Node = NodeTri_T(Point(0.0,0.0,-Rad),[0])
  push!(IcosahedronGrid.NodeList, Node)

  NodeTop = getnode(IcosahedronGrid.NodeList,1)
  NodeLayerTopFirst = getnode(IcosahedronGrid.NodeList,2)
  NodeLayerBottomFirst = getnode(IcosahedronGrid.NodeList,7)
  NodeBottom = getnode(IcosahedronGrid.NodeList,12)

# Edges
  NodeLayerTop = NodeLayerTopFirst
  for i = 1 : 5
    Edge = EdgeTri_T(NodeLayerTop,NodeTop,[0])
    push!(IcosahedronGrid.EdgeList, Edge)
    NodeLayerTop = NodeLayerTop.next
  end  
  NodeLayerTop = NodeLayerTopFirst
  for i = 1 : 4
    Edge = EdgeTri_T(NodeLayerTop,NodeLayerTop.next,[0])
    push!(IcosahedronGrid.EdgeList, Edge)
    NodeLayerTop = NodeLayerTop.next
  end  
  Edge = EdgeTri_T(NodeLayerTop,NodeLayerTopFirst,[0])
  push!(IcosahedronGrid.EdgeList, Edge)

  NodeLayerTop = NodeLayerTopFirst
  NodeLayerBottom = NodeLayerBottomFirst
  for i = 1 : 4
    Edge = EdgeTri_T(NodeLayerBottom,NodeLayerTop,[0])
    push!(IcosahedronGrid.EdgeList, Edge)
    Edge = EdgeTri_T(NodeLayerBottom,NodeLayerTop.next,[0])
    push!(IcosahedronGrid.EdgeList, Edge)
    NodeLayerTop = NodeLayerTop.next
    NodeLayerBottom = NodeLayerBottom.next
  end  
  Edge = EdgeTri_T(NodeLayerBottom,NodeLayerTop,[0])
  push!(IcosahedronGrid.EdgeList, Edge)
  Edge = EdgeTri_T(NodeLayerBottom,NodeLayerTopFirst,[0])
  push!(IcosahedronGrid.EdgeList, Edge)

  NodeLayerBottom = NodeLayerBottomFirst
  for i = 1 : 4
    Edge = EdgeTri_T(NodeLayerBottom,NodeLayerBottom.next,[0])
    push!(IcosahedronGrid.EdgeList, Edge)
    NodeLayerBottom = NodeLayerBottom.next
  end  
  Edge = EdgeTri_T(NodeLayerBottom,NodeLayerBottomFirst,[0])
  push!(IcosahedronGrid.EdgeList, Edge)

  NodeLayerBottom = NodeLayerBottomFirst
  for i = 1 : 5
    Edge = EdgeTri_T(NodeBottom,NodeLayerBottom,[0])
    push!(IcosahedronGrid.EdgeList, Edge)
    NodeLayerBottom = NodeLayerBottom.next
  end  

  EdgeTopFirst = getnode(IcosahedronGrid.EdgeList,1)
  EdgeLayerTopFirst = getnode(IcosahedronGrid.EdgeList,6)
  EdgeMidFirst = getnode(IcosahedronGrid.EdgeList,11)
  EdgeLayerBottomFirst = getnode(IcosahedronGrid.EdgeList,21)
  EdgeBottomFirst = getnode(IcosahedronGrid.EdgeList,26)

  #Faces
  EdgeTop = EdgeTopFirst
  EdgeLayerTop = EdgeLayerTopFirst
  for i = 1 : 4
    Face = FaceTri_T(EdgeLayerTop,EdgeTop.next,EdgeTop,1,1,-1)
    push!(IcosahedronGrid.FaceList, Face)
    EdgeTop = EdgeTop.next
    EdgeLayerTop = EdgeLayerTop.next
  end  
  Face = FaceTri_T(EdgeLayerTop,EdgeTopFirst,EdgeTop,1,1,-1)
  push!(IcosahedronGrid.FaceList, Face)

  EdgeLayerTop = EdgeLayerTopFirst
  EdgeMid = EdgeMidFirst
  EdgeLayerBottom = EdgeLayerBottomFirst
  for i = 1 : 4
    Face = FaceTri_T(EdgeMid,EdgeMid.next,EdgeLayerTop,-1,1,-1)
    push!(IcosahedronGrid.FaceList, Face)
    EdgeMid = EdgeMid.next
    Face = FaceTri_T(EdgeMid.next,EdgeMid,EdgeLayerBottom,1,-1,1)
    push!(IcosahedronGrid.FaceList, Face)
    EdgeLayerTop = EdgeLayerTop.next
    EdgeMid = EdgeMid.next
    EdgeLayerBottom = EdgeLayerBottom.next
  end  
  Face = FaceTri_T(EdgeMid,EdgeMid.next,EdgeLayerTop,-1,1,-1)
  push!(IcosahedronGrid.FaceList, Face)
  EdgeMid = EdgeMid.next
  Face = FaceTri_T(EdgeMidFirst,EdgeMid,EdgeLayerBottom,1,-1,1)
  push!(IcosahedronGrid.FaceList, Face)

  EdgeLayerBottom = EdgeLayerBottomFirst
  EdgeBottom = EdgeBottomFirst
  for i = 1 : 4
    Face = FaceTri_T(EdgeLayerBottom,EdgeBottom.next,EdgeBottom,1,-1,1)
    push!(IcosahedronGrid.FaceList, Face)
    EdgeLayerBottom = EdgeLayerBottom.next
    EdgeBottom = EdgeBottom.next
  end  
  Face = FaceTri_T(EdgeLayerBottom,EdgeBottomFirst,EdgeBottom,1,-1,1)
  push!(IcosahedronGrid.FaceList, Face)
  return IcosahedronGrid
      
end

function RefineEdge!(Edge,NodeList,EdgeList)
  P1 = Edge.data.Node1.data.P
  P2 = Edge.data.Node2.data.P
  NodeM = newnode(NodeList, NodeTri_T(MidPoint(P1,P2),[0]))
  insertafter!(NodeM, Edge.data.Node1)
  EdgeNew = newnode(EdgeList, EdgeTri_T(NodeM,Edge.data.Node2,[0]))
  insertafter!(EdgeNew, Edge)
  Edge.data.Node2 = NodeM
end  

function RefineEdgeTriangularGrid!(TriangularGrid)

  NodeList = TriangularGrid.NodeList
  EdgeList = TriangularGrid.EdgeList

  Edge = head(TriangularGrid.EdgeList)
  while ~attail(Edge)
    RefineEdge!(Edge,NodeList,EdgeList)
    Edge = Edge.next.next
  end  
end

function RefineFace!(Face,EdgeList,FaceList)

 Edge1 = Face.data.Edge1
 Edge1N = Edge1.next
 Edge2 = Face.data.Edge2
 Edge2N = Edge2.next
 Edge3 = Face.data.Edge3
 Edge3N = Edge3.next
 if Face.data.OrientE1  == 1 
   Face2Edge2 = Edge1
   Face2OrientE2 = 1
   Face3Edge1 = Edge1N
   Face3OrientE1 = 1
 else
   Face2Edge2 = Edge1N
   Face2OrientE2 = -1
   Face3Edge1 = Edge1
   Face3OrientE1 = -1
 end  
 if Face.data.OrientE2 == 1
   Face3Edge2 = Edge2
   Face3OrientE2 = 1
   Face1Edge1 = Edge2N
   Face1OrientE1 = 1
 else
   Face3Edge2 = Edge2N
   Face3OrientE2 = -1
   Face1Edge1 = Edge2
   Face1OrientE1 = -1
 end
 if Face.data.OrientE3 == 1
   Face1Edge2 = Edge3
   Face1OrientE2 = 1
   Face2Edge1 = Edge3N
   Face2OrientE1 = 1
 else
   Face1Edge2 = Edge3N
   Face1OrientE2 = -1
   Face2Edge1 = Edge3
   Face2OrientE1 = -1
 end
 EdgeI1 = newnode(EdgeList, EdgeTri_T(Edge2.data.Node2,Edge3.data.Node2,[0]))
 insertafter!(EdgeI1, Edge1N)
 EdgeI2 = newnode(EdgeList, EdgeTri_T(Edge3.data.Node2,Edge1.data.Node2,[0]))
 insertafter!(EdgeI2, Edge2N)
 EdgeI3 = newnode(EdgeList, EdgeTri_T(Edge1.data.Node2,Edge2.data.Node2,[0]))
 insertafter!(EdgeI3, Edge3N)

 Face1Edge3 = EdgeI1
 Face1OrientE3 = -1
 Face1 = newnode(FaceList,FaceTri_T(Face1Edge1,Face1Edge2,Face1Edge3,Face1OrientE1,
   Face1OrientE2,Face1OrientE3))
 insertafter!(Face1, Face)

 Face2Edge3 = EdgeI2
 Face2OrientE3 = -1
 Face2 = newnode(FaceList,FaceTri_T(Face2Edge1,Face2Edge2,Face2Edge3,Face2OrientE1,
   Face2OrientE2,Face2OrientE3))
 insertafter!(Face2, Face)

 Face3Edge3 = EdgeI3
 Face3OrientE3 = -1
 Face3 = newnode(FaceList,FaceTri_T(Face3Edge1,Face3Edge2,Face3Edge3,Face3OrientE1,
   Face3OrientE2,Face3OrientE3))
 insertafter!(Face3, Face)

 Face.data.Edge1 = EdgeI1
 Face.data.OrientE1 = 1
 Face.data.Edge2 = EdgeI2
 Face.data.OrientE2 = 1
 Face.data.Edge3 = EdgeI3
 Face.data.OrientE3 = 1

end

function RefineFaceTriangularGrid!(TriangularGrid)

  EdgeList = TriangularGrid.EdgeList
  FaceList = TriangularGrid.FaceList

  Face = head(FaceList)
  while ~attail(Face)
    RefineFace!(Face,EdgeList,FaceList)
    Face = Face.next
    Face = Face.next
    Face = Face.next
    Face = Face.next
  end
end


