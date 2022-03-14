function GridPolygon=PolygonToGrid(Polygon,Grid,OrientFace)
GridPolygon.nBar=[ 0  1   0   -1  
                   1  0  -1    0];
GridPolygon.Dim=3;
GridPolygon.Type='Quad';
NumNodesFace=size(Polygon.N,2);
NumNodes=2*NumNodesFace+1;
Nodes(1:NumNodes)=Node([0 0 0],0);
NodeNumber=1;
BoundaryNodeNumber=1;
BoundaryNode=zeros(1,2*NumNodesFace);
%BoundaryNode=zeros(1,1);
MidPoint=Grid.Nodes(Polygon.N(end)).P;
for i=1:NumNodesFace-1
  Nodes(NodeNumber)=Node(Grid.Nodes(Polygon.N(i)).P',NodeNumber);
  MidPoint=MidPoint+Grid.Nodes(Polygon.N(i)).P;
  BoundaryNode(BoundaryNodeNumber)=NodeNumber;
  BoundaryNodeNumber=BoundaryNodeNumber+1;
  NodeNumber=NodeNumber+1;
  Nodes(NodeNumber)=Node(0.5*(Grid.Nodes(Polygon.N(i)).P'+Grid.Nodes(Polygon.N(i+1)).P'),NodeNumber);
  BoundaryNode(BoundaryNodeNumber)=NodeNumber;
  BoundaryNodeNumber=BoundaryNodeNumber+1;
  NodeNumber=NodeNumber+1;
end
Nodes(NodeNumber)=Node(Grid.Nodes(Polygon.N(NumNodesFace)).P',NodeNumber);
BoundaryNode(BoundaryNodeNumber)=NodeNumber;
BoundaryNodeNumber=BoundaryNodeNumber+1;
NodeNumber=NodeNumber+1;
Nodes(NodeNumber)=Node(0.5*(Grid.Nodes(Polygon.N(NumNodesFace)).P'+Grid.Nodes(Polygon.N(1)).P'),NodeNumber);
BoundaryNode(BoundaryNodeNumber)=NodeNumber;
NodeNumber=NodeNumber+1;
MidPoint=(1/NumNodesFace)*MidPoint;
Nodes(NodeNumber)=Node(MidPoint',NodeNumber);

GridPolygon.Nodes=Nodes;
GridPolygon.NumNodes=NumNodes;
GridPolygon.BoundaryNode=BoundaryNode;
InteriorNode(1)=NodeNumber;
GridPolygon.InteriorNode=InteriorNode;


NumEdges=3*NumNodesFace;
Edges(1:NumEdges)=Edge([1 2],GridPolygon,0,0);
EdgeNumber=1;
BoundaryEdgeNumber=1;
BoundaryEdge=zeros(1,2*NumNodesFace);
InteriorEdgeNumber=1;
InteriorEdge=zeros(1,NumNodesFace);
N1=1;
N2=2;
EdgeNumberI=1;
EdgeNumberB=1;
for i=1:NumNodesFace-1
  Edges(EdgeNumber)=Edge([N1 N2],GridPolygon,EdgeNumber,EdgeNumberB+NumNodesFace);
  BoundaryEdge(BoundaryEdgeNumber)=EdgeNumberB+NumNodesFace;
  BoundaryEdgeNumber=BoundaryEdgeNumber+1;
  EdgeNumber=EdgeNumber+1;
  EdgeNumberB=EdgeNumberB+1;
  N1=N1+1;
  N2=N2+1;
  Edges(EdgeNumber)=Edge([N1 N2],GridPolygon,EdgeNumber,EdgeNumberB+NumNodesFace);
  BoundaryEdge(BoundaryEdgeNumber)=EdgeNumberB+NumNodesFace;
  BoundaryEdgeNumber=BoundaryEdgeNumber+1;
  EdgeNumber=EdgeNumber+1;
  EdgeNumberB=EdgeNumberB+1;
  N1=N1+1;
  N2=N2+1;
end
Edges(EdgeNumber)=Edge([N1 N2],GridPolygon,EdgeNumber,EdgeNumberB+NumNodesFace);
BoundaryEdge(BoundaryEdgeNumber)=EdgeNumberB+NumNodesFace;
BoundaryEdgeNumber=BoundaryEdgeNumber+1;
EdgeNumber=EdgeNumber+1;
EdgeNumberB=EdgeNumberB+1;
N1=N1+1;
N2=1;
Edges(EdgeNumber)=Edge([N1 N2],GridPolygon,EdgeNumber,EdgeNumberB+NumNodesFace);
BoundaryEdge(BoundaryEdgeNumber)=EdgeNumberB+NumNodesFace;
EdgeNumber=EdgeNumber+1;
EdgeNumberB=EdgeNumberB+1;

N1=2;
N2=2*NumNodesFace+1;
for i=1:NumNodesFace
  Edges(EdgeNumber)=Edge([N1 N2],GridPolygon,EdgeNumber,EdgeNumberI);
  InteriorEdge(InteriorEdgeNumber)=EdgeNumberI;
  InteriorEdgeNumber=InteriorEdgeNumber+1;
  EdgeNumber=EdgeNumber+1;
  EdgeNumberI=EdgeNumberI+1;
  N1=N1+2;
end
GridPolygon.Edges=Edges;
GridPolygon.NumEdges=NumEdges;
GridPolygon.BoundaryEdge=BoundaryEdge;
GridPolygon.InteriorEdge=InteriorEdge;
GridPolygon.NumEdgesI=EdgeNumberI-1;
GridPolygon.NumEdgesB=EdgeNumberB-1;


NumFaces=NumNodesFace;
Faces(1:NumFaces)=Face([0 0 0 0],GridPolygon,0,'x');
E1=2*NumNodesFace;
E2=3*NumNodesFace;
E3=2*NumNodesFace+1;
E4=1;
FaceNumber=1;
Type='o';
[Faces(FaceNumber),GridPolygon]=Face([E1 E2 E3 E4],GridPolygon,FaceNumber,Type,OrientFace);
FaceNumber=FaceNumber+1;
E1=2;
E2=2*NumNodesFace+1;
E3=2*NumNodesFace+2;
E4=3;
for i=2:NumNodesFace
  [Faces(FaceNumber),GridPolygon]=Face([E1 E2 E3 E4],GridPolygon,FaceNumber,Type,OrientFace);
  FaceNumber=FaceNumber+1;
  E1=E1+2;
  E2=E2+1;
  E3=E3+1;
  E4=E4+2;
end
GridPolygon.Faces=Faces;
GridPolygon.NumFaces=NumFaces;
GridPolygon=Orientation(GridPolygon);
GridPolygon=Renumbering(GridPolygon); 
end