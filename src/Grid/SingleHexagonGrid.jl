function Grid=HexagonGrid(NumNodes,OrientFace)

Grid.Dim=3;
Grid.Type='Hexa';

Nodes(1:NumNodes)=Node([0 0 0],0);
NodeNumber=1;
for k=1:NumNodes
  xx=cos(2*pi*k/NumNodes+pi/NumNodes);
  yy=sin(2*pi*k/NumNodes+pi/NumNodes);
  Nodes(NodeNumber)=Node([xx;yy;0],NodeNumber);
  NodeNumber=NodeNumber+1;
end
Grid.Nodes=Nodes;
Grid.NumNodes=NumNodes;

NumEdges=NumNodes;
Edges(1:NumEdges)=Edge([1 2],Grid,0,0);
N1=1;
N2=2;
EdgeNumber=1;
for i=1:NumEdges-1
  Edges(EdgeNumber)=Edge([N1 N2],Grid,EdgeNumber,EdgeNumber);
  EdgeNumber=EdgeNumber+1;
  N1=N1+1;
  N2=N2+1;
end
Edges(EdgeNumber)=Edge([N1 1],Grid,EdgeNumber,EdgeNumber);
Grid.Edges=Edges;
Grid.NumEdges=NumEdges;

NumFaces=1;
Faces(1:NumFaces)=Face([0 0 0 0],Grid,0,'x');
Type='o';
[Faces(1),Grid]=Face(1:NumEdges,Grid,1,Type,OrientFace);
Grid.Faces=Faces;
Grid.NumFaces=NumFaces;
end
