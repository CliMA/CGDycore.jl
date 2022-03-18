function [Grid] = PolygonGrid(np,OrientFace)
dphi=2*pi/np;

phi=0;
NumNodes=np;
Nodes(1:NumNodes)=Node([0 0 0],0);
NodeNumber=1;
for i=1:np
  x=sin(phi);
  y=-cos(phi);
  z=0;
  Nodes(NodeNumber)=Node([x;y;z],NodeNumber);
  NodeNumber=NodeNumber+1;
  phi=phi+dphi;
end
Grid.Nodes=Nodes;
Grid.NumNodes=NumNodes;

NumEdges=np;
Edges(1:NumEdges)=Edge([1 2],Grid,0,0);
EdgeNumber=1;
N1=1;
N2=2;
for i=1:np-1
  Edges(EdgeNumber)=Edge([N1 N2],Grid,EdgeNumber,EdgeNumber);
  EdgeNumber=EdgeNumber+1;
  N1=N2;
  N2=N2+1;
end
Edges(EdgeNumber)=Edge([N1 1],Grid,EdgeNumber,EdgeNumber);
Grid.Edges=Edges;
Grid.NumEdges=NumEdges;
Grid.NumEdgesI=0;
Grid.NumEdgesB=np;

NumFaces=1;
Faces(1:NumFaces)=Face([0 0 0 0],Grid,0,'x');
FaceNumber=1;
Type='o';
Faces(1)=Face(1:np,Grid,FaceNumber,Type,OrientFace);

Grid.Faces=Faces;
Grid.NumFaces=NumFaces;

end
