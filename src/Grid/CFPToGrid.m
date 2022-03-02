function [Grid] = CFPToGrid(Points,Faces,Cells,OrientFace)
Grid.Type='Hex';
NumNodes=size(Points,2);
Grid.Nodes(1:NumNodes)=Node([0 0 0],0);
NodeNumber=1;
for iN=1:NumNodes
  Grid.Nodes(NodeNumber)=Node([Points(iN).x;Points(iN).y;0],NodeNumber);
  NodeNumber=NodeNumber+1;
end
Grid.NumNodes=NumNodes;

NumEdges=size(Faces,2);
Grid.Edges(1:NumEdges)=Edge([1 2],Grid,0,0);

EdgeNumberI=0;
EdgeNumberB=0;
for iE=1:NumEdges
  if Faces(iE).Boundary==1
    EdgeNumberB=EdgeNumberB+1;
  else
    EdgeNumberI=EdgeNumberI+1;
  end
end
Grid.NumEdgesI=EdgeNumberI;
Grid.NumEdgesB=EdgeNumberB;

EdgeNumber=1;
EdgeNumberI=1;
EdgeNumberB=1;
EdgeOldToNew=zeros(1,NumEdges);
for iE=1:NumEdges
  if Faces(iE).Boundary==1
    Grid.Edges(EdgeNumber)=Edge(Faces(iE).P',Grid,EdgeNumber,...
      EdgeNumberB+Grid.NumEdgesI);
    EdgeOldToNew(iE)=EdgeNumberB+Grid.NumEdgesI;
    EdgeNumberB=EdgeNumberB+1;
  else
    Grid.Edges(EdgeNumber)=Edge(Faces(iE).P',Grid,EdgeNumber,...
      EdgeNumberI);
    EdgeOldToNew(iE)=EdgeNumberI;
    EdgeNumberI=EdgeNumberI+1;
  end
  EdgeNumber=EdgeNumber+1;
end
Grid.NumEdges=NumEdges;

NumFaces=size(Cells,2);
Faces1(1:NumFaces)=Face([0 0 0 0],Grid,0,'x');
FaceNumber=1;
Type='o';
for iF=1:NumFaces
  [Faces1(FaceNumber),Grid]=Face(Cells(iF).F,Grid,...
    FaceNumber,Type,OrientFace);
  FaceNumber=FaceNumber+1;
end
Grid.NumFaces=NumFaces;
Grid.Faces=Faces1;
end

