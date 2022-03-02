function [Grid]=GridToKiteGridNeu(Grid,OrientFace)
KiteGrid.nBar=[ 0  1   0   -1  
                1  0  -1    0];
KiteGrid.Dim=3;
KiteGrid.Type='Quad';

for iF=1:Grid.NumFaces
  Grid.Faces(iF).Kite.N=zeros(1,2*size(Grid.Faces(iF).N,2)+1);
  Grid.Faces(iF).Kite.N(1,1:size(Grid.Faces(iF).N,2))=...
    Grid.Faces(iF).N;
  Grid.Faces(iF).Kite.NumNodes=size(Grid.Faces(iF).N,2);
  Grid.Faces(iF).Kite.E=zeros(1,3*size(Grid.Faces(iF).N,2));
  Grid.Faces(iF).Kite.NumEdges=0;
  Grid.Faces(iF).Kite.F=zeros(1,size(Grid.Faces(iF).N,2));
  Grid.Faces(iF).Kite.NumFaces=0;
end
  

NumNodes=Grid.NumNodes+Grid.NumEdges+Grid.NumFaces;
Nodes(1:NumNodes)=Node([0 0 0],0);
NodeNumber=1;
for iN=1:Grid.NumNodes
  Nodes(NodeNumber)=Node(Grid.Nodes(iN).P,NodeNumber);
  NodeNumber=NodeNumber+1;
end
for iE=1:Grid.NumEdges
  P=0.5*(Grid.Nodes(Grid.Edges(iE).N(1)).P...
    +Grid.Nodes(Grid.Edges(iE).N(2)).P);
  Nodes(NodeNumber)=Node(P,NodeNumber);
  for iFLoc=1:size(Grid.Edges(iE).F,2)
    iF=Grid.Edges(iE).F(iFLoc);
    Grid.Faces(iF).Kite.NumNodes=Grid.Faces(iF).Kite.NumNodes+1;
    Grid.Faces(iF).Kite.N(1,Grid.Faces(iF).Kite.NumNodes)=NodeNumber;
  end
  NodeNumber=NodeNumber+1;
end
for iF=1:Grid.NumFaces
  P=[0 0 0];
  for iN=1:size(Grid.Faces(iF).N,2)
    P=P+Grid.Nodes(Grid.Faces(iF).N(iN)).P;
  end
  P=(1/size(Grid.Faces(iF).N,2))*P;
  Nodes(NodeNumber)=Node(P,NodeNumber);
  Grid.Faces(iF).Kite.NumNodes=Grid.Faces(iF).Kite.NumNodes+1;
  Grid.Faces(iF).Kite.N(1,Grid.Faces(iF).Kite.NumNodes)=NodeNumber;
  NodeNumber=NodeNumber+1;
end
KiteGrid.Nodes=Nodes;
KiteGrid.NumNodes=NumNodes;

NumEdges=2*Grid.NumEdges;
for iF=1:Grid.NumFaces
  NumEdges=NumEdges+size(Grid.Faces(iF).N,2);
end

Edges(1:NumEdges)=Edge([1 2],KiteGrid,0,0);
EdgeNumber=1;
for iE=1:Grid.NumEdges
  N1=Grid.Edges(iE).N(1);
  N2=iE+Grid.NumNodes;
  N3=Grid.Edges(iE).N(2);
  Edges(EdgeNumber)=Edge([N1 N2],KiteGrid,EdgeNumber,EdgeNumber);
  EdgeNumber=EdgeNumber+1;
  Edges(EdgeNumber)=Edge([N2 N3],KiteGrid,EdgeNumber,EdgeNumber);
%   for iFLoc=1:size(Grid.Edges(iE).F,2)
%     iF=Grid.Edges(iE).F(iFLoc);
%     Grid.Faces(iF).Kite.NumEdges=Grid.Faces(iF).Kite.NumEdges+1;
%     Grid.Faces(iF).Kite.E(1,Grid.Faces(iF).Kite.NumEdges)=EdgeNumber-1;
%     Grid.Faces(iF).Kite.NumEdges=Grid.Faces(iF).Kite.NumEdges+1;
%     Grid.Faces(iF).Kite.E(1,Grid.Faces(iF).Kite.NumEdges)=EdgeNumber;
%   end
  EdgeNumber=EdgeNumber+1;
end
for iF=1:Grid.NumFaces
  for i=1:size(Grid.Faces(iF).E,2)
    iE=Grid.Faces(iF).E(i);
    if Grid.Edges(iE).N(1)==Grid.Faces(iF).N(i);
      Grid.Faces(iF).Kite.E(2*i-1)=2*iE-1;
      Grid.Faces(iF).Kite.E(2*i)=2*iE;
    else
      Grid.Faces(iF).Kite.E(2*i-1)=2*iE;
      Grid.Faces(iF).Kite.E(2*i)=2*iE-1;
    end
  end
  Grid.Faces(iF).Kite.NumEdges=2*size(Grid.Faces(iF).E,2);
end
NumInnerNode=Grid.NumNodes+Grid.NumEdges+1;
for iF=1:Grid.NumFaces
  N2=NumInnerNode;
  for iE=1:size(Grid.Faces(iF).E,2)
    N1=Grid.NumNodes+Grid.Edges(Grid.Faces(iF).E(iE)).E;
    Edges(EdgeNumber)=Edge([N1 N2],KiteGrid,EdgeNumber,EdgeNumber);
    Grid.Faces(iF).Kite.NumEdges=Grid.Faces(iF).Kite.NumEdges+1;
    Grid.Faces(iF).Kite.E(1,Grid.Faces(iF).Kite.NumEdges)=EdgeNumber;
    EdgeNumber=EdgeNumber+1;
  end
  NumInnerNode=NumInnerNode+1;
  Grid.Faces(iF).Kite.NumEdgesI=Grid.Faces(iF).Kite.NumEdges;
  Grid.Faces(iF).Kite.NumEdgesB=0;
end

KiteGrid.Edges=Edges;
KiteGrid.NumEdges=NumEdges;
KiteGrid.NumEdgesI=NumEdges;
KiteGrid.NumEdgesB=0;
NumFaces=0;
for iF=1:Grid.NumFaces
  NumFaces=NumFaces+size(Grid.Faces(iF).N,2);
end
Faces(1:NumFaces)=Face([0 0 0 0],KiteGrid,0,'x');
FaceNumber=1;
NumInnerEdge=2*Grid.NumEdges+1;
Type='o';
for iF=1:Grid.NumFaces
  for iE=1:size(Grid.Faces(iF).E,2)-1
    if Grid.Edges(Grid.Faces(iF).E(iE)).N(1)...
        ==Grid.Faces(iF).N(iE+1)
      E1=2*Grid.Edges(Grid.Faces(iF).E(iE)).E-1;
    else
      E1=2*Grid.Edges(Grid.Faces(iF).E(iE)).E;
    end
    if Grid.Edges(Grid.Faces(iF).E(iE+1)).N(1)...
        ==Grid.Faces(iF).N(iE+1)
      E2=2*Grid.Edges(Grid.Faces(iF).E(iE+1)).E-1;
    else
      E2=2*Grid.Edges(Grid.Faces(iF).E(iE+1)).E;
    end
    E3=NumInnerEdge+1;
    E4=NumInnerEdge;
    [Faces(FaceNumber),KiteGrid]=Face([E1 E2 E3 E4],KiteGrid,...
      FaceNumber,Type,OrientFace);
    Grid.Faces(iF).Kite.NumFaces=Grid.Faces(iF).Kite.NumFaces+1;
    Grid.Faces(iF).Kite.F(1,Grid.Faces(iF).Kite.NumFaces)=FaceNumber;
    FaceNumber=FaceNumber+1;
    NumInnerEdge=NumInnerEdge+1;
  end
  iE=size(Grid.Faces(iF).E,2);
  if Grid.Edges(Grid.Faces(iF).E(iE)).N(1)...
      ==Grid.Faces(iF).N(1)
    E1=2*Grid.Edges(Grid.Faces(iF).E(iE)).E-1;
  else
    E1=2*Grid.Edges(Grid.Faces(iF).E(iE)).E;
  end
  if Grid.Edges(Grid.Faces(iF).E(1)).N(1)...
      ==Grid.Faces(iF).N(1)
    E2=2*Grid.Edges(Grid.Faces(iF).E(1)).E-1;
  else
    E2=2*Grid.Edges(Grid.Faces(iF).E(1)).E;
  end
  E3=NumInnerEdge+1-size(Grid.Faces(iF).E,2);
  E4=NumInnerEdge;
  [Faces(FaceNumber),KiteGrid]=Face([E1 E2 E3 E4],KiteGrid,...
    FaceNumber,Type,OrientFace);
  Grid.Faces(iF).Kite.NumFaces=Grid.Faces(iF).Kite.NumFaces+1;
  Grid.Faces(iF).Kite.F(1,Grid.Faces(iF).Kite.NumFaces)=FaceNumber;
  FaceNumber=FaceNumber+1;
  NumInnerEdge=NumInnerEdge+1;
end
KiteGrid.Faces=Faces;
KiteGrid.NumFaces=NumFaces;


KiteGrid=Orientation(KiteGrid);
KiteGrid=Renumbering(KiteGrid);
Grid.KiteGrid=KiteGrid;
end