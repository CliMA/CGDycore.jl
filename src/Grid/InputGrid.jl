function [Grid] = InputGrid(filename,OrientFace)
Grid.nBar=[ 0  1   0    1  
           -1  0  -1    0];
Grid.Dim=3;
Grid.Type='Quad';
fileID=fopen(filename);
tline=fgetl(fileID);
tline=fgetl(fileID);
NumNodes=str2num(tline);
Grid.NumNodes=NumNodes;
Nodes(1:NumNodes)=Node([0 0 0],0);
for i=1:NumNodes
  tline = fgetl(fileID);
  tt=str2num(tline);
  NodeNumber=tt(1);
  x=tt(2);
  y=tt(3);
  Nodes(NodeNumber)=Node([x;y;0],NodeNumber);
end
Grid.Nodes=Nodes;
tline=fgetl(fileID);
tline=fgetl(fileID);
NumEdges=str2num(tline);
Grid.NumEdges=NumEdges;
Edges(1:NumEdges)=Edge([1 2],Grid,0,0);
EdgeNumberB=1;
EdgeNumberI=1;
EdgeType=zeros(NumEdges,1);
for i=1:NumEdges
  tline=fgetl(fileID);
  tt=str2num(tline);
  EdgeNumber=tt(1);
  N1=tt(2);
  N2=tt(3);
  EdgeType(i)=tt(4);
  if EdgeType(i)==1
    Edges(EdgeNumber)=Edge([N1 N2],Grid,EdgeNumber,EdgeNumberB);
    EdgeNumberB=EdgeNumberB+1;
  else
    Edges(EdgeNumber)=Edge([N1 N2],Grid,EdgeNumber,EdgeNumberI);
    EdgeNumberI=EdgeNumberI+1;
  end
end
for i=1:NumEdges
  if EdgeType(i)==1
    Edges(i).EI=Edges(i).EI+EdgeNumberI-1;
  end
end
Grid.Edges=Edges;
Grid.NumEdgesI=EdgeNumberI-1;
Grid.NumEdgesB=EdgeNumberB-1;
tline=fgetl(fileID);
tline=fgetl(fileID);
NumFaces=str2num(tline);
Grid.NumFaces=NumFaces;
Faces(1:NumFaces)=Face([0 0 0 0],Grid,0,'x');
Type='o';
c=zeros(NumFaces,1);
for i=1:NumFaces
  c(i)=i;
  tline=fgetl(fileID);
  tt=str2num(tline);
  FaceNumber=tt(1);
  [Faces(FaceNumber),Grid]=Face(tt(2:end),Grid,FaceNumber,Type,OrientFace);
end
Grid.Faces=Faces;
%PlotFaceGrid(c,Grid,1)
fclose(fileID);
Grid=Orientation(Grid);
Grid=Renumbering(Grid);
end

