function [Grid] = HexagonGridNeu(nx,ny,r,OrientFace)
dphi=2*pi/6;
dx=sqrt(3)*r;
dy=2*r;
x0=sqrt(3)/2*r;
y0=r;


NumNodes=2*(2*nx+1);
Nodes(1:NumNodes)=Node([0 0 0],0);
NodeNumber=1;
yShift=y0;
for j=1:ny
  phi=0;
  xShift=x0;
  for i=1:nx
    phi=phi-dphi;
    x=sin(phi)+xShift;
    y=-cos(phi)+yShift;
    z=0;
    Nodes(NodeNumber)=Node([x;y;z],NodeNumber);
    NodeNumber=NodeNumber+1;
    phi=phi+dphi;
    x=sin(phi)+xShift;
    y=-cos(phi)+yShift;
    z=0;
    Nodes(NodeNumber)=Node([x;y;z],NodeNumber);
    NodeNumber=NodeNumber+1;
    xShift=xShift+dx;
  end
  phi=phi-dphi;
  x=sin(phi)+xShift;
  y=-cos(phi)+yShift;
  z=0;
  Nodes(NodeNumber)=Node([x;y;z],NodeNumber);
  NodeNumber=NodeNumber+1;
  
  xShift=x0;
  phi=pi;
  yShift=yShift;
  for i=1:nx
    phi=phi+dphi;
    x=sin(phi)+xShift;
    y=-cos(phi)+yShift;
    z=0;
    Nodes(NodeNumber)=Node([x;y;z],NodeNumber);
    NodeNumber=NodeNumber+1;
    phi=phi-dphi;
    x=sin(phi)+xShift;
    y=-cos(phi)+yShift;
    z=0;
    Nodes(NodeNumber)=Node([x;y;z],NodeNumber);
    NodeNumber=NodeNumber+1;
    xShift=xShift+dx;
  end
  phi=phi+dphi;
    x=sin(phi)+xShift;
    y=-cos(phi)+yShift;
    z=0;
    Nodes(NodeNumber)=Node([x;y;z],NodeNumber);
    NodeNumber=NodeNumber+1;
end


Grid.Nodes=Nodes;
Grid.NumNodes=NumNodes;

NumEdges=2*2*nx+nx+1;
Edges(1:NumEdges)=Edge([1 2],Grid,0,0);
EdgeNumber=1;

N1=1;
N2=N1+1;
for i=1:2*nx
  Edges(EdgeNumber)=Edge([N1 N2],Grid,EdgeNumber,EdgeNumber);
  EdgeNumber=EdgeNumber+1;
  N1=N2;
  N2=N2+1;
end
N1=2*nx+2;
N2=N1+1;
for j=1:ny
  for i=1:2*nx
    Edges(EdgeNumber)=Edge([N1 N2],Grid,EdgeNumber,EdgeNumber);
    EdgeNumber=EdgeNumber+1;
    N1=N2;
    N2=N2+1;
  end
end
for j=1:ny
  N1=1;
  N2=N1+2*nx+1;
  for i=1:nx+1
    Edges(EdgeNumber)=Edge([N1 N2],Grid,EdgeNumber,EdgeNumber);
    EdgeNumber=EdgeNumber+1;
    N1=N1+2;
    N2=N2+2;
  end
end
Grid.Edges=Edges;
Grid.NumEdges=NumEdges;
Grid.NumEdgesI=NumEdges;
Grid.NumEdgesB=0;

NumFaces=nx;
Faces(1:NumFaces)=Face([0 0 0 0],Grid,0,'x');
FaceNumber=1;
Type='o';
E1=1;
E2=E1+1;
E3=2*2*nx+2;
E4=E2+2*nx;
E5=E4-1;
E6=E3-1;
for i=1:nx
  Faces(FaceNumber)=Face([E1 E2 E3 E4 E5 E6],Grid,FaceNumber,Type,OrientFace);
  E1=E1+2;
  E2=E2+2;
  E3=E3+1;
  E4=E4+2;
  E5=E5+2;
  E6=E6+1;
  FaceNumber=FaceNumber+1;
end
Grid.Faces=Faces;
Grid.NumFaces=NumFaces;

Grid=FacesInEdges(Grid);
Grid=FacesInNodes(Grid);

end
