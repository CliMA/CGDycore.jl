function Grid=TriGrid(nx,ny,lx,ly,Case)
switch Case
  case 'Tri2'
    Grid=TriGrid2(nx,ny,lx,ly);
  case 'Tri4'
    Grid=TriGrid4(nx,ny,lx,ly);
end
Grid=RenumberingTri(Grid);
Grid=OrientationFaceTri(Grid);
Grid=OrdEdgeFacesTri(Grid);
Grid.Dim=3;
Grid.Type='Triangle';
Grid.nBar=[ 0   1   1
           -1   1   0];
end
function Grid=TriGrid2(nx,ny,lx,ly)
NumNodes=(nx+1)*(ny+1);
Nodes(1:NumNodes)=Node([0 0 0],0);
NodeNumber=1;
dx=lx/nx;
dy=ly/ny;
y=0;
for iy=1:ny+1
  x=0;
  for ix=1:nx+1
    Nodes(NodeNumber)=Node([x;y;0],NodeNumber);
    NodeNumber=NodeNumber+1;
    x=x+dx;
  end
  y=y+dy;
end
Grid.Nodes=Nodes;
Grid.NumNodes=NumNodes;

NumEdges=(nx+1)*ny+nx*(ny+1)+nx*ny;
Edges(1:NumEdges)=Edge([1 2],Grid,0);

N1=1;
N2=1+nx+1;
EdgeNumber=1;
for iy=1:ny
  for ix=1:nx+1
    Edges(EdgeNumber)=Edge([N1 N2],Grid,EdgeNumber);
    EdgeNumber=EdgeNumber+1;
    N1=N1+1;
    N2=N2+1;
  end
end
N1=1;
N2=2;
for iy=1:ny+1
  for ix=1:nx
    Edges(EdgeNumber)=Edge([N1 N2],Grid,EdgeNumber);
    EdgeNumber=EdgeNumber+1;
    N1=N1+1;
    N2=N2+1;
  end
  N1=N1+1;
  N2=N2+1;
end
N1=2;
N2=nx+2;
for iy=1:ny
  for ix=1:nx
    Edges(EdgeNumber)=Edge([N1 N2],Grid,EdgeNumber);
    EdgeNumber=EdgeNumber+1;
    N1=N1+1;
    N2=N2+1;
  end
  N1=N1+1;
  N2=N2+1;
end
Grid.Edges=Edges;
Grid.NumEdges=NumEdges;

NumFaces=2*nx*ny;
Faces(1:NumFaces)=Face([0 0 0 0],Grid,0,'x');
E1=(nx+1)*ny+1;
E2=2;
E3=(nx+1)*ny+1+nx;
E4=1;
E5=(nx+1)*ny+nx*(ny+1)+1;
FaceNumber=1;
Type='o';
for iy=1:ny
  for ix=1:nx
    [Faces(FaceNumber),Grid]=Face([E1 E5 E4],Grid,FaceNumber,Type);
    FaceNumber=FaceNumber+1;
    [Faces(FaceNumber),Grid]=Face([E2 E3 E5],Grid,FaceNumber,Type);
    FaceNumber=FaceNumber+1;
    E1=E1+1;
    E2=E2+1;
    E3=E3+1;
    E4=E4+1;
    E5=E5+1;
  end
  E2=E2+1;
  E4=E4+1;
end
Grid.Faces=Faces;
Grid.NumFaces=NumFaces;
end

function Grid=TriGrid4(nx,ny,lx,ly)
  
  NumNodes=(nx+1)*(ny+1)+nx*ny;
  Nodes(1:NumNodes)=Node([0 0 0],0);
  NodeNumber=1;
  dx=lx/nx;
  dy=ly/ny;
  y=0;
  for iy=1:ny+1
    x=0;
    for ix=1:nx+1
      Nodes(NodeNumber)=Node([x;y;0],NodeNumber);
      NodeNumber=NodeNumber+1;
      x=x+dx;
    end
    y=y+dy;
  end
  y=0.5*dy;
  for iy=1:ny
    x=0.5*dx;
    for ix=1:nx
      Nodes(NodeNumber)=Node([x;y;0],NodeNumber);
      NodeNumber=NodeNumber+1;
      x=x+dx;
    end
    y=y+dy;
  end
  Grid.Nodes=Nodes;
  Grid.NumNodes=NumNodes;
  
  NumEdges=(nx+1)*ny+nx*(ny+1)+4*nx*ny;
  Edges(1:NumEdges)=Edge([1 2],Grid,0);

  N1=1;
  N2=1+nx+1;
  EdgeNumber=1;
  for iy=1:ny
    for ix=1:nx+1
      Edges(EdgeNumber)=Edge([N1 N2],Grid,EdgeNumber);
      EdgeNumber=EdgeNumber+1;
      N1=N1+1;
      N2=N2+1;
    end
  end
  N1=1;
  N2=2;
  for iy=1:ny+1
    for ix=1:nx
      Edges(EdgeNumber)=Edge([N1 N2],Grid,EdgeNumber);
      EdgeNumber=EdgeNumber+1;
      N1=N1+1;
      N2=N2+1;
    end
    N1=N1+1;
    N2=N2+1;
  end
  
  N1=1;
  N2=2;
  N3=nx+3;
  N4=nx+2;
  NC=(nx+1)*(ny+1)+1;
  for iy=1:ny
    for ix=1:nx
      Edges(EdgeNumber)=Edge([N1 NC],Grid,EdgeNumber);
      EdgeNumber=EdgeNumber+1;
      Edges(EdgeNumber)=Edge([N2 NC],Grid,EdgeNumber);
      EdgeNumber=EdgeNumber+1;
      Edges(EdgeNumber)=Edge([N3 NC],Grid,EdgeNumber);
      EdgeNumber=EdgeNumber+1;
      Edges(EdgeNumber)=Edge([N4 NC],Grid,EdgeNumber);
      EdgeNumber=EdgeNumber+1;
      N1=N1+1;
      N2=N2+1;
      N3=N3+1;
      N4=N4+1;
      NC=NC+1;
    end
    N1=N1+1;
    N2=N2+1;
    N3=N3+1;
    N4=N4+1;
  end
  Grid.Edges=Edges;
  Grid.NumEdges=NumEdges;
  
  NumFaces=4*nx*ny;
  Faces(1:NumFaces)=Face([0 0 0 0],Grid,0,'x');
  E1=(nx+1)*ny+1;
  E2=2;
  E3=(nx+1)*ny+1+nx;
  E4=1;
  EC1=(nx+1)*ny+nx*(ny+1)+1;
  EC2=(nx+1)*ny+nx*(ny+1)+2;
  EC3=(nx+1)*ny+nx*(ny+1)+3;
  EC4=(nx+1)*ny+nx*(ny+1)+4;
  FaceNumber=1;
  Type='o';
  for iy=1:ny
    for ix=1:nx
      [Faces(FaceNumber),Grid]=Face([EC1 E1 EC2],Grid,FaceNumber,Type);
      FaceNumber=FaceNumber+1;
      [Faces(FaceNumber),Grid]=Face([EC2 E2 EC3],Grid,FaceNumber,Type);
      FaceNumber=FaceNumber+1;
      [Faces(FaceNumber),Grid]=Face([EC3 E3 EC4],Grid,FaceNumber,Type);
      FaceNumber=FaceNumber+1;
      [Faces(FaceNumber),Grid]=Face([EC4 E4 EC1],Grid,FaceNumber,Type);
      FaceNumber=FaceNumber+1;
      E1=E1+1;
      E2=E2+1;
      E3=E3+1;
      E4=E4+1;
      EC1=EC1+4;
      EC2=EC2+4;
      EC3=EC3+4;
      EC4=EC4+4;
    end
    E2=E2+1;
    E4=E4+1;
  end
  Grid.Faces=Faces;
  Grid.NumFaces=NumFaces;
  
end

function [Grid] = RenumberingTri(Grid)
for iF=1:Grid.NumFaces
  Grid.Faces(iF)=RenumberingFace3(Grid.Faces(iF),Grid);
end
for iE=1:Grid.NumEdges
  Grid.Edges(iE)=PosEdgeInFace(Grid.Edges(iE),Grid);
end
end

function Edge=PosEdgeInFace(Edge,Grid)
Edge.FE=zeros(2,1);
for iF=1:size(Edge.F,2)
  F=Edge.F(iF);
  for iE=1:3
    if Edge.E==Grid.Faces(F).E(iE)
      Edge.FE(iF)=iE;
      break
    end
  end
end
end

function Face=RenumberingFace3(Face,Grid)
Face.N=sort(Face.N);
N=[Face.N Face.N(1)];
E=Face.E;
N2=N(1);
for iN=1:3
  N1=N2;
  N2=N(iN+1);
  for iE=1:3
    if (N1==Grid.Edges(E(iE)).N(1) && N2==Grid.Edges(E(iE)).N(2)) || ...
        (N1==Grid.Edges(E(iE)).N(2) && N2==Grid.Edges(E(iE)).N(1))
      Face.E(iN)=E(iE);
      break
    end
  end
end
end

function Grid=OrientationFaceTri(Grid)
for iF=1:Grid.NumFaces
  Face=Grid.Faces(iF);
  r=cross(Grid.Nodes(Face.N(1)).P-Grid.Nodes(Face.N(2)).P,Grid.Nodes(Face.N(3)).P-Grid.Nodes(Face.N(2)).P);
  if r(3)<0
    Grid.Faces(iF).Ori=1;
  else
    Grid.Faces(iF).Ori=-1;
  end
end
end

function Grid=OrdEdgeFacesTri(Grid)
for iE=1:Grid.NumEdges
  Edge=Grid.Edges(iE);
  if size(Edge.F,2)>1
    if iE==Grid.Faces(Edge.F(1)).E(3) && Grid.Faces(Edge.F(1)).Ori==1
      Grid.Edges(iE).F(1)=Edge.F(2);
      Grid.Edges(iE).F(2)=Edge.F(1);
      Grid.Edges(iE).FE(1)=Edge.FE(2);
      Grid.Edges(iE).FE(2)=Edge.FE(1);
    end
  end
end
end


