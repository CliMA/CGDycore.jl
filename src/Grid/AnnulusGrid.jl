function Grid=AnnulusGrid(nx,ny,RadI,RadO,OrientFace,Boundary,Param)
Grid.nBar=[ 0  1   0   1
  -1  0  -1   0];
Grid.nBar3=[ 0  1   0   1
  -1  0  -1   0
  0   0  0    0];
Grid.Dim=3;
Grid.Type='Quad';
Grid.Form='Planar';
NumNodes=(nx)*(ny+1);
Nodes(1:NumNodes)=Node([0 0 0],0);
NodeNumber=1;
NodeNumberLoc=1;

for iy=1:ny+1
  for ix=1:nx
    Nodes(NodeNumber)=Node([AnnulusPoint(ix,iy,nx,ny,RadI,RadO,Param);0],NodeNumberLoc);
    NodeNumber=NodeNumber+1;
    NodeNumberLoc=NodeNumberLoc+1;
  end
end

Grid.Nodes=Nodes;
Grid.NumNodes=NumNodes;

NumEdges=(nx)*ny+nx*(ny+1);
Edges(1:NumEdges)=Edge([1 2],Grid,0,0,'');

EdgeNumber=1;
EdgeNumberI=1;
EdgeNumberB=1;
for iy=1:ny
  for ix=1:nx
    if strcmp(Boundary.WE,'Period') && ix==nx+1
      %Edges(EdgeNumber)=Edge([N1 N2],Grid,EdgeNumberI-nx);
    elseif strcmp(Boundary.WE,'FreeSlip') && (ix==nx+1 || ix==1)
      %Edges(EdgeNumber)=Edge([N1 N2],Grid,EdgeNumberB);
      EdgeNumberB=EdgeNumberB+1;
    else
      EdgeNumberI=EdgeNumberI+1;
    end
    EdgeNumber=EdgeNumber+1;
  end
end
for iy=1:ny+1
  for ix=1:nx
    if iy==1 || iy==ny+1
      EdgeNumberB=EdgeNumberB+1;
    else
      EdgeNumberI=EdgeNumberI+1;
    end
    EdgeNumber=EdgeNumber+1;
  end
end
Grid.NumEdges=EdgeNumber-1;
Grid.NumEdgesI=EdgeNumberI-1;
Grid.NumEdgesB=EdgeNumberB-1;

N1=1;
N2=1+nx;
EdgeNumber=1;
EdgeNumberI=1;
EdgeNumberB=Grid.NumEdgesI+1;
for iy=1:ny
  for ix=1:nx
    if strcmp(Boundary.WE,'Period') && ix==nx+1
      Edges(EdgeNumber)=Edge([N1 N2],Grid,EdgeNumber-nx,EdgeNumberI-nx,'');
    elseif strcmp(Boundary.WE,'FreeSlip') && (ix==nx+1 || ix==1)
      Edges(EdgeNumber)=Edge([N1 N2],Grid,EdgeNumber,EdgeNumberB,'');
      EdgeNumberB=EdgeNumberB+1;
    else
      Edges(EdgeNumber)=Edge([N1 N2],Grid,EdgeNumber,EdgeNumberI,'');
      EdgeNumberI=EdgeNumberI+1;
    end
    EdgeNumber=EdgeNumber+1;
    N1=N1+1;
    N2=N2+1;
  end
end
N1=1;
N2=2;
for iy=1:ny+1
  if iy==1
    BC='BCB';
  elseif iy==ny+1
    BC='BCT';
  else
    BC='';
  end
  for ix=1:nx
    if ix<nx
      if iy==1 || iy==ny+1
        Edges(EdgeNumber)=Edge([N1 N2],Grid,EdgeNumber,EdgeNumberB,BC);
        EdgeNumberB=EdgeNumberB+1;
      else
        Edges(EdgeNumber)=Edge([N1 N2],Grid,EdgeNumber,EdgeNumberI,BC);
        EdgeNumberI=EdgeNumberI+1;
      end
      EdgeNumber=EdgeNumber+1;
      N1=N1+1;
      N2=N2+1;
    else
      if iy==1 || iy==ny+1
        Edges(EdgeNumber)=Edge([N1 N2-nx],Grid,EdgeNumber,EdgeNumberB,BC);
        EdgeNumberB=EdgeNumberB+1;
      else
        Edges(EdgeNumber)=Edge([N1 N2-nx],Grid,EdgeNumber,EdgeNumberI,BC);
        EdgeNumberI=EdgeNumberI+1;
      end
      EdgeNumber=EdgeNumber+1;
    end
  end
  N1=N1+1;
  N2=N2+1;
end
Grid.Edges=Edges;


NumFaces=nx*ny;
Faces(1:NumFaces)=Face([0 0 0 0],Grid,0,'x');
E1=nx*ny+1;
E2=2;
E3=nx*ny+1+nx;
E4=1;
FaceNumber=1;
Type='o';
for iy=1:ny
  for ix=1:nx
    if ix<nx
      [Faces(FaceNumber),Grid]=Face([E1 E2 E3 E4],Grid,FaceNumber,Type,OrientFace);
      FaceNumber=FaceNumber+1;
      E1=E1+1;
      E2=E2+1;
      E3=E3+1;
      E4=E4+1;
    else
      [Faces(FaceNumber),Grid]=Face([E1 E2-nx E3 E4],Grid,FaceNumber,Type,OrientFace);
      FaceNumber=FaceNumber+1;
      E1=E1+1;
      E3=E3+1;
    end
  end
  E2=E2+1;
  E4=E4+1;
end
Grid.Faces=Faces;
Grid.NumFaces=NumFaces;
Grid=Orientation(Grid);
Grid=Renumbering(Grid);
end
function N=AnnulusPoint(i1,i2,n1,n2,RadI,RadO,Param)
hS=str2func(Param.hS);
phi=(i1-1)*2*pi/n1;
s=(i2-1)/n2;
r=s*RadO+(1-s)*(RadI+hS(phi,Param));
x=cos(phi)*r;
y=sin(phi)*r;
N=[x;y];
end
