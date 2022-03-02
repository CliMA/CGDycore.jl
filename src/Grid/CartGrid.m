function Grid=CartGrid(nx,ny,lx,ly,x0,y0,OrientFace,Boundary,Param)
Grid.nBar=[ 0  1   0   1  
           -1  0  -1   0];
Grid.nBar3=[ 0  1   0   1  
           -1  0  -1   0
           0   0  0    0];         
Grid.Dim=3;
Grid.Type='Quad';
Grid.Form='Planar';
Pert=0.0;
PertX=0.2;
PertY=0.2;
if strcmp(Boundary.WE,'Period') && strcmp(Boundary.BT,'Period')
  NumNodes=nx*ny;
elseif strcmp(Boundary.WE,'Period')
  NumNodes=nx*(ny+1);
elseif strcmp(Boundary.BT,'Period')
  NumNodes=(nx+1)*ny;
else
  NumNodes=(nx+1)*(ny+1);
end
Nodes(1:NumNodes)=Node([0 0 0],0);
NodeNumber=1;
dx=lx/nx;
dy=ly/ny;
y=y0;
P=zeros(3,nx+1,ny+1);
for iy=1:ny+1
  eta=(iy-1)/ny;
  x=x0;
  for ix=1:nx+1
    P(1,ix,iy)=x;
    P(2,ix,iy)=y;
    P(3,ix,iy)=0;
    P(2,ix,iy)=eta*(ly+y0)+(1-eta)*(hS(x,Param)+y0);
    x=x+dx;
  end
  y=y+dy;
end

y=y0;
for iy=1:ny+1
  x=x0;
  if iy==ny+1 && strcmp(Boundary.BT,'Period')
    for ix=1:nx+1
      if ix==nx+1 && strcmp(Boundary.WE,'Period')
      else
        Nodes(NodeNumber)=Node([x;y;0],NodeNumber);
        NodeNumber=NodeNumber+1;
      end
      x=x+dx;
    end
  else
    for ix=1:nx+1
      if ix==nx+1 && strcmp(Boundary.WE,'Period')
      else
        Nodes(NodeNumber)=Node([x;y;0],NodeNumber);
        NodeNumber=NodeNumber+1;
      end
      x=x+dx;
    end
  end
  y=y+dy;
end

Grid.Nodes=Nodes;
Grid.NumNodes=NumNodes;

if strcmp(Boundary.WE,'Period' ) && strcmp(Boundary.BT,'Period')
  NumEdges=2*nx*ny;
  NumEdgesX=nx*ny;
elseif strcmp(Boundary.WE,'Period') 
  NumEdges=nx*ny+nx*(ny+1);
  NumEdgesX=nx*(ny+1);
else
  NumEdges=(nx+1)*ny+nx*(ny+1);
  %NumEdgesX=nx*ny;
end

Edges(1:NumEdges)=Edge([1 2],Grid,0,0,'',0);

EdgeNumber=1;
EdgeNumberX=1;
EdgeNumberY=1;
BC='';
if strcmp(Boundary.WE,'Period')
  N1=1;
  N2=nx+1;
else
  N1=1;
  N2=nx+2;
end
for iy=1:ny
  for ix=1:nx+1
    if ix==nx+1 && strcmp(Boundary.WE,'Period')
    else
      if iy==ny && strcmp(Boundary.BT,'Period')
        Edges(EdgeNumber)=Edge([N1 1+(ix-1)],Grid,EdgeNumber,EdgeNumber,'Y',EdgeNumberY);
        EdgeNumber=EdgeNumber+1;
        EdgeNumberY=EdgeNumberY+1;
        N1=N1+1;
        N2=N2+1;
      else
        Edges(EdgeNumber)=Edge([N1 N2],Grid,EdgeNumber,EdgeNumber,'Y',EdgeNumberY);
        EdgeNumber=EdgeNumber+1;
        EdgeNumberY=EdgeNumberY+1;
        N1=N1+1;
        N2=N2+1;
      end
    end
  end
end
N1=1;
N2=2;
for iy=1:ny+1
  if iy==ny+1 && strcmp(Boundary.BT,'Period')
  else
    for ix=1:nx
      if ix==nx && strcmp(Boundary.WE,'Period')
        Edges(EdgeNumber)=Edge([N1 1+(iy-1)*nx],Grid,EdgeNumber,EdgeNumber,'X',EdgeNumberX);
        EdgeNumber=EdgeNumber+1;
        EdgeNumberX=EdgeNumberX+1;
        N1=N1+1;
        N2=N2+1;
      else
        Edges(EdgeNumber)=Edge([N1 N2],Grid,EdgeNumber,EdgeNumber,'X',EdgeNumberX);
        EdgeNumber=EdgeNumber+1;
        EdgeNumberX=EdgeNumberX+1;
        N1=N1+1;
        N2=N2+1;
        if ix==nx
          N1=N1+1;
          N2=N2+1;
        end
      end
    end
  end
end

Grid.Edges=Edges;
Grid.NumEdges=NumEdges;


NumFaces=nx*ny;
Faces(1:NumFaces)=Face([0 0 0 0],Grid,0,'x',[]);
if strcmp(Boundary.WE,'Period')
  E1=nx*ny+1;
  E3=nx*ny+1+nx;
else
  E1=(nx+1)*ny+1;
  E3=(nx+1)*ny+1+nx;
end
E2=2;
E4=1;
FaceNumber=1;
Type='o';
for iy=1:ny
  if iy==ny && strcmp(Boundary.BT,'Period')
    for ix=1:nx
      if ix==nx && strcmp(Boundary.WE,'Period')
        [Faces(FaceNumber),Grid]=Face([E1 1+(iy-1)*nx NumEdgesX+1+(ix-1) E4],Grid,FaceNumber,Type,OrientFace,...
          [P(:,ix,iy) P(:,ix+1,iy) P(:,ix+1,iy+1) P(:,ix,iy+1)]);
        FaceNumber=FaceNumber+1;
        E1=E1+1;
        E4=E4+1;
      else
        [Faces(FaceNumber),Grid]=Face([E1 E2 NumEdgesX+1+(ix-1) E4],Grid,FaceNumber,Type,OrientFace,...
          [P(:,ix,iy) P(:,ix+1,iy) P(:,ix+1,iy+1) P(:,ix,iy+1)]);
        FaceNumber=FaceNumber+1;
        E1=E1+1;
        E2=E2+1;
        E4=E4+1;
      end
    end
  else
    for ix=1:nx
      if ix==nx && strcmp(Boundary.WE,'Period')
        [Faces(FaceNumber),Grid]=Face([E1 1+(iy-1)*nx E3 E4],Grid,FaceNumber,Type,OrientFace,...
          [P(:,ix,iy) P(:,ix+1,iy) P(:,ix+1,iy+1) P(:,ix,iy+1)]);
        FaceNumber=FaceNumber+1;
        E1=E1+1;
        E3=E3+1;
      else
        [Faces(FaceNumber),Grid]=Face([E1 E2 E3 E4],Grid,FaceNumber,Type,OrientFace,...
          [P(:,ix,iy) P(:,ix+1,iy) P(:,ix+1,iy+1) P(:,ix,iy+1)]);
        FaceNumber=FaceNumber+1;
        E1=E1+1;
        E2=E2+1;
        E3=E3+1;
        E4=E4+1;
      end
    end
    E2=E2+1;
    E4=E4+1;
  end
end
Grid.Faces=Faces;
Grid.NumFaces=NumFaces;
Grid=Orientation(Grid);
Grid=Renumbering(Grid);
Grid=FacesInNodes(Grid);
end
