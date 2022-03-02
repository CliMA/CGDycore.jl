function Grid=CubedGrid(n,OrientFace,Param)
  Grid.nBar=[ 0  1   0   1 
             -1  0  -1   0];
  Grid.Dim=3;
  Grid.Type='Quad';
  Grid.Rad=Param.RadEarth;
  Grid.Form='Sphere';
  dd=2.0e0/n;
  NumNodes=(6*(n-1)*(n-1)+12*(n-1)+8);
  Nodes(1:NumNodes)=Node([0 0 0],0);
  NodeNumber=1;
  %Faces
  % West
  NodeNumberW=NodeNumber;
  x=-1.0;
  y=0.0;
  z=0.0;
  for k=1:n-1
    for j=1:n-1
      Nodes(NodeNumber)=Node(CubePoint(-1,j,k,n,x,y,z,Grid.Rad),NodeNumber);
      NodeNumber=NodeNumber+1;
    end
  end
  % ! East
  NodeNumberE=NodeNumber;
  x=1.0;
  y=0.0;
  z=0.0;
  for k=1:n-1
    for j=1:n-1
      Nodes(NodeNumber)=Node(CubePoint(-1,j,k,n,x,y,z,Grid.Rad),NodeNumber);
      NodeNumber=NodeNumber+1;
    end
  end
  % ! South
  NodeNumberS=NodeNumber;
  x=0.0;
  y=-1.0;
  z=0.0;
  for k=1:n-1
    for i=1:n-1
      Nodes(NodeNumber)=Node(CubePoint(i,-1,k,n,x,y,z,Grid.Rad),NodeNumber);
      NodeNumber=NodeNumber+1;
    end
  end
  % ! North
  NodeNumberN=NodeNumber;
  x=0.0;
  y=1.0;
  z=0.0;
  for k=1:n-1
    for i=1:n-1
      Nodes(NodeNumber)=Node(CubePoint(i,-1,k,n,x,y,z,Grid.Rad),NodeNumber);
      NodeNumber=NodeNumber+1;
    end
  end
  % ! Bottom
  NodeNumberB=NodeNumber;
  x=0.0;
  y=0.0;
  z=-1.0;
  for j=1:n-1
    for i=1:n-1
      Nodes(NodeNumber)=Node(CubePoint(i,j,-1,n,x,y,z,Grid.Rad),NodeNumber);
      NodeNumber=NodeNumber+1;
    end
  end
  % ! Top
  NodeNumberT=NodeNumber;
  x=0.0;
  y=0.0;
  z=1.0;
  for j=1:n-1
    for i=1:n-1
      Nodes(NodeNumber)=Node(CubePoint(i,j,-1,n,x,y,z,Grid.Rad),NodeNumber);
      NodeNumber=NodeNumber+1;
    end
  end
  % !Edges
  % !West East
  NodeNumberWEmm=NodeNumber;
  x=0.0;
  y=-1.0;
  z=-1.0;
  for i=1:n-1
    Nodes(NodeNumber)=Node(CubePoint(i,-1,-1,n,x,y,z,Grid.Rad),NodeNumber);
    NodeNumber=NodeNumber+1;
  end
  NodeNumberWEpm=NodeNumber;
  x=0.0;
  y=1.0;
  z=-1.0;
  for i=1:n-1
    Nodes(NodeNumber)=Node(CubePoint(i,-1,-1,n,x,y,z,Grid.Rad),NodeNumber);
    NodeNumber=NodeNumber+1;
  end
  NodeNumberWEmp=NodeNumber;
  x=0.0;
  y=-1.0;
  z=1.0;
  for i=1:n-1
    Nodes(NodeNumber)=Node(CubePoint(i,-1,-1,n,x,y,z,Grid.Rad),NodeNumber);
    NodeNumber=NodeNumber+1;
  end
  NodeNumberWEpp=NodeNumber;
  x=0.0;
  y=1.0;
  z=1.0;
  for i=1:n-1
    Nodes(NodeNumber)=Node(CubePoint(i,-1,-1,n,x,y,z,Grid.Rad),NodeNumber);
    NodeNumber=NodeNumber+1;
  end
  % South North
  NodeNumberSNmm=NodeNumber;
  x=-1.0;
  y=0.0;
  z=-1.0;
  for j=1:n-1
    Nodes(NodeNumber)=Node(CubePoint(-1,j,-1,n,x,y,z,Grid.Rad),NodeNumber);
    NodeNumber=NodeNumber+1;
  end
  NodeNumberSNpm=NodeNumber;
  x=1.0;
  y=0.0;
  z=-1.0;
  for j=1:n-1
    Nodes(NodeNumber)=Node(CubePoint(-1,j,-1,n,x,y,z,Grid.Rad),NodeNumber);
    NodeNumber=NodeNumber+1;
  end
  NodeNumberSNmp=NodeNumber;
  x=-1.0;
  y=0.0;
  z=1.0;
  for j=1:n-1
    Nodes(NodeNumber)=Node(CubePoint(-1,j,-1,n,x,y,z,Grid.Rad),NodeNumber);
    NodeNumber=NodeNumber+1;
  end
  NodeNumberSNpp=NodeNumber;
  x=1.0;
  y=0.0;
  z=1.0;
  for j=1:n-1
    Nodes(NodeNumber)=Node(CubePoint(-1,j,-1,n,x,y,z,Grid.Rad),NodeNumber);
    NodeNumber=NodeNumber+1;
  end
  % Bottom Top
   NodeNumberBTmm=NodeNumber;
   x=-1.0;
  y=-1.0;
  z=0.0;
  for k=1:n-1
    Nodes(NodeNumber)=Node(CubePoint(-1,-1,k,n,x,y,z,Grid.Rad),NodeNumber);
    NodeNumber=NodeNumber+1;
  end
  
   NodeNumberBTpm=NodeNumber;
    x=1.0;
  y=-1.0;
  z=0.0;
  for k=1:n-1
    Nodes(NodeNumber)=Node(CubePoint(-1,-1,k,n,x,y,z,Grid.Rad),NodeNumber);
    NodeNumber=NodeNumber+1;
  end
  
   NodeNumberBTmp=NodeNumber;
    x=-1.0;
  y=1.0;
  z=0.0;
  for k=1:n-1
    Nodes(NodeNumber)=Node(CubePoint(-1,-1,k,n,x,y,z,Grid.Rad),NodeNumber);
    NodeNumber=NodeNumber+1;
  end
  
  NodeNumberBTpp=NodeNumber;
   x=1.0;
  y=1.0;
  z=0.0;
  for k=1:n-1
    Nodes(NodeNumber)=Node(CubePoint(-1,-1,k,n,x,y,z,Grid.Rad),NodeNumber);
    NodeNumber=NodeNumber+1;
  end
  % Nodes
  NodeNumbermmm=NodeNumber;
  x=-1.0e0;
  y=-1.0e0;
  z=-1.0e0;
  Nodes(NodeNumber)=Node(CubePoint(-1,-1,-1,n,x,y,z,Grid.Rad),NodeNumber);
  NodeNumber=NodeNumber+1;
  NodeNumberpmm=NodeNumber;
  x=1.0e0;
  y=-1.0e0;
  z=-1.0e0;
  Nodes(NodeNumber)=Node(CubePoint(-1,-1,-1,n,x,y,z,Grid.Rad),NodeNumber);
  NodeNumber=NodeNumber+1;
  NodeNumbermpm=NodeNumber;
  x=-1.0e0;
  y=1.0e0;
  z=-1.0e0;
  Nodes(NodeNumber)=Node(CubePoint(-1,-1,-1,n,x,y,z,Grid.Rad),NodeNumber);
  NodeNumber=NodeNumber+1;
  NodeNumberppm=NodeNumber;
  x=1.0e0;
  y=1.0e0;
  z=-1.0e0;
  Nodes(NodeNumber)=Node(CubePoint(-1,-1,-1,n,x,y,z,Grid.Rad),NodeNumber);
  NodeNumber=NodeNumber+1;
  NodeNumbermmp=NodeNumber;
  x=-1.0e0;
  y=-1.0e0;
  z=1.0e0;
  Nodes(NodeNumber)=Node(CubePoint(-1,-1,-1,n,x,y,z,Grid.Rad),NodeNumber);
  NodeNumber=NodeNumber+1;
  NodeNumberpmp=NodeNumber;
  x=1.0e0;
  y=-1.0e0;
  z=1.0e0;
  Nodes(NodeNumber)=Node(CubePoint(-1,-1,-1,n,x,y,z,Grid.Rad),NodeNumber);
  NodeNumber=NodeNumber+1;
  
  NodeNumbermpp=NodeNumber;
  x=-1.0e0;
  y=1.0e0;
  z=1.0e0;
  Nodes(NodeNumber)=Node(CubePoint(-1,-1,-1,n,x,y,z,Grid.Rad),NodeNumber);
  NodeNumber=NodeNumber+1;
  NodeNumberppp=NodeNumber;
  x=1.0e0;
  y=1.0e0;
  z=1.0e0;
  Nodes(NodeNumber)=Node(CubePoint(-1,-1,-1,n,x,y,z,Grid.Rad),NodeNumber);
  NodeNumber=NodeNumber+1;
  Grid.Nodes=Nodes;
  
  NumEdges=12*(n-1)*n+12*n;
  NumEdgesI=12*(n-1)*n+12*n;
  Edges(1:NumEdges)=Edge([1 2],Grid,0,0,'',0);
  EdgeNumber=1;
  % West
  [Edges,EdgeNumber,EdgeNumberW1,EdgeNumberW2]=InsertFaceEdge(n,EdgeNumber,NodeNumberW...
    ,NodeNumberBTmm,NodeNumberBTmp...
    ,NodeNumberSNmm,NodeNumberSNmp...
    ,Edges,Grid);
  
  % East
  [Edges,EdgeNumber,EdgeNumberE1,EdgeNumberE2]=InsertFaceEdge(n,EdgeNumber,NodeNumberE...
                     ,NodeNumberBTpm,NodeNumberBTpp...
                     ,NodeNumberSNpm,NodeNumberSNpp...
                     ,Edges,Grid);
  % South
  [Edges,EdgeNumber,EdgeNumberS1,EdgeNumberS2]=InsertFaceEdge(n,EdgeNumber,NodeNumberS...
                     ,NodeNumberBTmm,NodeNumberBTpm...
                     ,NodeNumberWEmm,NodeNumberWEmp...
                     ,Edges,Grid);
  % North
  [Edges,EdgeNumber,EdgeNumberN1,EdgeNumberN2]=InsertFaceEdge(n,EdgeNumber,NodeNumberN...
                     ,NodeNumberBTmp,NodeNumberBTpp...
                     ,NodeNumberWEpm,NodeNumberWEpp...
                     ,Edges,Grid);
  % Bottom
  [Edges,EdgeNumber,EdgeNumberB1,EdgeNumberB2]=InsertFaceEdge(n,EdgeNumber,NodeNumberB...
                     ,NodeNumberSNmm,NodeNumberSNpm...
                     ,NodeNumberWEmm,NodeNumberWEpm...
                     ,Edges,Grid);
   % Top
   [Edges,EdgeNumber,EdgeNumberT1,EdgeNumberT2]=InsertFaceEdge(n,EdgeNumber,NodeNumberT...
     ,NodeNumberSNmp,NodeNumberSNpp...
     ,NodeNumberWEmp,NodeNumberWEpp...
     ,Edges,Grid);
   
% Edges
%West East
  [Edges,EdgeNumber,EdgeNumberWEmm]=InsertEdgeEdge(n,EdgeNumber,NodeNumberWEmm,NodeNumbermmm,NodeNumberpmm,Edges,Grid);
  [Edges,EdgeNumber,EdgeNumberWEpm]=InsertEdgeEdge(n,EdgeNumber,NodeNumberWEpm,NodeNumbermpm,NodeNumberppm,Edges,Grid);
  [Edges,EdgeNumber,EdgeNumberWEmp]=InsertEdgeEdge(n,EdgeNumber,NodeNumberWEmp,NodeNumbermmp,NodeNumberpmp,Edges,Grid);
  [Edges,EdgeNumber,EdgeNumberWEpp]=InsertEdgeEdge(n,EdgeNumber,NodeNumberWEpp,NodeNumbermpp,NodeNumberppp,Edges,Grid);
% South North
  [Edges,EdgeNumber,EdgeNumberSNmm]=InsertEdgeEdge(n,EdgeNumber,NodeNumberSNmm,NodeNumbermmm,NodeNumbermpm,Edges,Grid);
  [Edges,EdgeNumber,EdgeNumberSNpm]=InsertEdgeEdge(n,EdgeNumber,NodeNumberSNpm,NodeNumberpmm,NodeNumberppm,Edges,Grid);
  [Edges,EdgeNumber,EdgeNumberSNmp]=InsertEdgeEdge(n,EdgeNumber,NodeNumberSNmp,NodeNumbermmp,NodeNumbermpp,Edges,Grid);
  [Edges,EdgeNumber,EdgeNumberSNpp]=InsertEdgeEdge(n,EdgeNumber,NodeNumberSNpp,NodeNumberpmp,NodeNumberppp,Edges,Grid);
% Bottom Top
  [Edges,EdgeNumber,EdgeNumberBTmm]=InsertEdgeEdge(n,EdgeNumber,NodeNumberBTmm,NodeNumbermmm,NodeNumbermmp,Edges,Grid);
  [Edges,EdgeNumber,EdgeNumberBTpm]=InsertEdgeEdge(n,EdgeNumber,NodeNumberBTpm,NodeNumberpmm,NodeNumberpmp,Edges,Grid);
  [Edges,EdgeNumber,EdgeNumberBTmp]=InsertEdgeEdge(n,EdgeNumber,NodeNumberBTmp,NodeNumbermpm,NodeNumbermpp,Edges,Grid);
  [Edges,EdgeNumber,EdgeNumberBTpp]=InsertEdgeEdge(n,EdgeNumber,NodeNumberBTpp,NodeNumberppm,NodeNumberppp,Edges,Grid);
   
  Grid.Edges=Edges;
  
  NumFaces=6*n*n;
  Faces(1:NumFaces)=Face([0 0 0 0],Grid,0,'x');
  
  FaceNumber=1;
% Faces
% West
  [Faces,FaceNumber,Grid]=InsertFaceFace(n,FaceNumber,'W',EdgeNumberW1,EdgeNumberW2...
                     ,EdgeNumberSNmm,EdgeNumberSNmp...
                     ,EdgeNumberBTmm,EdgeNumberBTmp,Grid,Faces,OrientFace);
                   
% East
  [Faces,FaceNumber,Grid]=InsertFaceFace(n,FaceNumber,'E',EdgeNumberE1,EdgeNumberE2 ...
                     ,EdgeNumberSNpm,EdgeNumberSNpp ...
                     ,EdgeNumberBTpm,EdgeNumberBTpp,Grid,Faces,OrientFace);
% South
  [Faces,FaceNumber,Grid]=InsertFaceFace(n,FaceNumber,'S',EdgeNumberS1,EdgeNumberS2 ...
                     ,EdgeNumberWEmm,EdgeNumberWEmp ...
                     ,EdgeNumberBTmm,EdgeNumberBTpm,Grid,Faces,OrientFace);
% North
  [Faces,FaceNumber,Grid]=InsertFaceFace(n,FaceNumber,'N',EdgeNumberN1,EdgeNumberN2 ...
                     ,EdgeNumberWEpm,EdgeNumberWEpp ...
                     ,EdgeNumberBTmp,EdgeNumberBTpp,Grid,Faces,OrientFace);
% Bottom
  [Faces,FaceNumber,Grid]=InsertFaceFace(n,FaceNumber,'B',EdgeNumberB1,EdgeNumberB2 ...
                     ,EdgeNumberWEmm,EdgeNumberWEpm ...
                     ,EdgeNumberSNmm,EdgeNumberSNpm,Grid,Faces,OrientFace);
% Top                   
  [Faces,FaceNumber,Grid]=InsertFaceFace(n,FaceNumber,'T',EdgeNumberT1,EdgeNumberT2 ...
                     ,EdgeNumberWEmp,EdgeNumberWEpp ...
                     ,EdgeNumberSNmp,EdgeNumberSNpp,Grid,Faces,OrientFace);
                   
  Grid.Faces=Faces;
  
  Grid.NumNodes=size(Grid.Nodes,2);
  Grid.NumEdges=size(Grid.Edges,2);
  Grid.NumEdgesI=size(Grid.Edges,2);
  Grid.NumEdgesB=0;
  Grid.NumFaces=size(Grid.Faces,2);
  Grid.Dim=3;
  Grid=Orientation(Grid);
  Grid=Renumbering(Grid);
  Grid=FacesInNodes(Grid);

end

function N=CubePoint(i1,i2,i3,n,x,y,z,Rad)

if i1>0
  x=tan(i1*pi/(2*n)-0.25*pi);
end
if i2>0
  y=tan(i2*pi/(2*n)-0.25*pi);
end
if i3>0
  z=tan(i3*pi/(2*n)-0.25*pi);
end
N=[x;y;z];
N=N/norm(N)*Rad;
end

function [Edges,EdgeNumber,EdgeNumberStart1,EdgeNumberStart2]=InsertFaceEdge(n,EdgeNumber,NodeNumberStart...
                         ,NodeNumberE1Start1,NodeNumberE2Start1...
                         ,NodeNumberE1Start2,NodeNumberE2Start2,Edges,Grid)
  
  NodeNumber=NodeNumberStart;
  NodeNumberE1=NodeNumberE1Start1;
  NodeNumberE2=NodeNumberE2Start1;
  EdgeNumberStart1=EdgeNumber;
  for j=1:n-1
    for i=1:n
      N1=NodeNumber-1;
      N2=NodeNumber;
      if i==1 
        N1=NodeNumberE1;
        NodeNumberE1=NodeNumberE1+1;
      end  
      if i==n
        N2=NodeNumberE2;
        NodeNumberE2=NodeNumberE2+1;
        NodeNumber=NodeNumber-1;
      end
      Edges(EdgeNumber)=Edge([N1 N2],Grid,EdgeNumber,EdgeNumber,'X',EdgeNumber);
      EdgeNumber=EdgeNumber+1;
      NodeNumber=NodeNumber+1;
    end  
  end
  NodeNumber=NodeNumberStart;
  NodeNumberE1=NodeNumberE1Start2;
  NodeNumberE2=NodeNumberE2Start2;
  EdgeNumberStart2=EdgeNumber;
  for j=1:n
    for i=1:n-1
      N1=NodeNumber-(n-1);
      N2=NodeNumber;
      if j==1
        N1=NodeNumberE1;
        NodeNumberE1=NodeNumberE1+1;
      end  
      if j==n
        N2=NodeNumberE2;
        NodeNumberE2=NodeNumberE2+1;
      end
      Edges(EdgeNumber)=Edge([N1 N2],Grid,EdgeNumber,EdgeNumber,'X',EdgeNumber);
      EdgeNumber=EdgeNumber+1;
      NodeNumber=NodeNumber+1;
    end  
  end
end

function [Edges,EdgeNumber,EdgeNumberStart]=InsertEdgeEdge(n,EdgeNumber,NodeNumberStart,NodeNumberE1,NodeNumberE2,...
                         Edges,Grid)
  NodeNumber=NodeNumberStart;
  EdgeNumberStart=EdgeNumber;
  for i=1:n
    N1=NodeNumber-1;
    N2=NodeNumber;
    if i==1
      N1=NodeNumberE1;
    end  
    if i==n
      N2=NodeNumberE2;
    end 
    Edges(EdgeNumber)=Edge([N1 N2],Grid,EdgeNumber,EdgeNumber,'X',EdgeNumber);
    EdgeNumber=EdgeNumber+1;
    NodeNumber=NodeNumber+1;
  end  
end

function [Faces,FaceNumber,Grid]=InsertFaceFace(n,FaceNumber,Type,EdgeNumberStart1,EdgeNumberStart2...
                         ,EdgeNumberStartEW1,EdgeNumberStartEW2...
                         ,EdgeNumberStartSN1,EdgeNumberStartSN2,Grid,Faces,OrientFace)
  

  EdgeNumber1=EdgeNumberStart1;
  EdgeNumber2=EdgeNumberStart2;
  EdgeNumberEW1=EdgeNumberStartEW1;
  EdgeNumberEW2=EdgeNumberStartEW2;
  EdgeNumberSN1=EdgeNumberStartSN1;
  EdgeNumberSN2=EdgeNumberStartSN2;
  for j=1:n
    for i=1:n
      if (i==1 && j==1) || (i==1 && j==n)|| ...
          (i==n && j==1) || (i==n && j==n)
    
      end   
      %Face(FaceNumber)%Number=FaceNumber
      %Face(FaceNumber)%Boundary=1
      %ALLOCATE(Face(FaceNumber)%Edge(4))
      E1=EdgeNumber1-n;
      E2=EdgeNumber1;
      E3=EdgeNumber2-1;
      E4=EdgeNumber2;
      if j==1
        E1=EdgeNumberEW1;
        EdgeNumberEW1=EdgeNumberEW1+1;
      end
      if j==n
        E2=EdgeNumberEW2;
        EdgeNumberEW2=EdgeNumberEW2+1;
      end
      if i==1
        E3=EdgeNumberSN1;
        EdgeNumberSN1=EdgeNumberSN1+1;
      end
      if i==n
        E4=EdgeNumberSN2;
        EdgeNumberSN2=EdgeNumberSN2+1;
        EdgeNumber2=EdgeNumber2-1;
      end
      EdgeNumber1=EdgeNumber1+1;
      EdgeNumber2=EdgeNumber2+1;
      [Faces(FaceNumber),Grid]=Face([E1 E2 E3 E4],Grid,FaceNumber,Type,OrientFace);
      FaceNumber=FaceNumber+1;
    end
  end
end
