function CubedGrid(backend,FT,n,OrientFace,Rad,nz;order=true)
  nBar=[ 0  -1   0   1
        -1  0  1   0];
  Dim=3;
  Type=Quad()
  Rad=Rad;
  Form="Sphere";
  dd=2.0e0/n;

  NumNodes=(6*(n-1)*(n-1)+12*(n-1)+8);
  Nodes = map(1:NumNodes) do i
    Node()
  end
  NodeNumber=1;
  #Faces
  # West
  NodeNumberW=NodeNumber;
  x=-1.0;
  y=0.0;
  z=0.0;
  @inbounds for k=1:n-1
    @inbounds for j=1:n-1
      Nodes[NodeNumber]=Node(CubePoint(-1,j,k,n,x,y,z,Rad),NodeNumber,' ');
      NodeNumber=NodeNumber+1;
    end
  end
  # ! East
  NodeNumberE=NodeNumber;
  x=1.0;
  y=0.0;
  z=0.0;
  @inbounds for k=1:n-1
    @inbounds for j=1:n-1
      Nodes[NodeNumber]=Node(CubePoint(-1,j,k,n,x,y,z,Rad),NodeNumber,' ');
      NodeNumber=NodeNumber+1;
    end
  end
  # ! South
  NodeNumberS=NodeNumber;
  x=0.0;
  y=-1.0;
  z=0.0;
  @inbounds for k=1:n-1
    @inbounds for i=1:n-1
      Nodes[NodeNumber]=Node(CubePoint(i,-1,k,n,x,y,z,Rad),NodeNumber,' ');
      NodeNumber=NodeNumber+1;
    end
  end
  # ! North
  NodeNumberN=NodeNumber;
  x=0.0;
  y=1.0;
  z=0.0;
  @inbounds for k=1:n-1
    @inbounds for i=1:n-1
      Nodes[NodeNumber]=Node(CubePoint(i,-1,k,n,x,y,z,Rad),NodeNumber,' ');
      NodeNumber=NodeNumber+1;
    end
  end
  # ! Bottom
  NodeNumberB=NodeNumber;
  x=0.0;
  y=0.0;
  z=-1.0;
  @inbounds for j=1:n-1
    @inbounds for i=1:n-1
      Nodes[NodeNumber]=Node(CubePoint(i,j,-1,n,x,y,z,Rad),NodeNumber,' ');
      NodeNumber=NodeNumber+1;
    end
  end
  # ! Top
  NodeNumberT=NodeNumber;
  x=0.0;
  y=0.0;
  z=1.0;
  @inbounds for j=1:n-1
    @inbounds for i=1:n-1
      Nodes[NodeNumber]=Node(CubePoint(i,j,-1,n,x,y,z,Rad),NodeNumber,' ');
      NodeNumber=NodeNumber+1;
    end
  end
  # !Edges
  # !West East
  NodeNumberWEmm=NodeNumber;
  x=0.0;
  y=-1.0;
  z=-1.0;
  @inbounds for i=1:n-1
    Nodes[NodeNumber]=Node(CubePoint(i,-1,-1,n,x,y,z,Rad),NodeNumber,' ');
    NodeNumber=NodeNumber+1;
  end
  NodeNumberWEpm=NodeNumber;
  x=0.0;
  y=1.0;
  z=-1.0;
  @inbounds for i=1:n-1
    Nodes[NodeNumber]=Node(CubePoint(i,-1,-1,n,x,y,z,Rad),NodeNumber,' ');
    NodeNumber=NodeNumber+1;
  end
  NodeNumberWEmp=NodeNumber;
  x=0.0;
  y=-1.0;
  z=1.0;
  @inbounds for i=1:n-1
    Nodes[NodeNumber]=Node(CubePoint(i,-1,-1,n,x,y,z,Rad),NodeNumber,' ');
    NodeNumber=NodeNumber+1;
  end
  NodeNumberWEpp=NodeNumber;
  x=0.0;
  y=1.0;
  z=1.0;
  @inbounds for i=1:n-1
    Nodes[NodeNumber]=Node(CubePoint(i,-1,-1,n,x,y,z,Rad),NodeNumber,' ');
    NodeNumber=NodeNumber+1;
  end
  # South North
  NodeNumberSNmm=NodeNumber;
  x=-1.0;
  y=0.0;
  z=-1.0;
  @inbounds for j=1:n-1
    Nodes[NodeNumber]=Node(CubePoint(-1,j,-1,n,x,y,z,Rad),NodeNumber,' ');
    NodeNumber=NodeNumber+1;
  end
  NodeNumberSNpm=NodeNumber;
  x=1.0;
  y=0.0;
  z=-1.0;
  @inbounds for j=1:n-1
    Nodes[NodeNumber]=Node(CubePoint(-1,j,-1,n,x,y,z,Rad),NodeNumber,' ');
    NodeNumber=NodeNumber+1;
  end
  NodeNumberSNmp=NodeNumber;
  x=-1.0;
  y=0.0;
  z=1.0;
  @inbounds for j=1:n-1
    Nodes[NodeNumber]=Node(CubePoint(-1,j,-1,n,x,y,z,Rad),NodeNumber,' ');
    NodeNumber=NodeNumber+1;
  end
  NodeNumberSNpp=NodeNumber;
  x=1.0;
  y=0.0;
  z=1.0;
  @inbounds for j=1:n-1
    Nodes[NodeNumber]=Node(CubePoint(-1,j,-1,n,x,y,z,Rad),NodeNumber,' ');
    NodeNumber=NodeNumber+1;
  end
  # Bottom Top
   NodeNumberBTmm=NodeNumber;
   x=-1.0;
  y=-1.0;
  z=0.0;
  @inbounds for k=1:n-1
    Nodes[NodeNumber]=Node(CubePoint(-1,-1,k,n,x,y,z,Rad),NodeNumber,' ');
    NodeNumber=NodeNumber+1;
  end

   NodeNumberBTpm=NodeNumber;
    x=1.0;
  y=-1.0;
  z=0.0;
  @inbounds for k=1:n-1
    Nodes[NodeNumber]=Node(CubePoint(-1,-1,k,n,x,y,z,Rad),NodeNumber,' ');
    NodeNumber=NodeNumber+1;
  end

   NodeNumberBTmp=NodeNumber;
    x=-1.0;
  y=1.0;
  z=0.0;
  @inbounds for k=1:n-1
    Nodes[NodeNumber]=Node(CubePoint(-1,-1,k,n,x,y,z,Rad),NodeNumber,' ');
    NodeNumber=NodeNumber+1;
  end

  NodeNumberBTpp=NodeNumber;
   x=1.0;
  y=1.0;
  z=0.0;
  @inbounds for k=1:n-1
    Nodes[NodeNumber]=Node(CubePoint(-1,-1,k,n,x,y,z,Rad),NodeNumber,' ');
    NodeNumber=NodeNumber+1;
  end
  # Nodes
  NodeNumbermmm=NodeNumber;
  x=-1.0e0;
  y=-1.0e0;
  z=-1.0e0;
  Nodes[NodeNumber]=Node(CubePoint(-1,-1,-1,n,x,y,z,Rad),NodeNumber,' ');
  NodeNumber=NodeNumber+1;
  NodeNumberpmm=NodeNumber;
  x=1.0e0;
  y=-1.0e0;
  z=-1.0e0;
  Nodes[NodeNumber]=Node(CubePoint(-1,-1,-1,n,x,y,z,Rad),NodeNumber,' ');
  NodeNumber=NodeNumber+1;
  NodeNumbermpm=NodeNumber;
  x=-1.0e0;
  y=1.0e0;
  z=-1.0e0;
  Nodes[NodeNumber]=Node(CubePoint(-1,-1,-1,n,x,y,z,Rad),NodeNumber,' ');
  NodeNumber=NodeNumber+1;
  NodeNumberppm=NodeNumber;
  x=1.0e0;
  y=1.0e0;
  z=-1.0e0;
  Nodes[NodeNumber]=Node(CubePoint(-1,-1,-1,n,x,y,z,Rad),NodeNumber,' ');
  NodeNumber=NodeNumber+1;
  NodeNumbermmp=NodeNumber;
  x=-1.0e0;
  y=-1.0e0;
  z=1.0e0;
  Nodes[NodeNumber]=Node(CubePoint(-1,-1,-1,n,x,y,z,Rad),NodeNumber,' ');
  NodeNumber=NodeNumber+1;
  NodeNumberpmp=NodeNumber;
  x=1.0e0;
  y=-1.0e0;
  z=1.0e0;
  Nodes[NodeNumber]=Node(CubePoint(-1,-1,-1,n,x,y,z,Rad),NodeNumber,' ');
  NodeNumber=NodeNumber+1;

  NodeNumbermpp=NodeNumber;
  x=-1.0e0;
  y=1.0e0;
  z=1.0e0;
  Nodes[NodeNumber]=Node(CubePoint(-1,-1,-1,n,x,y,z,Rad),NodeNumber,' ');
  NodeNumber=NodeNumber+1;
  NodeNumberppp=NodeNumber;
  x=1.0e0;
  y=1.0e0;
  z=1.0e0;
  Nodes[NodeNumber]=Node(CubePoint(-1,-1,-1,n,x,y,z,Rad),NodeNumber,' ');
  NodeNumber=NodeNumber+1;

  NumEdges=12*(n-1)*n+12*n;
  NumEdgesI=12*(n-1)*n+12*n;
  Edges = map(1:NumEdges) do i
    Edge([1,2],Nodes,0,0,"",0);
  end
  EdgeNumber=1;
  # West
  (Edges,EdgeNumber,EdgeNumberW1,EdgeNumberW2)=InsertFaceEdge(n,EdgeNumber,NodeNumberW
    ,NodeNumberBTmm,NodeNumberBTmp
    ,NodeNumberSNmm,NodeNumberSNmp
    ,Edges,Nodes,Form,Rad);

  # East
  (Edges,EdgeNumber,EdgeNumberE1,EdgeNumberE2)=InsertFaceEdge(n,EdgeNumber,NodeNumberE
                     ,NodeNumberBTpm,NodeNumberBTpp
                     ,NodeNumberSNpm,NodeNumberSNpp
                     ,Edges,Nodes,Form,Rad);
  # South
  (Edges,EdgeNumber,EdgeNumberS1,EdgeNumberS2)=InsertFaceEdge(n,EdgeNumber,NodeNumberS
                     ,NodeNumberBTmm,NodeNumberBTpm
                     ,NodeNumberWEmm,NodeNumberWEmp
                     ,Edges,Nodes,Form,Rad);
  # North
  (Edges,EdgeNumber,EdgeNumberN1,EdgeNumberN2)=InsertFaceEdge(n,EdgeNumber,NodeNumberN
                     ,NodeNumberBTmp,NodeNumberBTpp
                     ,NodeNumberWEpm,NodeNumberWEpp
                     ,Edges,Nodes,Form,Rad);
  # Bottom
  (Edges,EdgeNumber,EdgeNumberB1,EdgeNumberB2)=InsertFaceEdge(n,EdgeNumber,NodeNumberB
                     ,NodeNumberSNmm,NodeNumberSNpm
                     ,NodeNumberWEmm,NodeNumberWEpm
                     ,Edges,Nodes,Form,Rad);
  # Top
  (Edges,EdgeNumber,EdgeNumberT1,EdgeNumberT2)=InsertFaceEdge(n,EdgeNumber,NodeNumberT
     ,NodeNumberSNmp,NodeNumberSNpp
     ,NodeNumberWEmp,NodeNumberWEpp
     ,Edges,Nodes,Form,Rad);

# Edges
#West East
  (Edges,EdgeNumber,EdgeNumberWEmm)=InsertEdgeEdge(n,EdgeNumber,NodeNumberWEmm,NodeNumbermmm,NodeNumberpmm,Edges,Nodes,Form,Rad);
  (Edges,EdgeNumber,EdgeNumberWEpm)=InsertEdgeEdge(n,EdgeNumber,NodeNumberWEpm,NodeNumbermpm,NodeNumberppm,Edges,Nodes,Form,Rad);
  (Edges,EdgeNumber,EdgeNumberWEmp)=InsertEdgeEdge(n,EdgeNumber,NodeNumberWEmp,NodeNumbermmp,NodeNumberpmp,Edges,Nodes,Form,Rad);
  (Edges,EdgeNumber,EdgeNumberWEpp)=InsertEdgeEdge(n,EdgeNumber,NodeNumberWEpp,NodeNumbermpp,NodeNumberppp,Edges,Nodes,Form,Rad);
# South North
  (Edges,EdgeNumber,EdgeNumberSNmm)=InsertEdgeEdge(n,EdgeNumber,NodeNumberSNmm,NodeNumbermmm,NodeNumbermpm,Edges,Nodes,Form,Rad);
  (Edges,EdgeNumber,EdgeNumberSNpm)=InsertEdgeEdge(n,EdgeNumber,NodeNumberSNpm,NodeNumberpmm,NodeNumberppm,Edges,Nodes,Form,Rad);
  (Edges,EdgeNumber,EdgeNumberSNmp)=InsertEdgeEdge(n,EdgeNumber,NodeNumberSNmp,NodeNumbermmp,NodeNumbermpp,Edges,Nodes,Form,Rad);
  (Edges,EdgeNumber,EdgeNumberSNpp)=InsertEdgeEdge(n,EdgeNumber,NodeNumberSNpp,NodeNumberpmp,NodeNumberppp,Edges,Nodes,Form,Rad);
# Bottom Top
  (Edges,EdgeNumber,EdgeNumberBTmm)=InsertEdgeEdge(n,EdgeNumber,NodeNumberBTmm,NodeNumbermmm,NodeNumbermmp,Edges,Nodes,Form,Rad);
  (Edges,EdgeNumber,EdgeNumberBTpm)=InsertEdgeEdge(n,EdgeNumber,NodeNumberBTpm,NodeNumberpmm,NodeNumberpmp,Edges,Nodes,Form,Rad);
  (Edges,EdgeNumber,EdgeNumberBTmp)=InsertEdgeEdge(n,EdgeNumber,NodeNumberBTmp,NodeNumbermpm,NodeNumbermpp,Edges,Nodes,Form,Rad);
  (Edges,EdgeNumber,EdgeNumberBTpp)=InsertEdgeEdge(n,EdgeNumber,NodeNumberBTpp,NodeNumberppm,NodeNumberppp,Edges,Nodes,Form,Rad);

  NumFaces=6*n*n;
  Faces = map(1:NumFaces) do i
    Face()
  end

  FaceNumber=1;
# Faces
# West
  (Faces,FaceNumber,Edges)=InsertFaceFace(n,FaceNumber,"W",EdgeNumberW1,EdgeNumberW2
                     ,EdgeNumberSNmm,EdgeNumberSNmp
                     ,EdgeNumberBTmm,EdgeNumberBTmp,Nodes,Edges,Faces,OrientFace,Form,Rad);

# East
  (Faces,FaceNumber,Edges)=InsertFaceFace(n,FaceNumber,"E",EdgeNumberE1,EdgeNumberE2
                     ,EdgeNumberSNpm,EdgeNumberSNpp
                     ,EdgeNumberBTpm,EdgeNumberBTpp,Nodes,Edges,Faces,OrientFace,Form,Rad);
# South
  (Faces,FaceNumber,Edges)=InsertFaceFace(n,FaceNumber,"S",EdgeNumberS1,EdgeNumberS2
                     ,EdgeNumberWEmm,EdgeNumberWEmp
                     ,EdgeNumberBTmm,EdgeNumberBTpm,Nodes,Edges,Faces,OrientFace,Form,Rad);
# North
  (Faces,FaceNumber,Edges)=InsertFaceFace(n,FaceNumber,"N",EdgeNumberN1,EdgeNumberN2
                     ,EdgeNumberWEpm,EdgeNumberWEpp
                     ,EdgeNumberBTmp,EdgeNumberBTpp,Nodes,Edges,Faces,OrientFace,Form,Rad);
# Bottom
  (Faces,FaceNumber,Edges)=InsertFaceFace(n,FaceNumber,"B",EdgeNumberB1,EdgeNumberB2
                     ,EdgeNumberWEmm,EdgeNumberWEpm
                     ,EdgeNumberSNmm,EdgeNumberSNpm,Nodes,Edges,Faces,OrientFace,Form,Rad);
# Top
  (Faces,FaceNumber,Edges)=InsertFaceFace(n,FaceNumber,"T",EdgeNumberT1,EdgeNumberT2
                     ,EdgeNumberWEmp,EdgeNumberWEpp
                     ,EdgeNumberSNmp,EdgeNumberSNpp,Nodes,Edges,Faces,OrientFace,Form,Rad);
  NumNodes=size(Nodes,1);
  NumEdges=size(Edges,1);
  NumEdgesI=size(Edges,1);
  NumEdgesB=0;
  NumFaces=size(Faces,1);
  Dim=3;
  if order
    Orientation!(Edges,Faces);
    Renumbering!(Edges,Faces);
  end  
  FacesInNodes!(Nodes,Faces)
  SortFacesInNodes!(Nodes,Faces)

  zP=zeros(nz)
  z=KernelAbstractions.zeros(backend,FT,nz+1)
  dzeta=zeros(nz)
  H=0.0
  nBar3 = zeros(0,0)
  NumNodesB = 0
  NumNodesG = 0
  NumEdgesB = 0
  NumEdgesG = 0
  NumFacesB = 0
  NumFacesG = 0
  AdaptGrid = ""
  EF=KernelAbstractions.zeros(backend,Int,0,0)
  FE=KernelAbstractions.zeros(backend,Int,0,0)
  return GridStruct{FT,
                    typeof(EF),
                    typeof(z)}(
    nz,
    zP,
    z,
    dzeta,
    H,
    NumFaces,
    NumFacesB,
    NumFacesG,
    Faces,
    NumEdges,
    NumEdgesB,
    NumEdgesG,
    Edges,
    NumNodes,
    NumNodesB,
    NumNodesG,
    Nodes,
    Form,
    Type,
    Dim,
    Rad,
    nBar3,
    nBar,
    AdaptGrid,
    EF,
    FE,
    )
end

function CubePoint(i1,i2,i3,n,x,y,z,Rad)

if i1>0
  x=tan(i1*pi/(2*n)-0.25*pi);
end
if i2>0
  y=tan(i2*pi/(2*n)-0.25*pi);
end
if i3>0
  z=tan(i3*pi/(2*n)-0.25*pi);
end
N=[x,y,z];
N=N/norm(N)*Rad;
return Point(N)
end

function InsertFaceEdge(n,EdgeNumber,NodeNumberStart,
                         NodeNumberE1Start1,NodeNumberE2Start1,
                         NodeNumberE1Start2,NodeNumberE2Start2,Edges,Nodes,Form,Rad)

  NodeNumber=NodeNumberStart;
  NodeNumberE1=NodeNumberE1Start1;
  NodeNumberE2=NodeNumberE2Start1;
  EdgeNumberStart1=EdgeNumber;
  @inbounds for j=1:n-1
    @inbounds for i=1:n
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
      Edges[EdgeNumber]=Edge([N1,N2],Nodes,EdgeNumber,EdgeNumber,"X",EdgeNumber;Form,Rad);
      EdgeNumber=EdgeNumber+1;
      NodeNumber=NodeNumber+1;
    end
  end
  NodeNumber=NodeNumberStart;
  NodeNumberE1=NodeNumberE1Start2;
  NodeNumberE2=NodeNumberE2Start2;
  EdgeNumberStart2=EdgeNumber;
  @inbounds for j=1:n
    @inbounds for i=1:n-1
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
      Edges[EdgeNumber]=Edge([N1,N2],Nodes,EdgeNumber,EdgeNumber,"X",EdgeNumber;Form,Rad);
      EdgeNumber=EdgeNumber+1;
      NodeNumber=NodeNumber+1;
    end
  end
  return (Edges,EdgeNumber,EdgeNumberStart1,EdgeNumberStart2)
end

function InsertEdgeEdge(n,EdgeNumber,
    NodeNumberStart,NodeNumberE1,NodeNumberE2,
     Edges,Nodes,Form,Rad)
  NodeNumber=NodeNumberStart;
  EdgeNumberStart=EdgeNumber;
  @inbounds for i=1:n
    N1=NodeNumber-1;
    N2=NodeNumber;
    if i==1
      N1=NodeNumberE1;
    end
    if i==n
      N2=NodeNumberE2;
    end
    Edges[EdgeNumber]=Edge([N1,N2],Nodes,EdgeNumber,EdgeNumber,"X",EdgeNumber;Form,Rad);
    EdgeNumber=EdgeNumber+1;
    NodeNumber=NodeNumber+1;
  end
  return (Edges,EdgeNumber,EdgeNumberStart)
end

function InsertFaceFace(n,FaceNumber,Type,EdgeNumberStart1,EdgeNumberStart2,
                         EdgeNumberStartEW1,EdgeNumberStartEW2,
                         EdgeNumberStartSN1,EdgeNumberStartSN2,
                         Nodes,Edges,Faces,OrientFace,Form,Rad)


  EdgeNumber1=EdgeNumberStart1;
  EdgeNumber2=EdgeNumberStart2;
  EdgeNumberEW1=EdgeNumberStartEW1;
  EdgeNumberEW2=EdgeNumberStartEW2;
  EdgeNumberSN1=EdgeNumberStartSN1;
  EdgeNumberSN2=EdgeNumberStartSN2;
  @inbounds for j=1:n
    @inbounds for i=1:n
      if (i==1 && j==1) || (i==1 && j==n)||
          (i==n && j==1) || (i==n && j==n)

      end
      #Face(FaceNumber)#Number=FaceNumber
      #Face(FaceNumber)#Boundary=1
      #ALLOCATE(Face(FaceNumber)#Edge(4))
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
      (Faces[FaceNumber],Edges)=Face([E1,E2,E3,E4],Nodes,Edges,FaceNumber,Type,OrientFace;Form=Form,Rad=Rad,P=zeros(Float64,0,0));
      FaceNumber=FaceNumber+1;
    end
  end
  return (Faces,FaceNumber,Edges)
end
