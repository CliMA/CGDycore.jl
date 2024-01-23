function InputGrid(backend,FT,filename,OrientFace,Rad,nz)

  coord = ncread(filename, "coord")
  connect1 = ncread(filename, "connect1")

  nBar=[ 0  1   0   1
             -1  0  -1   0]
  Dim = 3
  Type = Quad()
  Rad = Rad
  Form = "Sphere"
  NumNodes=size(coord,1)
  Nodes = map(1:NumNodes) do i
    Node()
  end
  NodeNumber = 1

  for i = 1 : NumNodes
    N = coord[i,:]
    N = N / norm(N) * Rad
    Nodes[NodeNumber] = Node(Point(N),NodeNumber)
    NodeNumber = NodeNumber+1
  end

  NumFaces = size(connect1,2)
  EdgeNumber = 0
  EdgeList = Dict()
  for i = 1 : NumFaces
    n1 = connect1[1,i]
    n2 = connect1[2,i]
    if n1 < n2 
      EdgeNumber += 1
      EdgeList[(n1,n2)] = EdgeNumber
    end  
    n1 = connect1[2,i]
    n2 = connect1[3,i]
    if n1 < n2 
      EdgeNumber += 1
      EdgeList[(n1,n2)] = EdgeNumber
    end  
    n1 = connect1[3,i]
    n2 = connect1[4,i]
    if n1 < n2 
      EdgeNumber += 1
      EdgeList[(n1,n2)] = EdgeNumber
    end  
    n1 = connect1[4,i]
    n2 = connect1[1,i]
    if n1 < n2 
      EdgeNumber += 1
      EdgeList[(n1,n2)] = EdgeNumber
    end  
  end
  NumEdges = EdgeNumber
  Edges = map(1:NumEdges) do i
    Edge()
  end

  EdgeNumber = 0
  for i = 1 : NumFaces
    n1 = connect1[1,i]
    n2 = connect1[2,i]
    if n1 < n2
      EdgeNumber += 1
      Edges[EdgeNumber]=Edge([n1,n2],Nodes,EdgeNumber,EdgeNumber,"",EdgeNumber);
    end
    n1 = connect1[2,i]
    n2 = connect1[3,i]
    if n1 < n2
      EdgeNumber += 1
      Edges[EdgeNumber]=Edge([n1,n2],Nodes,EdgeNumber,EdgeNumber,"",EdgeNumber);
    end
    n1 = connect1[3,i]
    n2 = connect1[4,i]
    if n1 < n2
      EdgeNumber += 1
      Edges[EdgeNumber]=Edge([n1,n2],Nodes,EdgeNumber,EdgeNumber,"",EdgeNumber);
    end
    n1 = connect1[4,i]
    n2 = connect1[1,i]
    if n1 < n2
      EdgeNumber += 1
      Edges[EdgeNumber]=Edge([n1,n2],Nodes,EdgeNumber,EdgeNumber,"",EdgeNumber);
    end
  end  

  Faces = map(1:NumFaces) do i
    Face()
  end
  FaceNumber = 0
  for i = 1 : NumFaces
    n1 = connect1[1,i]
    n2 = connect1[2,i]
    if n1 < n2
      e1 = EdgeList[(n1,n2)]  
    else
      e1 = EdgeList[(n2,n1)]  
    end  
    n1 = connect1[2,i]
    n2 = connect1[3,i]
    if n1 < n2
      e2 = EdgeList[(n1,n2)]  
    else
      e2 = EdgeList[(n2,n1)]  
    end  
    n1 = connect1[3,i]
    n2 = connect1[4,i]
    if n1 < n2
      e3 = EdgeList[(n1,n2)]  
    else
      e3 = EdgeList[(n2,n1)]  
    end  
    n1 = connect1[4,i]
    n2 = connect1[1,i]
    if n1 < n2
      e4 = EdgeList[(n1,n2)]  
    else
      e4 = EdgeList[(n2,n1)]  
    end  
    FaceNumber += 1
    (Faces[FaceNumber],Edges)=Face([e1,e2,e3,e4],Nodes,Edges,FaceNumber,"",OrientFace;P=zeros(Float64,0,0));
  end

  NumNodes=size(Nodes,1);
  NumEdges=size(Edges,1);
  NumEdgesI=size(Edges,1);
  NumEdgesB=0;
  NumFaces=size(Faces,1);
  Dim=3;

  Orientation!(Edges,Faces);
  Renumbering!(Edges,Faces);
  FacesInNodes!(Nodes,Faces)

  zP=zeros(nz)
  z=KernelAbstractions.zeros(backend,FT,nz+1)
  dzeta=zeros(nz)
  H=0.0
  colors=[[]]
  NumGhostFaces = 0
  nBar3 = zeros(0,0)
  nBar = zeros(0,0)
  NumBoundaryFaces = 0
  AdaptGrid = ""

  return GridStruct{FT,
                    typeof(z)}(
    nz,
    zP,
    z,
    dzeta,
    H,
    NumFaces,
    NumGhostFaces,
    Faces,
    NumEdges,
    Edges,
    NumNodes,
    Nodes,
    Form,
    Type,
    Dim,
    Rad,
    NumEdgesI,
    NumEdgesB,
    nBar3,
    nBar,
    colors,
    NumBoundaryFaces,
    AdaptGrid,
    )
end

function InputGridH(backend,FT,filename,OrientFace,Rad,nz)

  Vertices = ncread(filename, "Vertices")
  NumNodes = size(Vertices,2)
  ListEdges = ncread(filename, "Edges")
  NumEdges = size(ListEdges,2)
  ListFaces = ncread(filename, "Cells")
  NumFaces = size(ListFaces,2)

  nBar=[ 0  1   0   1
             -1  0  -1   0]
  Dim = 3
  Type = "Quad"
  Rad = Rad
  Form = "Sphere"

  Nodes = map(1:NumNodes) do i
    Node()
  end
  NodeNumber = 1
  for i = 1 : NumNodes
    N = Grids.sphereDeg2cart(Vertices[1,i],Vertices[2,i],Rad)
    Nodes[NodeNumber] = Node(Point(N),NodeNumber)
    NodeNumber = NodeNumber+1
  end
  Nodes = Nodes

  Edges = map(1:NumEdges) do i
    Edge()
  end
  EdgeNumber = 1
  for i = 1 : NumEdges
    Edges[EdgeNumber] = Edge(ListEdges[:,i],Nodes,EdgeNumber,EdgeNumber,"",EdgeNumber)  
    EdgeNumber = EdgeNumber+1
  end

  Faces = map(1:NumFaces) do i
    Face()
  end
  FaceNumber = 0
  for i = 1 : NumFaces
    FaceNumber += 1
    e1=Int(ListFaces[1,i])
    e2=Int(ListFaces[2,i])
    e3=Int(ListFaces[3,i])
    e4=Int(ListFaces[4,i])
    (Faces[FaceNumber],Edges)=Face([e1,e2,e3,e4],Nodes,Edges,FaceNumber,"",OrientFace;P=zeros(Float64,0,0));
  end

  NumNodes=size(Nodes,1);
  NumEdges=size(Edges,1);
  NumEdgesI=size(Edges,1);
  NumEdgesB=0;
  NumFaces=size(Faces,1);
  Dim=3;
  Orientation!(Edges,Faces);
  Renumbering!(Edges,Faces);
  FacesInNodes!(Nodes,Faces)

  zP=zeros(nz)
  z=KernelAbstractions.zeros(backend,FT,nz+1)
  dzeta=zeros(nz)
  H=0.0
  colors=[[]]
  NumGhostFaces = 0
  nBar3 = zeros(0,0)
  nBar = zeros(0,0)
  NumBoundaryFaces = 0

  return GridStruct{FT,
                    typeof(z)}(
    nz,
    zP,
    z,
    dzeta,
    H,
    NumFaces,
    NumGhostFaces,
    Faces,
    NumEdges,
    Edges,
    NumNodes,
    Nodes,
    Form,
    Type,
    Dim,
    Rad,
    NumEdgesI,
    NumEdgesB,
    nBar3,
    nBar,
    colors,
    NumBoundaryFaces,
    )
end

function InputGridMsh(backend,FT,filename,OrientFace,Rad,nz)
  nBar=[ 0  1   0   1
             -1  0  -1   0]
  Dim = 3
  Type = Quad()
  Rad = Rad
  Form = "Sphere"
  f = open(filename)
  lines = readlines(f)
  close(f)
  s1,s2,s3,s4=split(lines[25]," ")
  NumNodes = parse(Int,s4)
  Nodes = map(1:NumNodes) do i
    Node()
  end
  NodeNumber = 1
  for i = 26 + NumNodes : 25 + 2 * NumNodes
    s1,s2,s3 = split(lines[i]," ")
    r1 = parse(Float64,s1)
    r2 = parse(Float64,s2)
    r3 = parse(Float64,s3)
    N = [r1,r2,r3]
    N = N / norm(N) * Rad
    Nodes[NodeNumber] = Node(Point(N),NodeNumber)
    NodeNumber = NodeNumber+1
  end
  s1,s2,s3,s4 = split(lines[34 + 2*NumNodes]," ")
  NumFaces = parse(Int,s4)
  NodesF = zeros(Int,NumFaces,4)
  iF = 1
  EdgeList = Any[]
  for i = 35 + 2*NumNodes : 34 + 2*NumNodes + NumFaces
    s1,s2,s3,s4,s5 = split(lines[i]," ")
    NodesF[iF,1] = parse(Int,s2) - 2
    NodesF[iF,2] = parse(Int,s3) - 2
    NodesF[iF,3] = parse(Int,s4) - 2
    NodesF[iF,4] = parse(Int,s5) - 2
    if NodesF[iF,1] < NodesF[iF,2]
      e = [NodesF[iF,1],NodesF[iF,2]]
    else
      e = [NodesF[iF,2],NodesF[iF,1]]
    end    
    push!(EdgeList,e)
    if NodesF[iF,2] < NodesF[iF,3]
      e = [NodesF[iF,2],NodesF[iF,3]]
    else
      e = [NodesF[iF,3],NodesF[iF,2]]
    end    
    push!(EdgeList,e)
    if NodesF[iF,3] < NodesF[iF,4]
      e = [NodesF[iF,3],NodesF[iF,4]]
    else
      e = [NodesF[iF,4],NodesF[iF,3]]
    end    
    push!(EdgeList,e)
    if NodesF[iF,4] < NodesF[iF,1]
      e = [NodesF[iF,4],NodesF[iF,1]]
    else
      e = [NodesF[iF,1],NodesF[iF,4]]
    end    
    push!(EdgeList,e)
    iF += 1
  end 
  EdgeList = unique(EdgeList)
  NumEdges = length(EdgeList)
  EdgeDict = Dict()
  iE = 1
  for i in EdgeList
    EdgeDict[i] = iE
    iE += 1
  end  
  Edges = map(1:NumEdges) do i
    Edge()
  end
  EdgeNumber = 1
  for i = 1 : NumEdges
    Edges[EdgeNumber] = Edge(EdgeList[i],Nodes,EdgeNumber,EdgeNumber,"",EdgeNumber)
    EdgeNumber += 1
  end

  Faces = map(1:NumFaces) do i
    Face()
  end
  FaceNumber = 1
  for iF = 1 : NumFaces
    if NodesF[iF,1] < NodesF[iF,2]
      e = [NodesF[iF,1],NodesF[iF,2]]
    else
      e = [NodesF[iF,2],NodesF[iF,1]]
    end
    e1 = EdgeDict[e]
    if NodesF[iF,2] < NodesF[iF,3]
      e = [NodesF[iF,2],NodesF[iF,3]]
    else
      e = [NodesF[iF,3],NodesF[iF,2]]
    end
    e2 = EdgeDict[e]
    if NodesF[iF,3] < NodesF[iF,4]
      e = [NodesF[iF,3],NodesF[iF,4]]
    else
      e = [NodesF[iF,4],NodesF[iF,3]]
    end
    e3 = EdgeDict[e]
    if NodesF[iF,4] < NodesF[iF,1]
      e = [NodesF[iF,4],NodesF[iF,1]]
    else
      e = [NodesF[iF,1],NodesF[iF,4]]
    end
    e4 = EdgeDict[e]
    (Faces[FaceNumber],Edges)=Face([e1,e2,e3,e4],Nodes,Edges,FaceNumber,"",OrientFace;P=zeros(Float64,0,0));
    FaceNumber += 1
  end 

  NumNodes=size(Nodes,1);
  NumEdges=size(Edges,1);
  NumEdgesI=size(Edges,1);
  NumEdgesB=0;
  NumFaces=size(Faces,1);
  Dim=3;
  Orientation!(Edges,Faces);
  Renumbering!(Edges,Faces);
  FacesInNodes!(Nodes,Faces)

  zP=zeros(nz)
  z=KernelAbstractions.zeros(backend,FT,nz+1)
  dzeta=zeros(nz)
  H=0.0
  colors=[[]]
  NumGhostFaces = 0
  nBar3 = zeros(0,0)
  nBar = zeros(0,0)
  NumBoundaryFaces = 0
  AdaptGrid = ""

  return GridStruct{FT,
                    typeof(z)}(
    nz,
    zP,
    z,
    dzeta,
    H,
    NumFaces,
    NumGhostFaces,
    Faces,
    NumEdges,
    Edges,
    NumNodes,
    Nodes,
    Form,
    Type,
    Dim,
    Rad,
    NumEdgesI,
    NumEdgesB,
    nBar3,
    nBar,
    colors,
    NumBoundaryFaces,
    AdaptGrid,
    )


end
