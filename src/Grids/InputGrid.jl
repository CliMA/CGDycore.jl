function TestOrientation(Faces,Nodes)
  NumFaces = length(Faces)
  for iF = 1 : NumFaces
    lenN = length(Faces[iF].N)
    Mid = Faces[iF].Mid
    for iN = 1 : lenN - 1
      P1 = Nodes[Faces[iF].N[iN]].P
      P2 = Nodes[Faces[iF].N[iN+1]].P
      O = dot(cross(P1 - Mid,P2 -Mid),Mid)
      if O <= 0
        @show iF,iN
      end
    end
    P1 = Nodes[Faces[iF].N[lenN]].P
    P2 = Nodes[Faces[iF].N[1]].P
    O = dot(cross(P1 - Mid,P2 -Mid),Mid)
    if O <= 0
      @show iF,lenN
    end
  end
  return nothing
end

function InputGridMPASO(backend,FT,filename,OrientFace,Rad,nz)
  xVertex = ncread(filename, "xVertex")
  yVertex = ncread(filename, "yVertex")
  zVertex = ncread(filename, "zVertex")

  verticesOnEdge = ncread(filename, "verticesOnEdge")
  edgesOnCell = ncread(filename, "edgesOnCell")
  nEdgesOnCell = ncread(filename, "nEdgesOnCell")
  nBar=[ 0  1   0   1
             -1  0  -1   0]
  Dim = 3
  Type = Polygonal()
  Form = "Sphere"
  NumNodes=size(xVertex,1)
  Nodes = map(1:NumNodes) do i
    Node()
  end
  NodeNumber = 1
  for i = 1 : NumNodes
    P = Point(xVertex[i],yVertex[i],zVertex[i]) 
    P = P / norm(P) * Rad
    Nodes[NodeNumber] = Node(P,NodeNumber,' ')
    NodeNumber += 1
  end
  NumEdges = size(verticesOnEdge,2)
  Edges = map(1:NumEdges) do i
    Edge()
  end
  EdgeNumber = 1
  for i = 1 : NumEdges
    Edges[EdgeNumber] = Edge(verticesOnEdge[:,i],Nodes,EdgeNumber,EdgeNumber,"",EdgeNumber)
    EdgeNumber += 1
  end
  NumFaces = size(edgesOnCell,2)
  Faces = map(1:NumFaces) do i
    Face()
  end
  FaceNumber = 1
  e = zeros(Int,size(edgesOnCell,1))
  for i = 1 : NumFaces
    e .= edgesOnCell[:,i]  
   (Faces[FaceNumber],Edges)=Face(e[1:nEdgesOnCell[i]],Nodes,Edges,FaceNumber,"",OrientFace;P=zeros(Float64,0,0));
    FaceNumber += 1
  end
  NumEdgesI = size(Edges,1)

  NumEdgesB = 0

  FacesInNodes!(Nodes,Faces)

  for iE = 1 : NumEdges
    if Edges[iE].F[2] == 0
      Edges[iE].Type = "B"
      Nodes[Edges[iE].N[1]].Type = 'B'
      Nodes[Edges[iE].N[2]].Type = 'B'
    end
  end  
  TestOrientation(Faces,Nodes)

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
    Nodes[NodeNumber] = Node(Point(N),NodeNumber,' ')
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
    N = sphereDeg2cart(Vertices[1,i],Vertices[2,i],Rad)
    Nodes[NodeNumber] = Node(Point(N),NodeNumber,' ')
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
  iLMeshFormat = 0  
  iLNodes = 0
  iLEntities = 0  
  iLElements = 0  
  iLPhysicalNames = 0  
  for iline = 1 : size(lines,1)
    lineLoc = lines[iline]  
    if lineLoc[2:end] == "MeshFormat"
      iLMeshFormat = iline  
    elseif lineLoc[2:end] == "Entities"
      iLEntities = iline  
    elseif lineLoc[2:end] == "Elements"
      iLElements = iline  
    elseif lineLoc[2:end] == "Nodes"
      iLNodes = iline  
    elseif lineLoc[2:end] == "PhysicalNames"
      iLPhysicalNames = iline  
    end
  end  
  s1,s2,s3,s4 = split(lines[iLNodes+1]," ")
  numEntityBlocks = parse(Int,s1)
  NumNodes = parse(Int,s2)
  minNodeTag = parse(Int,s3)
  maxNodeTag = parse(Int,s4)
  iline = iLNodes + 2
  NumNodes = 0
  NumNodesR = 0
  for iEntBl = 1 : numEntityBlocks
    s1,s2,s3,s4 = split(lines[iline]," ")
    entityDim = parse(Int,s1)
    entityTag = parse(Int,s2)
    parametric = parse(Int,s3)
    nodeTag = parse(Int,s4)
    iline = iline + nodeTag + 1
    for iNoTg = 1 : nodeTag
      if entityDim > 1
        NumNodes = NumNodes + 1  
      else  
        NumNodesR = NumNodesR + 1  
      end  
      iline = iline + 1
    end
  end

  Nodes = map(1:NumNodes) do i
    Node()
  end
  iline = iLNodes + 2
  NodeNumber = 1
  for iEntBl = 1 : numEntityBlocks
    s1,s2,s3,s4 = split(lines[iline]," ")
    entityDim = parse(Int,s1)
    entityTag = parse(Int,s2)
    parametric = parse(Int,s3)
    nodeTag = parse(Int,s4)
    iline = iline + nodeTag + 1
    for iNoTg = 1 : nodeTag
      if entityDim > 1  
        s1,s2,s3 = split(lines[iline]," ")
        r1 = parse(Float64,s1)
        r2 = parse(Float64,s2)
        r3 = parse(Float64,s3)
        N = [r1,r2,r3]
        N = N / norm(N) * Rad
        Nodes[NodeNumber] = Node(Point(N),NodeNumber,' ')
        NodeNumber = NodeNumber + 1
      end
      iline = iline + 1
    end  
  end    
  s1,s2,s3,s4 = split(lines[iLElements+1]," ")
  numEntityBlocks = parse(Int,s1)
  numElements = parse(Int,s2)
  minElementTag = parse(Int,s3)
  minElementTag = parse(Int,s4)

  iline = iLElements + 2
  NumFaces = 0
  for iEntBl = 1 : numEntityBlocks
    s1,s2,s3,s4 = split(lines[iline]," ")
    entityDim = parse(Int,s1)
    entityTag = parse(Int,s2)
    elementType = parse(Int,s3)
    numElementsInBlock = parse(Int,s4)
    iline = iline + 1
    for iElemTag = 1 : numElementsInBlock
      if elementType > 1
        NumFaces = NumFaces + 1  
      end
      iline = iline + 1
    end
  end  

  NumFacesL = 0
  NumFacesT = 0
  NumFacesQ = 0
  NodesF = zeros(Int,NumFaces,4)
  TypeF = zeros(Int,NumFaces)
  iF = 1
  EdgeList = Any[]
  iline = iLElements + 2
  for iEntBl = 1 : numEntityBlocks
    s1,s2,s3,s4 = split(lines[iline]," ")
    entityDim = parse(Int,s1)
    entityTag = parse(Int,s2)
    elementType = parse(Int,s3)
    numElementsInBlock = parse(Int,s4)  
    iline = iline + 1
    for iElemTag = 1 : numElementsInBlock
      if elementType == -1
        TypeF[iF] = 1  
        NumFacesL += 1
        s1,s2,s3 = split(lines[iline]," ")   
        NodesF[iF,1] = parse(Int,s2) - NumNodesR
        NodesF[iF,2] = parse(Int,s3) - NumNodesR
        if NodesF[iF,1] < NodesF[iF,2]
          e = [NodesF[iF,1],NodesF[iF,2]]
        else
          e = [NodesF[iF,2],NodesF[iF,1]]
        end
#       push!(EdgeList,e)
        iF += 1
      elseif elementType == 2
        TypeF[iF] = 2  
        NumFacesT += 1
        s1,s2,s3,s4 = split(lines[iline]," ")   
        NodesF[iF,1] = parse(Int,s2) - NumNodesR
        NodesF[iF,2] = parse(Int,s3) - NumNodesR
        NodesF[iF,3] = parse(Int,s4) - NumNodesR
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
        if NodesF[iF,3] < NodesF[iF,1]
          e = [NodesF[iF,3],NodesF[iF,1]]
        else
          e = [NodesF[iF,1],NodesF[iF,3]]
        end
        push!(EdgeList,e)
        iF += 1
      elseif elementType == 3
        NumFacesQ += 1
        TypeF[iF] = 3  
        s1,s2,s3,s4,s5 = split(lines[iline]," ")   
        NodesF[iF,1] = parse(Int,s2) - NumNodesR
        NodesF[iF,2] = parse(Int,s3) - NumNodesR
        NodesF[iF,3] = parse(Int,s4) - NumNodesR
        NodesF[iF,4] = parse(Int,s5) - NumNodesR
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
      iline += 1
    end
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

  NumFaces = NumFacesT + NumFacesQ
  Faces = map(1:NumFaces) do i
    Face()
  end
  FaceNumber = 1
  for iF = 1 : NumFaces
     if TypeF[iF] == 2
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
      if NodesF[iF,3] < NodesF[iF,1]
        e = [NodesF[iF,3],NodesF[iF,1]]
      else
        e = [NodesF[iF,1],NodesF[iF,3]]
      end
      e3 = EdgeDict[e]
      (Faces[FaceNumber],Edges)=Face([e1,e2,e3],Nodes,Edges,FaceNumber,"Tri",OrientFace;P=zeros(Float64,0,0));
      FaceNumber += 1
    elseif TypeF[iF] == 3   
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
      (Faces[FaceNumber],Edges)=Face([e1,e2,e3,e4],Nodes,Edges,FaceNumber,"Quad",OrientFace;P=zeros(Float64,0,0));
      FaceNumber += 1
    end
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
