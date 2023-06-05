function InputGrid(filename,OrientFace,Rad,Grid)

  coord = ncread(filename, "coord")
  connect1 = ncread(filename, "connect1")

  Grid.nBar=[ 0  1   0   1
             -1  0  -1   0]
  Grid.Dim = 3
  Grid.Type = "Quad"
  Grid.Rad = Rad
  Grid.Form = "Sphere"
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

  Grid.Nodes = Nodes

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
    Edge([1,2],Grid,0,0,"",0)
  end

  EdgeNumber = 0
  for i = 1 : NumFaces
    n1 = connect1[1,i]
    n2 = connect1[2,i]
    if n1 < n2
      EdgeNumber += 1
      Edges[EdgeNumber]=Edge([n1,n2],Grid,EdgeNumber,EdgeNumber,"",EdgeNumber);
    end
    n1 = connect1[2,i]
    n2 = connect1[3,i]
    if n1 < n2
      EdgeNumber += 1
      Edges[EdgeNumber]=Edge([n1,n2],Grid,EdgeNumber,EdgeNumber,"",EdgeNumber);
    end
    n1 = connect1[3,i]
    n2 = connect1[4,i]
    if n1 < n2
      EdgeNumber += 1
      Edges[EdgeNumber]=Edge([n1,n2],Grid,EdgeNumber,EdgeNumber,"",EdgeNumber);
    end
    n1 = connect1[4,i]
    n2 = connect1[1,i]
    if n1 < n2
      EdgeNumber += 1
      Edges[EdgeNumber]=Edge([n1,n2],Grid,EdgeNumber,EdgeNumber,"",EdgeNumber);
    end
  end  
  Grid.Edges = Edges

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
    (Faces[FaceNumber],Grid)=Face([e1,e2,e3,e4],Grid,FaceNumber,"",OrientFace;P=zeros(Float64,0,0));
  end
  Grid.Faces = Faces

  Grid.NumNodes=size(Grid.Nodes,1);
  Grid.NumEdges=size(Grid.Edges,1);
  Grid.NumEdgesI=size(Grid.Edges,1);
  Grid.NumEdgesB=0;
  Grid.NumFaces=size(Grid.Faces,1);
  Grid.Dim=3;
  Grid=Orientation(Grid);
  Grid=Renumbering(Grid);
  Grid=FacesInNodes(Grid);

  #Boundary/Interior faces
  BoundaryFaces = zeros(Int,0)
  for iE = 1 : Grid.NumEdges
    if Grid.Edges[iE].F[1] == 0 || Grid.Edges[iE].F[2] == 0
      for iN in Grid.Edges[iE].N
        for iF in Grid.Nodes[iN].F
          push!(BoundaryFaces,iF)
        end
      end
    end
  end
  BoundaryFaces = unique(BoundaryFaces)
  Grid.BoundaryFaces = BoundaryFaces
  Grid.InteriorFaces = setdiff(collect(UnitRange(1,Grid.NumFaces)),Grid.BoundaryFaces)
  return Grid
end

function InputGridH(filename,OrientFace,Rad,Grid)

  Vertices = ncread(filename, "Vertices")
  NumNodes = size(Vertices,2)
  ListEdges = ncread(filename, "Edges")
  NumEdges = size(ListEdges,2)
  ListFaces = ncread(filename, "Cells")
  NumFaces = size(ListFaces,2)

  Grid.nBar=[ 0  1   0   1
             -1  0  -1   0]
  Grid.Dim = 3
  Grid.Type = "Quad"
  Grid.Rad = Rad
  Grid.Form = "Sphere"

  Nodes = map(1:NumNodes) do i
    Node()
  end
  NodeNumber = 1
  for i = 1 : NumNodes
    N = sphereDeg2cart(Vertices[1,i],Vertices[2,i],Rad)
    Nodes[NodeNumber] = Node(Point(N),NodeNumber)
    NodeNumber = NodeNumber+1
  end
  Grid.Nodes = Nodes

  Edges = map(1:NumEdges) do i
    Edge()
  end
  EdgeNumber = 1
  for i = 1 : NumEdges
    Edges[EdgeNumber] = Edge(ListEdges[:,i],Grid,EdgeNumber,EdgeNumber,"",EdgeNumber)  
    EdgeNumber = EdgeNumber+1
  end
  Grid.Edges = Edges

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
    (Faces[FaceNumber],Grid)=Face([e1,e2,e3,e4],Grid,FaceNumber,"",OrientFace;P=zeros(Float64,0,0));
  end
  Grid.Faces = Faces

  Grid.NumNodes=size(Grid.Nodes,1);
  Grid.NumEdges=size(Grid.Edges,1);
  Grid.NumEdgesI=size(Grid.Edges,1);
  Grid.NumEdgesB=0;
  Grid.NumFaces=size(Grid.Faces,1);
  Grid.Dim=3;
  Grid=Orientation(Grid);
  Grid=Renumbering(Grid);
  Grid=FacesInNodes(Grid);

  #Boundary/Interior faces
  BoundaryFaces = zeros(Int,0)
  for iE = 1 : Grid.NumEdges
    if Grid.Edges[iE].F[1] == 0 || Grid.Edges[iE].F[2] == 0
      for iN in Grid.Edges[iE].N
        for iF in Grid.Nodes[iN].F
          push!(BoundaryFaces,iF)
        end
      end
    end
  end
  BoundaryFaces = unique(BoundaryFaces)
  Grid.BoundaryFaces = BoundaryFaces
  Grid.InteriorFaces = setdiff(collect(UnitRange(1,Grid.NumFaces)),Grid.BoundaryFaces)
  return Grid
end

function InputGridMsh(filename,OrientFace,Rad,Grid)
  Grid.nBar=[ 0  1   0   1
             -1  0  -1   0]
  Grid.Dim = 3
  Grid.Type = "Quad"
  Grid.Rad = Rad
  Grid.Form = "Sphere"
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
  Grid.Nodes = Nodes
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
    Edges[EdgeNumber] = Edge(EdgeList[i],Grid,EdgeNumber,EdgeNumber,"",EdgeNumber)
    EdgeNumber += 1
  end
  Grid.Edges = Edges

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
    (Faces[FaceNumber],Grid)=Face([e1,e2,e3,e4],Grid,FaceNumber,"",OrientFace;P=zeros(Float64,0,0));
    FaceNumber += 1
  end 

  Grid.Faces = Faces

  Grid.NumNodes=size(Grid.Nodes,1);
  Grid.NumEdges=size(Grid.Edges,1);
  Grid.NumEdgesI=size(Grid.Edges,1);
  Grid.NumEdgesB=0;
  Grid.NumFaces=size(Grid.Faces,1);
  Grid.Dim=3;
  Grid=Orientation(Grid);
  Grid=Renumbering(Grid);
  Grid=FacesInNodes(Grid);

  #Boundary/Interior faces
  BoundaryFaces = zeros(Int,0)
  for iE = 1 : Grid.NumEdges
    if Grid.Edges[iE].F[1] == 0 || Grid.Edges[iE].F[2] == 0
      for iN in Grid.Edges[iE].N
        for iF in Grid.Nodes[iN].F
          push!(BoundaryFaces,iF)
        end
      end
    end
  end
  BoundaryFaces = unique(BoundaryFaces)
  Grid.BoundaryFaces = BoundaryFaces
  Grid.InteriorFaces = setdiff(collect(UnitRange(1,Grid.NumFaces)),Grid.BoundaryFaces)
  return Grid

end
