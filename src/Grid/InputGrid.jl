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
