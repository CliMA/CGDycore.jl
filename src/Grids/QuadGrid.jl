function QuadGrid(backend,FT,Fac)
  NumNodes = 4
  Rad = 1.0
  Form = "Planar"
  Dim = 3
  Type = Quad()
  Nodes = map(1:NumNodes) do i
    Node()
  end
  Nodes[1]=Node(Point([-Fac,-Fac,0.0]),1,' ')
  Nodes[2]=Node(Point([ Fac,-Fac,0.0]),2,' ')
  Nodes[3]=Node(Point([ Fac, Fac,0.0]),3,' ')
  Nodes[4]=Node(Point([-Fac, Fac,0.0]),4,' ')
  NumEdges = 4
  Edges= map(1:NumEdges) do i
    Edge([1,2],Nodes,0,0,"",0)
  end
  Edges[1] = Edge([1;2],Nodes,1,1,"",1;Form,Rad)
  Edges[2] = Edge([2;3],Nodes,2,2,"",2;Form,Rad)
  Edges[3] = Edge([3;4],Nodes,3,3,"",3;Form,Rad)
  Edges[4] = Edge([4;1],Nodes,4,4,"",4;Form,Rad)
  NumFaces = 1
  Faces = map(1:NumFaces) do i
    Face()
  end
  (Faces[1], Edges) = Face([1;2;3;4],Nodes,Edges,1,"Sphere",Grids.OrientFaceCart;
       P=zeros(Float64,0,0),Form=Form,Rad=Rad)
  NumNodes = size(Nodes,1)
  NumEdges = size(Edges,1)
  NumFaces = size(Faces,1)
  NumEdgesI = size(Edges,1)
  NumEdgesB = 0

  FacesInNodes!(Nodes,Faces)

  nz = 1
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
