function Grid2KiteGrid(backend,FT,Grid,OrientFace)

  Type=Quad()
  Form = SphericalGrid()

  NumNodes = Grid.NumNodes + Grid.NumEdges + Grid.NumFaces
  Nodes = map(1:NumNodes) do i
    Node()
  end
  NodeNumber = 1
  for iN = 1 : Grid.NumNodes
     P = Grid.Nodes[iN].P 
     if Grid.Nodes[iN].Type == 'B'
       Nodes[NodeNumber]=Node(P,NodeNumber,'B') 
     else  
       Nodes[NodeNumber]=Node(P,NodeNumber,'N') 
     end    
     NodeNumber += 1
  end
  for iE = 1 : Grid.NumEdges
     P = Grid.Edges[iE].Mid
     Nodes[NodeNumber]=Node(P,NodeNumber,'E')
     NodeNumber += 1
  end
  for iF = 1 : Grid.NumFaces
     P = Grid.Faces[iF].Mid
     Nodes[NodeNumber]=Node(P,NodeNumber,'F')
     NodeNumber += 1
  end

  NumEdges = 2 * Grid.NumEdges
  for iF = 1 : Grid.NumFaces
    NumEdges += length(Grid.Faces[iF].E)
  end  
  Edges = map(1:NumEdges) do i
    Edge([1,2],Nodes,0,0,"",0);
  end
  EdgeNumber = 1
  for iE = 1 : Grid.NumEdges
    N1 = Grid.Edges[iE].N[1]  
    N2 = Grid.NumNodes + iE
    Edges[EdgeNumber] = Edge([N1,N2],Nodes,EdgeNumber,EdgeNumber,"E1",EdgeNumber)
    EdgeNumber += 1
    N1 = Grid.Edges[iE].N[2]  
    N2 = Grid.NumNodes + iE
    Edges[EdgeNumber] = Edge([N1,N2],Nodes,EdgeNumber,EdgeNumber,"E2",EdgeNumber)
    EdgeNumber += 1
  end  
  for iF = 1 : Grid.NumFaces
    for i = 1 : length(Grid.Faces[iF].E)  
      N1 = Grid.NumNodes + Grid.Faces[iF].E[i]
      N2 = Grid.NumNodes + Grid.NumEdges + iF
      Edges[EdgeNumber] = Edge([N1,N2],Nodes,EdgeNumber,EdgeNumber,"EM",EdgeNumber)
      EdgeNumber += 1
    end
  end  


  NumFaces = 0
  for iF = 1 : Grid.NumFaces
    NumFaces += length(Grid.Faces[iF].N)
  end   
  Faces = map(1:NumFaces) do i
    Face()
  end
  ie1 = 1
  FaceNumber = 1
  NumFaces = 0
  for iF = 1 : Grid.NumFaces
    for i = 1 : length(Grid.Faces[iF].N)  
      E2 = 2 * Grid.NumEdges + NumFaces + i
      iN = Grid.Faces[iF].N[i]  
      iEP = Grid.Faces[iF].E[i]
      if i == 1 
        E3 = 2 * Grid.NumEdges + NumFaces + length(Grid.Faces[iF].E)    
        iEM = Grid.Faces[iF].E[length(Grid.Faces[iF].N)]  
      else
        iEM = Grid.Faces[iF].E[i-1] 
        E3 = 2 * Grid.NumEdges + NumFaces + i - 1    
      end  
      if iN == Grid.Edges[iEP].N[1]
        E1 = 2 * iEP - 1
      else  
        E1 = 2 * iEP
      end  
      if iN == Grid.Edges[iEM].N[1]
        E4 = 2 * iEM - 1
      else  
        E4 = 2 * iEM
      end  
      (Faces[FaceNumber],Edges)=Face([E1,E2,E3,E4],Nodes,Edges,FaceNumber,"Quad",OrientFace;P=zeros(Float64,0,0))
      FaceNumber += 1
    end
    NumFaces += length(Grid.Faces[iF].E)
  end  

# Orientation!(Edges,Faces)
# Renumbering!(Edges,Faces)
  FacesInNodes!(Nodes,Faces)
  SortFacesInNodes!(Nodes,Faces)

  colors=[[]]
  NumGhostFaces = 0
  NumEdgesI = NumEdges
  NumEdgesB = 0
  NumBoundaryFaces = 0
  z = Grid.z
  NumNodesB = 0
  NumNodesG = 0
  NumEdgesB = 0
  NumEdgesG = 0
  NumFacesB = 0
  NumFacesG = 0
  EF=KernelAbstractions.zeros(backend,Int,0,0)
  FE=KernelAbstractions.zeros(backend,Int,0,0)
  return GridStruct{FT,
                    typeof(EF),
                    typeof(z)}(
    Grid.nz,
    Grid.zP,
    z,
    Grid.dzeta,
    Grid.H,
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
    Grid.Dim,
    Grid.Rad,
    Grid.nBar3,
    Grid.nBar,
    Grid.AdaptGrid,
    EF,
    FE,
    )
end

