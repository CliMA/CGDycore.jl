function Grid2KiteGrid(backend,FT,Grid,OrientFace)

  Type=Quad()

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
    N1 = Grid.NumNodes + iE
    N2 = Grid.Edges[iE].N[2]  
    Edges[EdgeNumber] = Edge([N1,N2],Nodes,EdgeNumber,EdgeNumber,"E2",EdgeNumber)
    EdgeNumber += 1
  end  
  for iF = 1 : Grid.NumFaces
    for i = 1 : length(Grid.Faces[iF].E)  
      N1 = Grid.NumNodes + Grid.NumEdges + iF
      N2 = Grid.NumNodes + Grid.Faces[iF].E[i]
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
  for iF = 1 : Grid.NumFaces
    for i = 1 : length(Grid.Faces[iF].E)  
       E4 = ie1 + 2 * Grid.NumEdges 
       if i == 1
         E1 = ie1 + 2 * Grid.NumEdges + length(Grid.Faces[iF].N) - 1
         if Grid.Faces[iF].OrientE[end] * Grid.Faces[iF].Orientation == 1
           E2 = 2 * Grid.Faces[iF].E[end] 
         else
           E2 = 2 * Grid.Faces[iF].E[end] - 1  
         end  
       else
         E1 = E4 - 1
         if Grid.Faces[iF].OrientE[i-1] * Grid.Faces[iF].Orientation == 1
           E2 = 2 * Grid.Faces[iF].E[i-1] 
         else
           E2 = 2 * Grid.Faces[iF].E[i-1] - 1 
         end  
       end
       if Grid.Faces[iF].OrientE[i] * Grid.Faces[iF].Orientation == 1
         E3 = 2 * Grid.Faces[iF].E[i] - 1
       else  
         E3 = 2 * Grid.Faces[iF].E[i]
       end  
       (Faces[FaceNumber],Edges)=Face([E1,E2,E3,E4],Nodes,Edges,FaceNumber,"Quad",OrientFace;P=zeros(Float64,0,0))
       FaceNumber += 1
       ie1 += 1
    end
  end  

  Orientation!(Edges,Faces)
# Renumbering!(Edges,Faces)
  FacesInNodes!(Nodes,Faces)
  SortFacesInNodes!(Nodes,Faces)

  colors=[[]]
  NumGhostFaces = 0
  NumEdgesI = NumEdges
  NumEdgesB = 0
  NumBoundaryFaces = 0
  z = Grid.z
  return GridStruct{FT,
                    typeof(z)}(
    Grid.nz,
    Grid.zP,
    z,
    Grid.dzeta,
    Grid.H,
    NumFaces,
    NumGhostFaces,
    Faces,
    NumEdges,
    Edges,
    NumNodes,
    Nodes,
    Grid.Form,
    Type,
    Grid.Dim,
    Grid.Rad,
    NumEdgesI,
    NumEdgesB,
    Grid.nBar3,
    Grid.nBar,
    colors,
    NumBoundaryFaces,
    Grid.AdaptGrid,
    )
end

