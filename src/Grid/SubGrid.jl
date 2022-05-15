function ConstructSubGrid(GlobalGrid,Proc,ProcNumber)
  SubGrid = Grid(GlobalGrid.nz,GlobalGrid.Topography)

# Number of faces
  NumFaces = 0
  FaceNumbers = zeros(Int,0)
  EdgeNumbers = zeros(Int,0)
  NodeNumbers = zeros(Int,0)
  for iF = 1 : GlobalGrid.NumFaces
    if Proc[iF] == ProcNumber
      NumFaces += 1
      push!(FaceNumbers,iF)
      for i = 1 : 4
        push!(EdgeNumbers,GlobalGrid.Faces[iF].E[i])
        push!(NodeNumbers,GlobalGrid.Faces[iF].N[i])
      end
    end
  end
  EdgeNumbers = unique(EdgeNumbers)
  NodeNumbers = unique(NodeNumbers)

  NumNodes = size(NodeNumbers,1)
  Nodes = map(1:NumNodes) do i
    Node()
  end
  for i = 1:NumNodes
    Nodes[i] = deepcopy(GlobalGrid.Nodes[NodeNumbers[i]])  
    GlobalGrid.Nodes[NodeNumbers[i]].N = i
    Nodes[i].NG = Nodes[i].N
    Nodes[i].N = i
  end  
  @show Nodes[1]
  NumEdges = size(EdgeNumbers,1)
  Edges = map(1:NumEdges) do i
    Edge()
  end
  for i = 1:NumEdges
    Edges[i] = deepcopy(GlobalGrid.Edges[EdgeNumbers[i]])  
    GlobalGrid.Edges[EdgeNumbers[i]].E = i
    Edges[i].EG = Edges[i].E
    Edges[i].E = i
    Edges[i].N[1] = GlobalGrid.Nodes[Edges[i].N[1]].N
    Edges[i].N[2] = GlobalGrid.Nodes[Edges[i].N[2]].N
  end  
  @show Edges[1]

  Faces = map(1:NumFaces) do i
    Face()
  end
  for i = 1:NumFaces
    Faces[i] = deepcopy(GlobalGrid.Faces[FaceNumbers[i]])
    GlobalGrid.Faces[FaceNumbers[i]].F = i
    Faces[i].FG = Faces[i].F
    Faces[i].F = i
    Faces[i].E[1] = GlobalGrid.Edges[Faces[i].E[1]].E
    Faces[i].E[2] = GlobalGrid.Edges[Faces[i].E[2]].E
    Faces[i].E[3] = GlobalGrid.Edges[Faces[i].E[3]].E
    Faces[i].E[4] = GlobalGrid.Edges[Faces[i].E[4]].E
    Faces[i].N[1] = GlobalGrid.Nodes[Faces[i].N[1]].N
    Faces[i].N[2] = GlobalGrid.Nodes[Faces[i].N[2]].N
    Faces[i].N[3] = GlobalGrid.Nodes[Faces[i].N[3]].N
    Faces[i].N[4] = GlobalGrid.Nodes[Faces[i].N[4]].N
  end
  @show Faces[1]
  for i = 1:NumEdges
    Edges[i].F[1] = GlobalGrid.Faces[Edges[i].F[1]].F
    Edges[i].F[2] = GlobalGrid.Faces[Edges[i].F[2]].F
  end

  return SubGrid
end
