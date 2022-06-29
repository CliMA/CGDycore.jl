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
    Nodes[i].FG = similar(Nodes[i].F)
    Nodes[i].FP = similar(Nodes[i].F)
    Nodes[i].FG .= Nodes[i].F
    Nodes[i].FP .= Proc[Nodes[i].F]
    Nodes[i].N = i
  end  
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
    Edges[i].FG .= Edges[i].F
    Edges[i].FP .= Proc[Edges[i].F]
  end  

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
  # Physischer Rand und Prozessor Rand 
  NumInBoundEdges = 0
  for i = 1:NumEdges
    if Proc[Edges[i].F[1]] == ProcNumber  
      Edges[i].F[1] = GlobalGrid.Faces[Edges[i].F[1]].F
    else
      Edges[i].F[1] = 0  
      NumInBoundEdges += 1
    end  
    if Proc[Edges[i].F[2]] == ProcNumber  
      Edges[i].F[2] = GlobalGrid.Faces[Edges[i].F[2]].F
    else
      Edges[i].F[2] = 0  
      NumInBoundEdges += 1
    end  
  end
  InBoundEdges = zeros(Int,NumInBoundEdges)
  InBoundEdgesP = zeros(Int,NumInBoundEdges)
  NumInBoundEdges = 0
  for i = 1:NumEdges
    if Proc[Edges[i].FG[1]] == ProcNumber  
    else
      NumInBoundEdges += 1
      InBoundEdges[NumInBoundEdges] = i
      InBoundEdgesP[NumInBoundEdges] = Proc[Edges[i].FG[1]]
    end  
    if Proc[Edges[i].FG[2]] == ProcNumber  
    else
      NumInBoundEdges += 1
      InBoundEdges[NumInBoundEdges] = i
      InBoundEdgesP[NumInBoundEdges] = Proc[Edges[i].FG[2]]
    end  
  end
  SubGrid.NumNeiProc = size(unique(InBoundEdgesP),1)
  SubGrid.NeiProc = unique(InBoundEdgesP)

  SubGrid.NumFaces = NumFaces
  SubGrid.Faces = Faces
  SubGrid.NumEdges = NumEdges
  SubGrid.Edges = Edges
  SubGrid.NumNodes = NumNodes
  SubGrid.Nodes = Nodes
  SubGrid.NumInBoundEdges = NumInBoundEdges
  SubGrid.InBoundEdges = InBoundEdges
  SubGrid.InBoundEdgesP = InBoundEdgesP

  SubGrid.Dim=3;
  SubGrid=Renumbering(SubGrid);
  SubGrid=FacesInNodes(SubGrid);
  SubGrid.Form = GlobalGrid.Form

  return SubGrid
end

function Decompose(Grid,NumProc)

  NumFaces=Grid.NumFaces
  LocalNumfaces=zeros(Int,NumProc)
  LocalNumfaces.=NumFaces / NumProc
  Rest = mod(NumFaces, NumProc)
  for iP=1:Rest
    LocalNumfaces[iP]+=1
  end
  CellToProc = zeros(Int,NumFaces)
  for ic = 1 : NumProc
    CellToProc[sum(LocalNumfaces[1:ic-1])+1:sum(LocalNumfaces[1:ic])] .= ic
  end  
  return CellToProc
end
