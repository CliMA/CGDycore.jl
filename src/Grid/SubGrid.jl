function ConstructSubGrid(GlobalGrid,Proc,ProcNumber)
  SubGrid = Grid(GlobalGrid.nz,GlobalGrid.Topography)

# Number of faces
  DictF = Dict()
  NumFaces = 0
  FaceNumbers = zeros(Int,0)
  GhostFaceNumbers = zeros(Int,0)
  EdgeNumbers = zeros(Int,0)
  NodeNumbers = zeros(Int,0)
  for iF = 1 : GlobalGrid.NumFaces
    if Proc[iF] == ProcNumber
      NumFaces += 1
      DictF[iF] = NumFaces 
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
    Nodes[i].MasterSlave = 1
    for j in eachindex(Nodes[i].FP)
      if ProcNumber > Nodes[i].FP[j]
        Nodes[i].MasterSlave = 0
        exit
      end
    end  
    for j in eachindex(Nodes[i].FP)
      if ProcNumber != Nodes[i].FP[j]
        push!(GhostFaceNumbers,Nodes[i].FG[j])
      end  
    end  
  end  

  GhostFaceNumbers = unique(GhostFaceNumbers)
  i = 1
  NumGhostFaces = 0
  for iGF in eachindex(GhostFaceNumbers)
    DictF[GhostFaceNumbers[iGF]] = i + NumFaces
    i += 1
    NumGhostFaces += 1
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
    Edges[i].FP .= 0
    for j = 1:2
      iF = Edges[i].F[j]
      if iF > 0
        Edges[i].FP[j] = Proc[iF]
      end
    end  
    Edges[i].MasterSlave = 1
    for j in eachindex(Edges[i].FP)
      if Edges[i].FP[j] >0  
        if ProcNumber > Edges[i].FP[j] 
          Edges[i].MasterSlave = 0
          exit
        end
      end
    end  
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
  for i = 1:NumEdges
    if Edges[i].F[1] > 0  
      if Proc[Edges[i].F[1]] == ProcNumber  
        Edges[i].F[1] = GlobalGrid.Faces[Edges[i].F[1]].F
      else
        Edges[i].F[1] = 0  
      end  
    end
    if Edges[i].F[2] > 0  
      if Proc[Edges[i].F[2]] == ProcNumber  
        Edges[i].F[2] = GlobalGrid.Faces[Edges[i].F[2]].F
      else
        Edges[i].F[2] = 0  
      end  
    end
  end

  SubGrid.NumFaces = NumFaces
  SubGrid.NumGhostFaces = NumGhostFaces
  SubGrid.Faces = Faces
  SubGrid.NumEdges = NumEdges
  SubGrid.Edges = Edges
  SubGrid.NumNodes = NumNodes
  SubGrid.Nodes = Nodes

  SubGrid.Dim=3;
  SubGrid=Renumbering(SubGrid);
  SubGrid=FacesInNodes(SubGrid);
  SubGrid.Form = GlobalGrid.Form

  #Boundary/Interior faces
  BoundaryFaces = zeros(Int,0)
  for iE = 1 : SubGrid.NumEdges
    if SubGrid.Edges[iE].F[1] == 0 || SubGrid.Edges[iE].F[2] == 0
      for iN in SubGrid.Edges[iE].N
        for iF in SubGrid.Nodes[iN].F  
          push!(BoundaryFaces,iF)
        end
      end
    end
  end  
  BoundaryFaces = unique(BoundaryFaces)
  SubGrid.BoundaryFaces = BoundaryFaces
  SubGrid.InteriorFaces = setdiff(collect(UnitRange(1,SubGrid.NumFaces)),SubGrid.BoundaryFaces)

  SubGrid.Rad = GlobalGrid.Rad
# Add Ghost Faces
  for i = 1:NumNodes
    Nodes[i].F = zeros(Int,0)
    for j in eachindex(Nodes[i].FG)
      push!(Nodes[i].F,DictF[Nodes[i].FG[j]])  
    end  
  end  


  return SubGrid
end

function Decompose(Grid,NumProc)

  NumFaces=Grid.NumFaces
  LocalNumfaces=zeros(Int,NumProc)
  LocalNumfaces.= floor(NumFaces / NumProc)
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

function DecomposeEqualArea(Grid,NumProc)

  NumFaces=Grid.NumFaces
  LocalNumfaces=zeros(Int,NumProc)
  LocalNumfaces.= floor(NumFaces / NumProc)
  Rest = mod(NumFaces, NumProc)
  for iP=1:Rest
    LocalNumfaces[iP]+=1
  end
  (n_regions,s_cap) = eq_caps(NumProc)
  n_collars = size(n_regions,1) - 2
  CellToProc = zeros(Int,NumFaces)
  Faces=Grid.Faces
  coord=zeros(Float64,NumFaces,3)
  for iF=1:NumFaces
    P = Faces[iF].Mid  
    coord[iF,3]=iF
    (coord[iF,1],coord[iF,2],r) = cart2sphere(P.x,P.y,P.z)  
  end  
  p=sortslices(coord, dims=1, lt=compare_NS_WE)
  NumNodeA = 1
  NumNodeE = LocalNumfaces[1]
  iP = zeros(Int,NumNodeE-NumNodeA+1)
  @. iP = round(Int,p[NumNodeA:NumNodeE,3])
  CellToProc[iP] .= 1
  if NumProc > 1
    region_n = 2
    NumNodeECol = NumNodeE   
    for collar_n = 1 : n_collars
      NumNodeACol = NumNodeECol + 1
      region_nCol = region_n
      for region_ew = 1 : n_regions[collar_n + 1]
        NumNodeECol += LocalNumfaces[region_nCol]  
        region_nCol += 1
      end
      pCol = sortslices(p[NumNodeACol:NumNodeECol,:], dims=1, lt=compare_WE_NS)
      NumNodeE = 0
      for region_ew = 1 : n_regions[collar_n + 1]
        NumNodeA = NumNodeE + 1  
        NumNodeE +=  LocalNumfaces[region_n] 
        iP = zeros(Int,NumNodeE-NumNodeA+1)
        @. iP = round(Int,pCol[NumNodeA:NumNodeE,3])
        CellToProc[iP] .= region_n
        region_n += 1
      end
    end
    NumNodeA = NumNodeECol + 1
    NumNodeE =  NumNodeECol + LocalNumfaces[end]
    iP = zeros(Int,NumNodeE-NumNodeA+1)
    @. iP = round(Int,p[NumNodeA:NumNodeE,3])
    CellToProc[iP] .= NumProc
  end   
  return CellToProc
end
