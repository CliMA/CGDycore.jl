function ConstructSubGrid(GlobalGrid,Proc,ProcNumber;order=true)
  backend = get_backend(GlobalGrid.z)
  FT = eltype(GlobalGrid.z)


# Number of faces
  DictF = Dict()
  NumFaces = 0
  FaceNumbersI = zeros(Int,0)
  FaceNumbersB = zeros(Int,0)
  GhostFaceNumbers = zeros(Int,0)
  EdgeNumbers = zeros(Int,0)
  NodeNumbers = zeros(Int,0)
  NumBoundaryFaces = 0
  NumInteriorFaces = 0
  @inbounds for iF = 1 : GlobalGrid.NumFaces
    if Proc[iF] == ProcNumber
      Inner = true  
      NumFaces += 1
      @inbounds for i = 1 : length(GlobalGrid.Faces[iF].N)
        push!(EdgeNumbers,GlobalGrid.Faces[iF].E[i])
        push!(NodeNumbers,GlobalGrid.Faces[iF].N[i])
        N = GlobalGrid.Nodes[GlobalGrid.Faces[iF].N[i]]
        for j = 1 : length(N.F)
          if Proc[N.F[j]] != ProcNumber
            Inner = false
          end
        end  
      end
      if Inner 
        NumInteriorFaces += 1  
        push!(FaceNumbersI,iF)
      else
        NumBoundaryFaces += 1  
        push!(FaceNumbersB,iF)
      end  
    end
  end
  FaceNumbers =[FaceNumbersB; FaceNumbersI]
  for i = 1 : length(FaceNumbers)
     iF = FaceNumbers[i]
     DictF[iF] = i
  end   
  EdgeNumbers = unique(EdgeNumbers)
  NodeNumbers = unique(NodeNumbers)

  NumNodes = size(NodeNumbers,1)
  Nodes = map(1:NumNodes) do i
    Node()
  end
  @inbounds for i = 1:NumNodes
    Nodes[i] = deepcopy(GlobalGrid.Nodes[NodeNumbers[i]])  
    GlobalGrid.Nodes[NodeNumbers[i]].N = i
    Nodes[i].NG = Nodes[i].N
    Nodes[i].FG = similar(Nodes[i].F)
    Nodes[i].FP = similar(Nodes[i].F)
    Nodes[i].FG .= Nodes[i].F
    Nodes[i].FP .= Proc[Nodes[i].F]
    Nodes[i].N = i
    Nodes[i].MasterSlave = 1
    @inbounds for j in eachindex(Nodes[i].FP)
      if ProcNumber > Nodes[i].FP[j]
        Nodes[i].MasterSlave = 0
        exit
      end
    end  
    @inbounds for j in eachindex(Nodes[i].FP)
      if ProcNumber != Nodes[i].FP[j]
        push!(GhostFaceNumbers,Nodes[i].FG[j])
      end  
    end  
  end  

  GhostFaceNumbers = unique(GhostFaceNumbers)
  i = 1
  NumGhostFaces = 0
  @inbounds for iGF in eachindex(GhostFaceNumbers)
    DictF[GhostFaceNumbers[iGF]] = i + NumFaces
    i += 1
    NumGhostFaces += 1
  end  

  NumEdges = size(EdgeNumbers,1)
  Edges = map(1:NumEdges) do i
    Edge()
  end
  @inbounds for i = 1:NumEdges
    Edges[i] = deepcopy(GlobalGrid.Edges[EdgeNumbers[i]])  
    GlobalGrid.Edges[EdgeNumbers[i]].E = i
    Edges[i].EG = Edges[i].E
    Edges[i].E = i
    Edges[i].N[1] = GlobalGrid.Nodes[Edges[i].N[1]].N
    Edges[i].N[2] = GlobalGrid.Nodes[Edges[i].N[2]].N
    Edges[i].FG .= Edges[i].F
    Edges[i].FP .= 0
    @inbounds for j = 1:2
      iF = Edges[i].F[j]
      if iF > 0
        Edges[i].FP[j] = Proc[iF]
      end
    end  
    Edges[i].MasterSlave = 1
    @inbounds for j in eachindex(Edges[i].FP)
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
  @inbounds for i = 1:NumFaces
    Faces[i] = deepcopy(GlobalGrid.Faces[FaceNumbers[i]])
    GlobalGrid.Faces[FaceNumbers[i]].F = i
    Faces[i].FG = Faces[i].F
    Faces[i].F = i
    for j = 1 : length(Faces[i].E)
      Faces[i].E[j] = GlobalGrid.Edges[Faces[i].E[j]].E
      Faces[i].N[j] = GlobalGrid.Nodes[Faces[i].N[j]].N
    end
  end
  # Physischer Rand und Prozessor Rand 
  @inbounds for i = 1:NumEdges
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

  Dim=3;
  if GlobalGrid.Type == Quad()
    if order  
      Renumbering!(Edges,Faces);
    end  
  end
  FacesInNodes!(Nodes,Faces)
  Form = GlobalGrid.Form
  Rad = GlobalGrid.Rad
# Stencil  
  @inbounds for iF=1:NumFaces
    StencilLoc=zeros(Int, 16,1);
    StencilLoc[:] .= iF;
    iS=0;
    @inbounds for i=1:length(Faces[iF].N)
      iN=Faces[iF].N[i];
      @inbounds for j=1:size(Nodes[iN].F,1)
        jF=Nodes[iN].F[j];
        inside=false;
        @inbounds for jS=1:iS
          if StencilLoc[jS]==jF
            inside=true;
            break
          end
        end
        if !inside
          iS=iS+1;
          StencilLoc[iS]=jF;
        end
      end
    end
    Faces[iF].Stencil = zeros(Int,iS)
    Faces[iF].Stencil .= StencilLoc[1:iS]
  end

# Add Ghost Faces
  @inbounds for i = 1:NumNodes
    Nodes[i].F = zeros(Int,0)
    @inbounds for j in eachindex(Nodes[i].FG)
      push!(Nodes[i].F,DictF[Nodes[i].FG[j]])  
    end  
  end  
  nz = GlobalGrid.nz
  zP=zeros(nz)
  z=KernelAbstractions.zeros(backend,FT,nz+1)
  dzeta=zeros(nz)
  H = GlobalGrid.H
  NumEdgesI = 0
  NumEdgesB = 0
  for iE = 1 : NumEdges
    if Edges[iE].F[2] == 0
      if ProcNumber == 2 && Edges[iE].Type == "B" 
        NumEdgesB += 1
      elseif ProcNumber == 2 && Edges[iE].Type == ""
        NumEdgesI += 1
      end
    else
      NumEdgesI += 1  
    end  
  end    
  colors=[[]]
  nBar3=zeros(0,0)
  nBar=zeros(0,0)
  Type = GlobalGrid.Type
  AdaptGrid = GlobalGrid.AdaptGrid
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

function Decompose(Grid,NumProc)

  NumFaces=Grid.NumFaces
  LocalNumfaces=zeros(Int,NumProc)
  LocalNumfaces.= floor(NumFaces / NumProc)
  Rest = mod(NumFaces, NumProc)
  @inbounds for iP=1:Rest
    LocalNumfaces[iP]+=1
  end
  CellToProc = zeros(Int,NumFaces)
  @inbounds for ic = 1 : NumProc
    CellToProc[sum(LocalNumfaces[1:ic-1])+1:sum(LocalNumfaces[1:ic])] .= ic
  end  
  return CellToProc
end

function DecomposeEqualArea(Grid,NumProc)

  NumFaces=Grid.NumFaces
  LocalNumfaces=zeros(Int,NumProc)
  LocalNumfaces.= floor(NumFaces / NumProc)
  Rest = mod(NumFaces, NumProc)
  @inbounds for iP=1:Rest
    LocalNumfaces[iP]+=1
  end
  (n_regions,s_cap) = Parallels.eq_caps(NumProc)
  n_collars = size(n_regions,1) - 2
  CellToProc = zeros(Int,NumFaces)
  Faces=Grid.Faces
  coord=zeros(Float64,NumFaces,3)
  @inbounds for iF=1:NumFaces
    P = Faces[iF].Mid  
    coord[iF,3]=iF
    (coord[iF,1],coord[iF,2],r) = cart2sphere(P.x,P.y,P.z)  
  end  
  p=sortslices(coord, dims=1, lt=Parallels.compare_NS_WE)
  NumNodeA = 1
  NumNodeE = LocalNumfaces[1]
  iP = zeros(Int,NumNodeE-NumNodeA+1)
  @. iP = round(Int,p[NumNodeA:NumNodeE,3])
  CellToProc[iP] .= 1
  if NumProc > 1
    region_n = 2
    NumNodeECol = NumNodeE   
    @inbounds for collar_n = 1 : n_collars
      NumNodeACol = NumNodeECol + 1
      region_nCol = region_n
      @inbounds for region_ew = 1 : n_regions[collar_n + 1]
        NumNodeECol += LocalNumfaces[region_nCol]  
        region_nCol += 1
      end
      pCol = sortslices(p[NumNodeACol:NumNodeECol,:], dims=1, lt=Parallels.compare_WE_NS)
      NumNodeE = 0
      @inbounds for region_ew = 1 : n_regions[collar_n + 1]
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
