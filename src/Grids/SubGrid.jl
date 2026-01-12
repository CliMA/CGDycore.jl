function Decompose(Grid,nx,ny,NumProc)

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
  LenLocalX = div(nx,NumProc) 
  RestX = mod(nx,NumProc)
  iProc1 = 1
  iProc2 = LenLocalX + 1
  Proc = 0
  for i = 1 : RestX
    Proc  += 1  
    for iy = 1 : ny
      @. CellToProc[iProc1+(iy-1)*nx:iProc2+(iy-1)*nx] = Proc  
    end
    iProc1 += LenLocalX + 1
    iProc2 += LenLocalX + 1
  end  
  iProc2 -= 1
  for i = RestX + 1 : NumProc
    Proc  += 1  
    for iy = 1 : ny
      @. CellToProc[iProc1+(iy-1)*nx:iProc2+(iy-1)*nx] = Proc  
    end
    iProc1 += LenLocalX 
    iProc2 += LenLocalX
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

function ConstructSubGridGhost(GlobalGrid,Proc,ProcNumber;order=true)
  backend = get_backend(GlobalGrid.z)
  FT = eltype(GlobalGrid.z)


# Number of faces
  DictF = Dict()
  NumFaces = 0
  FaceNumbersI = zeros(Int,0)
  FaceNumbersB = zeros(Int,0)
  FaceNumbersG = zeros(Int,0)
  EdgeNumbers = zeros(Int,0)
  EdgeNumbersG = zeros(Int,0)
  NodeNumbers = zeros(Int,0)
  NodeNumbersG = zeros(Int,0)

  NumFacesB = 0
  NumFacesI = 0
  @inbounds for iF = 1 : GlobalGrid.NumFaces
    if Proc[iF] == ProcNumber
      Inner = true  
      NumFaces += 1
      for i = 1 : length(GlobalGrid.Faces[iF].N)
        push!(EdgeNumbers,GlobalGrid.Faces[iF].E[i])
        push!(NodeNumbers,GlobalGrid.Faces[iF].N[i])
        N = GlobalGrid.Nodes[GlobalGrid.Faces[iF].N[i]]
        for j = 1 : length(N.F)
          if Proc[N.F[j]] != ProcNumber
            Inner = false
            iFG = N.F[j]
            push!(FaceNumbersG,N.F[j])
            for k = 1 : length(GlobalGrid.Faces[iFG].N)
              push!(EdgeNumbersG,GlobalGrid.Faces[iFG].E[k])
              push!(NodeNumbersG,GlobalGrid.Faces[iFG].N[k])
            end  
          end
        end  
      end
      if Inner 
        NumFacesI += 1  
        push!(FaceNumbersI,iF)
      else
        NumFacesB += 1  
        push!(FaceNumbersB,iF)
      end  
    end
  end
  FaceNumbersG = unique(FaceNumbersG)
  NumFacesG = length(FaceNumbersG)

  NodeNumbers = unique(NodeNumbers)
  NodeNumbersG = unique(NodeNumbersG)
  NodeNumbersI = setdiff(NodeNumbers,NodeNumbersG)
  NodeNumbersB = Base.intersect(NodeNumbers,NodeNumbersG)
  NodeNumbersG = setdiff(NodeNumbersG,NodeNumbers)
  NumNodesI = length(NodeNumbersI)
  NumNodesB = length(NodeNumbersB)
  NumNodesG = length(NodeNumbersG)
  NumNodes = NumNodesI + NumNodesB

  EdgeNumbers = unique(EdgeNumbers)
  EdgeNumbersG = unique(EdgeNumbersG)
  EdgeNumbersI = setdiff(EdgeNumbers,EdgeNumbersG)
  EdgeNumbersB = Base.intersect(EdgeNumbers,EdgeNumbersG)
  EdgeNumbersG = setdiff(EdgeNumbersG,EdgeNumbers)
  NumEdgesI = length(EdgeNumbersI)
  NumEdgesB = length(EdgeNumbersB)
  NumEdgesG = length(EdgeNumbersG)
  NumEdges = NumEdgesI + NumEdgesB

  FaceNumbers =[FaceNumbersB; FaceNumbersI; FaceNumbersG]
  EdgeNumbers =[EdgeNumbersB; EdgeNumbersI; EdgeNumbersG]
  NodeNumbers =[NodeNumbersB; NodeNumbersI; NodeNumbersG]

  for i = 1 : length(FaceNumbers)
     iF = FaceNumbers[i]
     DictF[iF] = i
  end   

  Nodes = map(1:NumNodes + NumNodesG) do i
    Node()
  end
  @inbounds for i = 1 : NumNodes + NumNodesG
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
  end  

  Edges = map(1 : NumEdges + NumEdgesG) do i
    Edge()
  end
  @inbounds for i = 1 : NumEdges + NumEdgesG
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

  Faces = map(1:NumFaces + NumFacesG) do i
    Face()
  end
  @inbounds for i = 1 : NumFaces + NumFacesG
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
  @inbounds for i = 1 : NumEdges + NumEdgesG
    if Edges[i].F[1] in FaceNumbers
      Edges[i].F[1] = GlobalGrid.Faces[Edges[i].F[1]].F
    else
      Edges[i].F[1] = 0
    end
    if Edges[i].F[2] in FaceNumbers
      Edges[i].F[2] = GlobalGrid.Faces[Edges[i].F[2]].F
    else
      Edges[i].F[2] = 0
    end
  end

  Dim=3;
  if GlobalGrid.Type == Quad()
    if order  
      Renumbering!(Edges,Faces);
    end  
  end
# FacesInNodes!(Nodes,Faces)
  FG2F = Dict()
  for iF = 1 : NumFaces + NumFacesG
     iFG = Faces[iF].FG 
     FG2F[iFG] = iF 
  end
  for iN = 1 : NumNodes
    for iF = 1 : length(Nodes[iN].FG)
      iFG = Nodes[iN].FG[iF]  
      Nodes[iN].F[iF] = FG2F[iFG]
    end
  end  

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
  for iE = 1 : length(Edges)
    PosEdgeInFace!(Edges[iE],Edges,Faces)
  end
  @show NumNodes,NumNodesB,NumNodesG
  @views EdgesInNodes!(Nodes[1:NumNodes],Edges,Faces)
 

  H = GlobalGrid.H
  nz = GlobalGrid.nz
  zP=zeros(nz)
  z=KernelAbstractions.zeros(backend,FT,nz+1)
  dzeta=zeros(nz)
  nBar3=zeros(0,0)
  Type = GlobalGrid.Type
  AdaptGrid = GlobalGrid.AdaptGrid
  EF = KernelAbstractions.zeros(backend,Int,0,0)
  FE = KernelAbstractions.zeros(backend,Int,0,0)
  return GridStruct{FT,
                    typeof(EF),
                    typeof(z)}(
    nz,
    zP,
    z,
    dzeta,
    H,
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
    Dim,
    Rad,
    nBar3,
    GlobalGrid.nBar,
    AdaptGrid,
    EF,
    FE,
    )
end
