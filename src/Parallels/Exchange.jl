mutable struct ExchangeStruct{FT<:AbstractFloat,
                              IT1<:AbstractArray,
                              AT3<:AbstractArray}
  IndSendBuffer::Dict{Int,IT1}
  IndSendBufferF::Dict{Int,Array{Int,1}}
  IndRecvBuffer::Dict{Int,IT1}
  IndRecvBufferF::Dict{Int,Array{Int,1}}
  NeiProc::Array{Int, 1}
  Proc::Int
  ProcNumber::Int
  InitSendBuffer::Bool
  InitSendBufferF::Bool
  SendBuffer::Dict
  #SendBuffer3::Dict{Int,Array{FT, 3}}
  SendBuffer3::Dict{Int,AT3}
  SendBufferF::Dict{Int,Array{FT, 4}}
  InitRecvBuffer::Bool
  InitRecvBufferF::Bool
  RecvBuffer::Dict
# RecvBuffer3::Dict{Int,Array{FT, 3}}
  RecvBuffer3::Dict{Int,AT3}
  RecvBufferF::Dict{Int,Array{FT, 4}}
  sreq::MPI.UnsafeMultiRequest
  rreq::MPI.UnsafeMultiRequest
end

function ExchangeStruct{FT}(backend) where FT<:AbstractFloat
  IndSendBuffer = Dict()
  IndSendBufferF = Dict()
  IndRecvBuffer = Dict()
  IndRecvBufferF = Dict()
  NeiProcN = zeros(Int,0)
  Proc = 0
  ProcNumber = 0
  InitSendBuffer = false
  InitSendBufferF = false
  SendBuffer = Dict()
  SendBuffer3 = Dict()
  SendBufferF = Dict()
  InitRecvBuffer = false
  InitRecvBufferF = false
  RecvBuffer = Dict()
  RecvBuffer3 = Dict()
  RecvBufferF = Dict()
  sreq = MPI.UnsafeMultiRequest(0)
  rreq = MPI.UnsafeMultiRequest(0)
  AT3 = KernelAbstractions.zeros(backend,FT,0,0,0)
  IT1 = KernelAbstractions.zeros(backend,Int,0)
  return ExchangeStruct{FT,
                        typeof(IT1),
                        typeof(AT3)}(
    IndSendBuffer,
    IndSendBufferF,
    IndRecvBuffer,
    IndRecvBufferF,
    NeiProcN,
    Proc,
    ProcNumber,
    InitSendBuffer,
    InitSendBufferF,
    SendBuffer,
    SendBuffer3,
    SendBufferF,
    InitRecvBuffer,
    InitRecvBufferF,
    RecvBuffer,
    RecvBuffer3,
    RecvBufferF,
    rreq,
    sreq,
   )
end  


function ExchangeStruct{FT}(backend,SubGrid,OrdPoly,CellToProc,Proc,ProcNumber,HorLimit) where FT<:AbstractFloat

  # Inner Nodes on Edges  
  NumInBoundEdges = 0
  InBoundEdges = zeros(Int,NumInBoundEdges)
  InBoundEdgesP = zeros(Int,NumInBoundEdges)
  NumInBoundEdges = 0
  @inbounds for i = 1:SubGrid.NumEdges
    if CellToProc[SubGrid.Edges[i].FG[1]] == Proc  
    else
      NumInBoundEdges += 1
      push!(InBoundEdges, i)
      push!(InBoundEdgesP, CellToProc[SubGrid.Edges[i].FG[1]])
    end  
    if CellToProc[SubGrid.Edges[i].FG[2]] == Proc  
    else
      NumInBoundEdges += 1
      push!(InBoundEdges, i)
      push!(InBoundEdgesP, CellToProc[SubGrid.Edges[i].FG[1]])
    end  
  end

  NeiProcE = unique(InBoundEdgesP)
  NumNeiProcE = length(NeiProcE)
  DictE=Dict()
  @inbounds for iE = 1 : NumInBoundEdges
    DictE[SubGrid.Edges[InBoundEdges[iE]].EG] = (SubGrid.Edges[InBoundEdges[iE]].E,
      InBoundEdgesP[iE])
  end  

  GlobBuffer = Dict()
  SendBufferE = Dict()
  @inbounds for iP in eachindex(NeiProcE)
    LocTemp = zeros(Int,0)  
    GlobTemp = zeros(Int,0)  
    @inbounds for iEB = 1 : NumInBoundEdges
      if InBoundEdgesP[iEB] == NeiProcE[iP]
        iE = InBoundEdges[iEB] 
        push!(GlobTemp,SubGrid.Edges[iE].EG)
        @inbounds for k = 1 : OrdPoly - 1
          push!(LocTemp,k + SubGrid.NumNodes + (iE - 1)*(OrdPoly - 1))
        end
      end  
    end
    SendBufferE[NeiProcE[iP]] = LocTemp
    GlobBuffer[NeiProcE[iP]] = GlobTemp
  end  

  GlobRecvBuffer = Dict()
  rreq = MPI.UnsafeMultiRequest(length(NeiProcE))
  @inbounds for iP in eachindex(NeiProcE)
    GlobRecvBuffer[NeiProcE[iP]] = similar(GlobBuffer[NeiProcE[iP]])
    tag = Proc + ProcNumber*NeiProcE[iP]
    MPI.Irecv!(GlobRecvBuffer[NeiProcE[iP]], NeiProcE[iP] - 1, tag, MPI.COMM_WORLD, rreq[iP])
  end  
  sreq = MPI.UnsafeMultiRequest(length(NeiProcE))
  @inbounds for iP in eachindex(NeiProcE)
    tag = NeiProcE[iP] + ProcNumber*Proc
    MPI.Isend(GlobBuffer[NeiProcE[iP]], NeiProcE[iP] - 1, tag, MPI.COMM_WORLD, sreq[iP])
  end  

  stats = MPI.Waitall(rreq)
  stats = MPI.Waitall(sreq)

  MPI.Barrier(MPI.COMM_WORLD)
  RecvBufferE = Dict()
  @inbounds for iP in eachindex(NeiProcE)
    GlobInd = GlobRecvBuffer[NeiProcE[iP]]  
    LocTemp=zeros(Int,0)
    iEB1 = 0
    @inbounds for iEB = 1 : NumInBoundEdges
      if InBoundEdgesP[iEB] == NeiProcE[iP]
        iEB1 += 1  
        (iE,) = DictE[GlobInd[iEB1]]  
        @inbounds for k = 1 : OrdPoly - 1
          push!(LocTemp,k + SubGrid.NumNodes + (iE - 1)*(OrdPoly - 1))
        end
      end
    end
    RecvBufferE[NeiProcE[iP]] = LocTemp
  end

  # Nodes 
  NumInBoundNodes = 0
  InBoundNodes = zeros(Int,0)
  InBoundNodesP = zeros(Int,0)
  InBoundNodesFG = zeros(Int,0)
  @inbounds for i = 1 : SubGrid.NumNodes
    NeiProcLoc = zeros(Int,0)  
    NeiFG = zeros(Int,0)  
    @inbounds for iF in eachindex(SubGrid.Nodes[i].FG)
      if CellToProc[SubGrid.Nodes[i].FG[iF]] != Proc
        push!(NeiProcLoc,CellToProc[SubGrid.Nodes[i].FG[iF]])  
        push!(NeiFG,SubGrid.Nodes[i].FG[iF])  
      end
    end
    NeiProcLoc = unique(NeiProcLoc)
    @inbounds for iC in eachindex(NeiProcLoc)
      push!(InBoundNodes,SubGrid.Nodes[i].N)  
      push!(InBoundNodesP,NeiProcLoc[iC])
      push!(InBoundNodesFG,NeiFG[iC])
      NumInBoundNodes += 1
    end  
  end   
  NeiProcN = unique(InBoundNodesP)
  NeiFG = unique(InBoundNodesFG)
  NumNeiProcN = length(NeiProcN)
  NumNeiFG = length(NeiFG)
  DictN=Dict()
  @inbounds for iN = 1 : NumInBoundNodes
    DictN[(SubGrid.Nodes[InBoundNodes[iN]].NG,InBoundNodesP[iN])] = (SubGrid.Nodes[InBoundNodes[iN]].N,
      InBoundNodesP[iN])
  end
  GlobBuffer = Dict()
  SendBufferN = Dict()
  @inbounds for iP in eachindex(NeiProcN)
    LocTemp = zeros(Int,0)
    GlobTemp = zeros(Int,0)
    @inbounds for iNB = 1 : NumInBoundNodes
      if InBoundNodesP[iNB] == NeiProcN[iP]
        iN = InBoundNodes[iNB]
        push!(GlobTemp,SubGrid.Nodes[iN].NG)
        push!(LocTemp,iN)
      end
    end
    SendBufferN[NeiProcN[iP]] = LocTemp
    GlobBuffer[NeiProcN[iP]] = GlobTemp
  end
  GlobRecvBuffer = Dict()
  rreq = MPI.UnsafeMultiRequest(length(NeiProcN))
  @inbounds for iP in eachindex(NeiProcN)
    GlobRecvBuffer[NeiProcN[iP]] = similar(GlobBuffer[NeiProcN[iP]])
    tag = Proc + ProcNumber*NeiProcN[iP]
    MPI.Irecv!(GlobRecvBuffer[NeiProcN[iP]], NeiProcN[iP] - 1, tag, MPI.COMM_WORLD, rreq[iP])
  end  
  sreq = MPI.UnsafeMultiRequest(length(NeiProcN))
  @inbounds for iP in eachindex(NeiProcN)
    tag = NeiProcN[iP] + ProcNumber*Proc
    MPI.Isend(GlobBuffer[NeiProcN[iP]], NeiProcN[iP] - 1, tag, MPI.COMM_WORLD, sreq[iP])
  end  

  stats = MPI.Waitall(rreq)
  stats = MPI.Waitall(sreq)

  MPI.Barrier(MPI.COMM_WORLD)

  RecvBufferN = Dict()
  @inbounds for iP in eachindex(NeiProcN)
    GlobInd = GlobRecvBuffer[NeiProcN[iP]]
    LocTemp=zeros(Int,0)
    iNB1 = 0
    @inbounds for iNB = 1 : NumInBoundNodes
      if InBoundNodesP[iNB] == NeiProcN[iP]
        iNB1 += 1
        (iN,) = DictN[(GlobInd[iNB1],NeiProcN[iP])]
        push!(LocTemp,iN)
      end
    end
    RecvBufferN[NeiProcN[iP]] = LocTemp
  end
  @inbounds for iP in eachindex(NeiProcE)
    RecvBufferN[NeiProcE[iP]] = [RecvBufferN[NeiProcE[iP]];RecvBufferE[NeiProcE[iP]]] 
    SendBufferN[NeiProcE[iP]] = [SendBufferN[NeiProcE[iP]];SendBufferE[NeiProcE[iP]]] 
  end  

  if HorLimit 
#   Ghost faces for horizontal limiter
    NumInBoundFaces = 0
    InBoundFaces = zeros(Int,0)
    InBoundFacesP = zeros(Int,0)
    InBoundFacesFG = zeros(Int,0)
    InBoundFacesF = zeros(Int,0)
    DictF = Dict()
    @inbounds for i in eachindex(SubGrid.Faces)
      InF = zeros(Int,0)
      NeiProcLoc = zeros(Int,0)
      NeiFG = zeros(Int,0)
      NeiF = zeros(Int,0)
      N = SubGrid.Faces[i].N
      for iN in eachindex(SubGrid.Nodes[N])
        FG = SubGrid.Nodes[N[iN]].FG
        F = SubGrid.Nodes[N[iN]].F
        for iF in eachindex(FG)
          if CellToProc[FG[iF]] != Proc
            DictF[(FG[iF],CellToProc[FG[iF]])] = F[iF]  
            push!(InF,i)
            push!(NeiProcLoc,CellToProc[FG[iF]])
            push!(NeiFG,SubGrid.Faces[i].FG)
            push!(NeiF,F[iF])
          end
        end
      end
      @inbounds for iC in eachindex(NeiProcLoc)
        push!(InBoundFaces,InF[iC])
        push!(InBoundFacesP,NeiProcLoc[iC])
        push!(InBoundFacesFG,NeiFG[iC])
        push!(InBoundFacesF,NeiF[iC])
        NumInBoundFaces += 1
      end
    end

    GlobBuffer = Dict()
    IndSendBufferF = Dict()
    IndRecvBufferF = Dict()
    NeiProcF = unique(InBoundFacesP)
    @inbounds for iP in eachindex(NeiProcF)
      LocTemp=zeros(Int,0)
      GlobTemp=zeros(Int,0)
      GlobFTemp=zeros(Int,0)
      @inbounds for iFB = 1 : NumInBoundFaces
        if  InBoundFacesP[iFB] == NeiProcF[iP]
          push!(GlobTemp,InBoundFacesFG[iFB])
          push!(GlobFTemp,InBoundFacesF[iFB])
          push!(LocTemp,InBoundFaces[iFB])
        end
      end
      LocTemp = unique(LocTemp)
      GlobTemp = unique(GlobTemp)
      IndSendBufferF[NeiProcF[iP]] = LocTemp
      IndRecvBufferF[NeiProcF[iP]] = unique(GlobFTemp)
      GlobBuffer[NeiProcF[iP]] = GlobTemp
    end

    GlobRecvBuffer = Dict()
    rreq = MPI.UnsafeMultiRequest(length(NeiProcF))
    @inbounds for iP in eachindex(NeiProcF)
      GlobRecvBuffer[NeiProcF[iP]] = similar(IndRecvBufferF[NeiProcF[iP]])
      tag = Proc + ProcNumber*NeiProcF[iP]
      MPI.Irecv!(GlobRecvBuffer[NeiProcF[iP]], NeiProcF[iP] - 1, tag, MPI.COMM_WORLD, rreq[iP])
    end  
    sreq = MPI.UnsafeMultiRequest(length(NeiProcF))
    @inbounds for iP in eachindex(NeiProcF)
      tag = NeiProcF[iP] + ProcNumber*Proc
      MPI.Isend(GlobBuffer[NeiProcF[iP]], NeiProcF[iP] - 1, tag, MPI.COMM_WORLD, sreq[iP])
    end  

    stats = MPI.Waitall(rreq)
    stats = MPI.Waitall(sreq)
    IndRecvBufferF = Dict()

    @inbounds for iP in eachindex(NeiProcF)
      GlobInd = GlobRecvBuffer[NeiProcF[iP]]
      LocTemp=zeros(Int,0)
      iFB1 = 0
      @inbounds for iFB in eachindex(GlobInd)
        iFB1 += 1
        (iF,) = DictF[(GlobInd[iFB1],NeiProcF[iP])]
        push!(LocTemp,iF)
      end
      IndRecvBufferF[NeiProcF[iP]] = LocTemp
    end
  else
    IndSendBufferF = Dict()
    IndRecvBufferF=Dict()  
  end  

  InitSendBuffer = true
  InitSendBufferF = true
  SendBuffer = Dict()
  SendBuffer3 = Dict{Int,Array{Int,3}}()
  SendBufferF = Dict{Int,Array{Int,4}}()
  @inbounds for iP in eachindex(NeiProcN)
    SendBuffer[iP] = zeros(0)
    SendBuffer3[iP] = zeros(0,0,0)
    SendBufferF[iP] = zeros(0,0,0,0)
  end
  InitRecvBuffer = true
  InitRecvBufferF = true
  RecvBuffer = Dict()
  RecvBuffer3 = Dict{Int,Array{Int,3}}()
  RecvBufferF = Dict{Int,Array{Int,4}}()
  @inbounds for iP in NeiProcN
    RecvBuffer[iP] = zeros(0)
    RecvBuffer3[iP] = zeros(0,0,0)
    RecvBufferF[iP] = zeros(0,0,0,0)
  end

  sreq = MPI.UnsafeMultiRequest(length(NeiProcN))
  rreq = MPI.UnsafeMultiRequest(length(NeiProcN))

  # Copy from CPU to device
  AT3 = KernelAbstractions.zeros(backend,FT,0,0,0)
  IT1 = KernelAbstractions.zeros(backend,Int,0)

  SendBuffer = Dict()
  for (key,) in SendBufferN
    SendBuffer[key] = KernelAbstractions.zeros(backend,Int,size(SendBufferN[key]))
    copyto!(SendBuffer[key],SendBufferN[key])
  end  
  RecvBuffer = Dict()
  for (key,) in RecvBufferN
    RecvBuffer[key] = KernelAbstractions.zeros(backend,Int,size(RecvBufferN[key]))
    copyto!(RecvBuffer[key],RecvBufferN[key])
  end  


  return ExchangeStruct{FT,
                        typeof(IT1),
                        typeof(AT3)}(
    SendBuffer,
    IndSendBufferF,
    RecvBuffer,
    IndRecvBufferF,
    NeiProcN,
    Proc,
    ProcNumber,
    InitSendBuffer,
    InitSendBufferF,
    SendBuffer,
    SendBuffer3,
    SendBufferF,
    InitRecvBuffer,
    InitRecvBufferF,
    RecvBuffer,
    RecvBuffer3,
    RecvBufferF,
    rreq,
    sreq,
   )   
end

function InitExchangeDG(SubGrid,OrdPoly,CellToProc,Proc,ProcNumber,Parallel)

  if Parallel
    NumInBoundEdges = 0
    InBoundEdges = zeros(Int,NumInBoundEdges)
    InBoundEdgesP = zeros(Int,NumInBoundEdges)
    NumInBoundEdges = 0
    @inbounds for i = 1:SubGrid.NumEdges
      if CellToProc[SubGrid.Edges[i].FG[1]] == Proc  
      else
        NumInBoundEdges += 1
        push!(InBoundEdges, i)
        push!(InBoundEdgesP, CellToProc[SubGrid.Edges[i].FG[1]])
      end  
      if CellToProc[SubGrid.Edges[i].FG[2]] == Proc  
      else
        NumInBoundEdges += 1
        push!(InBoundEdges, i)
        push!(InBoundEdgesP, CellToProc[SubGrid.Edges[i].FG[1]])
      end  
    end

    NeiProcE = unique(InBoundEdgesP)
    NumNeiProcE = length(NeiProcE)
    DictE=Dict()
    @inbounds for iE = 1 : NumInBoundEdges
      DictE[SubGrid.Edges[InBoundEdges[iE]].EG] = (SubGrid.Edges[InBoundEdges[iE]].E,
        InBoundEdgesP[iE])
    end  

    GlobBuffer = Dict()
    SendBufferE = Dict()
    @inbounds for iP in eachindex(NeiProcE)
      LocTemp=zeros(Int,0)  
      GlobTemp=zeros(Int,0)  
      @inbounds for iEB = 1 : NumInBoundEdges
        if InBoundEdgesP[iEB] == NeiProcE[iP]
          iE = InBoundEdges[iEB] 
          push!(GlobTemp,SubGrid.Edges[iE].EG)
          @inbounds for k = 1 : OrdPoly - 1
            push!(LocTemp,k + SubGrid.NumNodes + (iE - 1)*(OrdPoly + 1))
          end
        end  
      end
      SendBufferE[NeiProcE[iP]] = LocTemp
      GlobBuffer[NeiProcE[iP]] = GlobTemp
    end  

    GlobRecvBuffer = Dict()
    rreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProcE .- 1)]
    @inbounds for iP in eachindex(NeiProcE)
      GlobRecvBuffer[NeiProcE[iP]] = similar(GlobBuffer[NeiProcE[iP]])
      tag = Proc + ProcNumber*NeiProcE[iP]
      rreq[iP] = MPI.Irecv!(GlobRecvBuffer[NeiProcE[iP]], NeiProcE[iP] - 1, tag, MPI.COMM_WORLD)
    end  
    sreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProcE .- 1)]
    @inbounds for iP in eachindex(NeiProcE)
      tag = NeiProcE[iP] + ProcNumber*Proc
      sreq[iP] = MPI.Isend(GlobBuffer[NeiProcE[iP]], NeiProcE[iP] - 1, tag, MPI.COMM_WORLD)
    end  

    stats = MPI.Waitall!(rreq)
    stats = MPI.Waitall!(sreq)

    MPI.Barrier(MPI.COMM_WORLD)
    RecvBufferE = Dict()
    @inbounds for iP in eachindex(NeiProcE)
      GlobInd = GlobRecvBuffer[NeiProcE[iP]]  
      LocTemp=zeros(Int,0)
      iEB1 = 0
      @inbounds for iEB = 1 : NumInBoundEdges
        if InBoundEdgesP[iEB] == NeiProcE[iP]
          iEB1 += 1  
          (iE,) = DictE[GlobInd[iEB1]]  
          @inbounds for k = 1 : OrdPoly - 1
            push!(LocTemp,k + SubGrid.NumNodes + (iE - 1)*(OrdPoly + 1))
          end
        end
      end
      RecvBufferE[NeiProcE[iP]] = LocTemp
    end

    RecvBufferN = Dict()
    SendBufferN = Dict()
    @inbounds for iP in eachindex(NeiProcE)
      RecvBufferN[NeiProcE[iP]] = [RecvBufferE[NeiProcE[iP]]] 
      SendBufferN[NeiProcE[iP]] = [SendBufferE[NeiProcE[iP]]] 
    end  
  else
    SendBufferN=Dict()  
    RecvBufferN=Dict()  
  end  

  NeiProcN = deepcopy(NeiProcE)
  InitSendBuffer = true
  InitSendBufferF = true
  SendBuffer = Dict()
  SendBuffer3 = Dict{Int,Array{Int,3}}()
  SendBufferF = Dict{Int,Array{Int,4}}()
  @inbounds for iP in eachindex(NeiProcN)
    SendBuffer[iP] = zeros(0)
    SendBuffer3[iP] = zeros(0,0,0)
    SendBufferF[iP] = zeros(0,0,0,0)
  end
  InitRecvBuffer = true
  InitRecvBufferF = true
  RecvBuffer = Dict()
  RecvBuffer3 = Dict{Int,Array{Int,3}}()
  RecvBufferF = Dict{Int,Array{Int,4}}()
  @inbounds for iP in NeiProcN
    RecvBuffer[iP] = zeros(0)
    RecvBuffer3[iP] = zeros(0,0,0)
    RecvBufferF[iP] = zeros(0,0,0,0)
  end
  rreq = MPI.RequestSet(MPI.Request[MPI.REQUEST_NULL for _ in (NeiProcN .- 1)])
  sreq = MPI.RequestSet(MPI.Request[MPI.REQUEST_NULL for _ in (NeiProcN .- 1)])


  return ExchangeStruct(
    SendBufferN,
    RecvBufferN,
    NeiProcN,
    Proc,
    ProcNumber,
    Parallel,
    InitSendBuffer,
    SendBuffer,
    SendBuffer3,
    InitRecvBuffer,
    RecvBuffer,
    RecvBuffer3,
    rreq,
    sreq,
   )   
end


function ExchangeData!(U::Array{Float64,3},Exchange)

  IndSendBuffer = Exchange.IndSendBuffer
  IndRecvBuffer = Exchange.IndRecvBuffer
  NeiProc = Exchange.NeiProc
  Proc = Exchange.Proc
  ProcNumber = Exchange.ProcNumber
  nz = size(U,1)
  nT = size(U,3)
  SendBuffer = Dict()
  @inbounds for iP in eachindex(NeiProc)
    SendBuffer[NeiProc[iP]] = deepcopy(U[:,IndSendBuffer[NeiProc[iP]],:])
  end

  RecvBuffer = Dict()
  rreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProc .- 1)]
  @inbounds for iP in eachindex(NeiProc)
    RecvBuffer[NeiProc[iP]] = zeros(nz,length(IndRecvBuffer[NeiProc[iP]]),nT)
    tag = Proc + ProcNumber*NeiProc[iP]
    rreq[iP] = MPI.Irecv!(RecvBuffer[NeiProc[iP]], NeiProc[iP] - 1, tag, MPI.COMM_WORLD)
  end  
  sreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProc .- 1)]
  @inbounds for iP in eachindex(NeiProc)
    tag = NeiProc[iP] + ProcNumber*Proc
    sreq[iP] = MPI.Isend(SendBuffer[NeiProc[iP]], NeiProc[iP] - 1, tag, MPI.COMM_WORLD)
  end

  stats = MPI.Waitall!(rreq)
  stats = MPI.Waitall!(sreq)

  MPI.Barrier(MPI.COMM_WORLD)
  #Receive
  for iP in eachindex(NeiProc)
    @views @. U[:,IndRecvBuffer[NeiProc[iP]],:] += RecvBuffer[NeiProc[iP]]
  end
end  

function ExchangeData3D!(U,Exchange)

  IndSendBuffer = Exchange.IndSendBuffer
  IndRecvBuffer = Exchange.IndRecvBuffer
  NeiProc = Exchange.NeiProc
  Proc = Exchange.Proc
  ProcNumber = Exchange.ProcNumber
  nz = size(U,1)
  nT = size(U,3)
  if Exchange.InitRecvBuffer
    @inbounds for iP in eachindex(NeiProc)
      Exchange.RecvBuffer[NeiProc[iP]] = zeros(nz,length(IndRecvBuffer[NeiProc[iP]]),nT)
      Exchange.SendBuffer[NeiProc[iP]] = zeros(nz,length(IndRecvBuffer[NeiProc[iP]]),nT)
    end  
    RecvBuffer = Exchange.RecvBuffer
    SendBuffer = Exchange.SendBuffer
    Exchange.InitRecvBuffer = false
    Exchange.InitSendBuffer = false
    Exchange.rreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProc .- 1)]
    Exchange.sreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProc .- 1)]
    rreq = Exchange.rreq
    sreq = Exchange.sreq
  else
    RecvBuffer = Exchange.RecvBuffer
    SendBuffer = Exchange.SendBuffer
    rreq = Exchange.rreq
    sreq = Exchange.sreq
  end    

  @inbounds for iP in NeiProc
    @views @. SendBuffer[iP] = U[:,IndSendBuffer[iP],:]
  end

  i = 0
  @inbounds for iP in NeiProc
    tag = Proc + ProcNumber*iP
    i += 1
    rreq[i] = MPI.Irecv!(RecvBuffer[iP], iP - 1, tag, MPI.COMM_WORLD)
  end  
  i = 0
  @inbounds for iP in NeiProc
    tag = iP + ProcNumber*Proc
    i += 1
    sreq[i] = MPI.Isend(SendBuffer[iP], iP - 1, tag, MPI.COMM_WORLD)
  end

  stats = MPI.Waitall!(rreq)
  stats = MPI.Waitall!(sreq)

  MPI.Barrier(MPI.COMM_WORLD)
  #Receive
  @inbounds for iP in NeiProc
    @views @. U[:,IndRecvBuffer[iP],:] += RecvBuffer[iP]
  end
end    

function ExchangeDataFSend(cFMin,cFMax,Exchange)
  IndSendBufferF = Exchange.IndSendBufferF
  IndRecvBufferF = Exchange.IndRecvBufferF
  NeiProc = Exchange.NeiProc
  Proc = Exchange.Proc
  ProcNumber = Exchange.ProcNumber
  nz = size(cFMin,1)
  nT = size(cFMin,3)
  if Exchange.InitRecvBufferF
    @inbounds for iP in NeiProc
      Exchange.RecvBufferF[iP] = zeros(nz,length(IndRecvBufferF[iP]),nT,2)
      Exchange.SendBufferF[iP] = zeros(nz,length(IndRecvBufferF[iP]),nT,2)
    end
    RecvBufferF = Exchange.RecvBufferF
    SendBufferF = Exchange.SendBufferF
    Exchange.InitRecvBufferF = false
    Exchange.InitSendBufferF = false
    rreq = Exchange.rreq
    sreq = Exchange.sreq
  else
    RecvBufferF = Exchange.RecvBufferF
    SendBufferF = Exchange.SendBufferF
    rreq = Exchange.rreq
    sreq = Exchange.sreq
  end

  @inbounds for iP in NeiProc
    i = 0
    @views @inbounds for Ind in IndSendBufferF[iP]
      i += 1
      @views @. SendBufferF[iP][:,i,:,1] = cFMin[:,Ind,:]
      @views @. SendBufferF[iP][:,i,:,2] = cFMax[:,Ind,:]
    end
  end
  i = 0
  @inbounds for iP in NeiProc
    tag = Proc + ProcNumber*iP
    i += 1
    @views MPI.Irecv!(RecvBufferF[iP], iP - 1, tag, MPI.COMM_WORLD, rreq[i])
  end
  i = 0
  @inbounds for iP in NeiProc
    tag = iP + ProcNumber*Proc
    i += 1
    @views MPI.Isend(SendBufferF[iP], iP - 1, tag, MPI.COMM_WORLD, sreq[i])
  end
end

function ExchangeDataFRecv!(cFMin,cFMax,Exchange)

  IndRecvBufferF = Exchange.IndRecvBufferF
  NeiProc = Exchange.NeiProc
  RecvBufferF = Exchange.RecvBufferF
  rreq = Exchange.rreq
  sreq = Exchange.sreq

  stats = MPI.Waitall(rreq)
  stats = MPI.Waitall(sreq)
  MPI.Barrier(MPI.COMM_WORLD)
  #Receive
  @inbounds for iP in NeiProc
    i = 0
    @inbounds for Ind in IndRecvBufferF[iP]
      i += 1
      @views @. cFMin[:,Ind,:] = RecvBufferF[iP][:,i,:,1]
      @views @. cFMax[:,Ind,:] = RecvBufferF[iP][:,i,:,2]
    end
  end
end

function InitExchangeData3D(backend,FT,nz,nT,Exchange) 
  IndSendBuffer = Exchange.IndSendBuffer
  IndRecvBuffer = Exchange.IndRecvBuffer
  NeiProc = Exchange.NeiProc
  @inbounds for iP in NeiProc
    Exchange.RecvBuffer3[iP] = KernelAbstractions.zeros(backend,FT,nz,length(IndRecvBuffer[iP]),nT)
    Exchange.SendBuffer3[iP] = KernelAbstractions.zeros(backend,FT,nz,length(IndSendBuffer[iP]),nT)
  end  
end  

function ExchangeData3DSend(U,p,Exchange)

  IndSendBuffer = Exchange.IndSendBuffer
  IndRecvBuffer = Exchange.IndRecvBuffer
  NeiProc = Exchange.NeiProc
  Proc = Exchange.Proc
  ProcNumber = Exchange.ProcNumber
  nz = size(U,1)
  nT = size(U,3)
  RecvBuffer3 = Exchange.RecvBuffer3
  SendBuffer3 = Exchange.SendBuffer3
  rreq = Exchange.rreq
  sreq = Exchange.sreq

  @inbounds for iP in NeiProc
    i = 0
    @views @inbounds for Ind in IndSendBuffer[iP]
      i += 1
      @views @. SendBuffer3[iP][:,i,1:nT] = U[:,Ind,:]
      @views @. SendBuffer3[iP][:,i,nT + 1] = p[:,Ind]
    end
  end
  i = 0
  @inbounds for iP in NeiProc
    tag = Proc + ProcNumber*iP
    i += 1
    @views MPI.Irecv!(RecvBuffer3[iP], iP - 1, tag, MPI.COMM_WORLD, rreq[i])
  end  
  i = 0
  @inbounds for iP in NeiProc
    tag = iP + ProcNumber*Proc
    i += 1
    @views MPI.Isend(SendBuffer3[iP], iP - 1, tag, MPI.COMM_WORLD, sreq[i])
  end
end

function ExchangeData3DSend(U,Exchange)

  IndSendBuffer = Exchange.IndSendBuffer
  IndRecvBuffer = Exchange.IndRecvBuffer
  NeiProc = Exchange.NeiProc
  Proc = Exchange.Proc
  ProcNumber = Exchange.ProcNumber
  nz = size(U,1)
  nT = size(U,3)
  RecvBuffer3 = Exchange.RecvBuffer3
  SendBuffer3 = Exchange.SendBuffer3
  rreq = Exchange.rreq
  sreq = Exchange.sreq

  @inbounds for iP in NeiProc
    i = 0
    @views @inbounds for Ind in IndSendBuffer[iP]
      i += 1
      @views @. SendBuffer3[iP][1:nz,i,1:nT] = U[:,Ind,:]
    end
  end
  i = 0
  @inbounds for iP in NeiProc
    tag = Proc + ProcNumber*iP
    i += 1
    @views MPI.Irecv!(RecvBuffer3[iP][1:nz,:,1:nT], iP - 1, tag, MPI.COMM_WORLD, rreq[i])
  end  
  i = 0
  @inbounds for iP in NeiProc
    tag = iP + ProcNumber*Proc
    i += 1
    @views MPI.Isend(SendBuffer3[iP][1:nz,:,1:nT], iP - 1, tag, MPI.COMM_WORLD, sreq[i])
  end
end

function ExchangeData3DSendGPU(U,Exchange)
  backend = get_backend(U)
  FT = eltype(U)

  IndSendBuffer = Exchange.IndSendBuffer
  IndRecvBuffer = Exchange.IndRecvBuffer
  NeiProc = Exchange.NeiProc
  Proc = Exchange.Proc
  ProcNumber = Exchange.ProcNumber
  Nz = size(U,1)
  nT = size(U,3)
  RecvBuffer3 = Exchange.RecvBuffer3
  SendBuffer3 = Exchange.SendBuffer3
  rreq = Exchange.rreq
  sreq = Exchange.sreq

  group = (Nz,5,1)
  KExchangeData3DSendKernel! = ExchangeData3DSendKernel!(backend,group)
  @inbounds for iP in NeiProc
    ndrange = (Nz,length(IndSendBuffer[iP]),nT)
    KExchangeData3DSendKernel!(U,SendBuffer3[iP],IndSendBuffer[iP],ndrange=ndrange)
  end

  i = 0
  @inbounds for iP in NeiProc
    tag = Proc + ProcNumber*iP
    i += 1
    @views MPI.Irecv!(RecvBuffer3[iP][1:Nz,:,1:nT], iP - 1, tag, MPI.COMM_WORLD, rreq[i])
  end  
  i = 0
  @inbounds for iP in NeiProc
    tag = iP + ProcNumber*Proc
    i += 1
    @views MPI.Isend(SendBuffer3[iP][1:Nz,:,1:nT], iP - 1, tag, MPI.COMM_WORLD, sreq[i])
  end
end

@kernel function ExchangeData3DSendKernel!(U,SendBuffer,IndSendBuffer)

  Iz,I,IT = @index(Global, NTuple)
  NumInd = @uniform @ndrange()[2]
  NT = @uniform @ndrange()[3]

  if I <= NumInd && IT <= NT
    @inbounds Ind = IndSendBuffer[I]    
    @inbounds SendBuffer[Iz,I,IT] = U[Iz,Ind,IT]
  end  
end

function ExchangeData3DRecv!(U,p,Exchange)

  nT = size(U,3)
  IndRecvBuffer = Exchange.IndRecvBuffer
  NeiProc = Exchange.NeiProc
  RecvBuffer3 = Exchange.RecvBuffer3
  rreq = Exchange.rreq
  sreq = Exchange.sreq

  stats = MPI.Waitall(rreq)
  stats = MPI.Waitall(sreq)
  MPI.Barrier(MPI.COMM_WORLD)
  #Receive
  @inbounds for iP in NeiProc
    i = 0
    @inbounds for Ind in IndRecvBuffer[iP]
      i += 1
      @views @. U[:,Ind,:] += RecvBuffer3[iP][:,i,1:nT]
      @views @. p[:,Ind] += RecvBuffer3[iP][:,i,nT+1]
    end
  end
end  

function ExchangeData3DRecv!(U,Exchange)

  nz = size(U,1)
  nT = size(U,3)
  IndRecvBuffer = Exchange.IndRecvBuffer
  NeiProc = Exchange.NeiProc
  RecvBuffer3 = Exchange.RecvBuffer3
  rreq = Exchange.rreq
  sreq = Exchange.sreq

  stats = MPI.Waitall(rreq)
  stats = MPI.Waitall(sreq)
  MPI.Barrier(MPI.COMM_WORLD)
  #Receive
  @inbounds for iP in NeiProc
    i = 0
    @inbounds for Ind in IndRecvBuffer[iP]
      i += 1
      @views @. U[1:nz,Ind,:] += RecvBuffer3[iP][1:nz,i,1:nT]
    end
  end
end  

function ExchangeData3DRecvGPU!(U,Exchange)
  backend = get_backend(U)
  FT = eltype(U)

  Nz = size(U,1)
  nT = size(U,3)
  IndRecvBuffer = Exchange.IndRecvBuffer
  NeiProc = Exchange.NeiProc
  Proc = Exchange.Proc
  RecvBuffer3 = Exchange.RecvBuffer3
  rreq = Exchange.rreq
  sreq = Exchange.sreq

  stats = MPI.Waitall(rreq)
  stats = MPI.Waitall(sreq)
  MPI.Barrier(MPI.COMM_WORLD)

  group = (Nz,5,1)
  KExchangeData3DRecvKernel! = ExchangeData3DRecvKernel!(backend,group)

  #Receive
  @inbounds for iP in NeiProc
    ndrange = (Nz,length(IndRecvBuffer[iP]),nT)
    KExchangeData3DRecvKernel!(U,RecvBuffer3[iP],IndRecvBuffer[iP],ndrange=ndrange)
  end
end  

@kernel function ExchangeData3DRecvKernel!(U,RecvBuffer,IndRecvBuffer)

  Iz,I,IT = @index(Global, NTuple)
  NumInd = @uniform @ndrange()[2]
  NT = @uniform @ndrange()[3]
  if I <= NumInd && IT <= NT
    @inbounds Ind = IndRecvBuffer[I]    
    @inbounds @atomic U[Iz,Ind,IT] += RecvBuffer[Iz,I,IT]
  end  
end  

function ExchangeData!(U::AbstractArray{FT,2},Exchange) where FT<:AbstractFloat

  nz = size(U,1)

  IndSendBuffer = Exchange.IndSendBuffer
  IndRecvBuffer = Exchange.IndRecvBuffer
  NeiProc = Exchange.NeiProc
  Proc = Exchange.Proc
  ProcNumber = Exchange.ProcNumber
  SendBuffer = Dict()
  @inbounds for iP in eachindex(NeiProc)
    SendBuffer[NeiProc[iP]] = U[:,IndSendBuffer[NeiProc[iP]]]
  end

  RecvBuffer = Dict()
  rreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProc .- 1)]
  @inbounds for iP in eachindex(NeiProc)
    RecvBuffer[NeiProc[iP]] = zeros(nz,length(IndRecvBuffer[NeiProc[iP]]))
    tag = Proc + ProcNumber*NeiProc[iP]
    rreq[iP] = MPI.Irecv!(RecvBuffer[NeiProc[iP]], NeiProc[iP] - 1, tag, MPI.COMM_WORLD)
  end  
  sreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProc .- 1)]
  @inbounds for iP in eachindex(NeiProc)
    tag = NeiProc[iP] + ProcNumber*Proc
    sreq[iP] = MPI.Isend(SendBuffer[NeiProc[iP]], NeiProc[iP] - 1, tag, MPI.COMM_WORLD)
  end

  stats = MPI.Waitall!(rreq)
  stats = MPI.Waitall!(sreq)

  MPI.Barrier(MPI.COMM_WORLD)
  #Receive
  @inbounds for iP in eachindex(NeiProc)
    U[:,IndRecvBuffer[NeiProc[iP]]] .+= RecvBuffer[NeiProc[iP]]
  end
end

function ExchangeData!(U::AbstractArray{FT,1},Exchange) where FT<:AbstractFloat

  IndSendBuffer = Exchange.IndSendBuffer
  IndRecvBuffer = Exchange.IndRecvBuffer
  NeiProc = Exchange.NeiProc
  Proc = Exchange.Proc
  ProcNumber = Exchange.ProcNumber
  SendBuffer = Dict()
  @inbounds for iP in eachindex(NeiProc)
    SendBuffer[NeiProc[iP]] = U[IndSendBuffer[NeiProc[iP]]]
  end

  RecvBuffer = Dict()
  rreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProc .- 1)]
  @inbounds for iP in eachindex(NeiProc)
    RecvBuffer[NeiProc[iP]] = zeros(length(IndRecvBuffer[NeiProc[iP]]))
    tag = Proc + ProcNumber*NeiProc[iP]
    rreq[iP] = MPI.Irecv!(RecvBuffer[NeiProc[iP]], NeiProc[iP] - 1, tag, MPI.COMM_WORLD)
  end  
  sreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProc .- 1)]
  @inbounds for iP in eachindex(NeiProc)
    tag = NeiProc[iP] + ProcNumber*Proc
    sreq[iP] = MPI.Isend(SendBuffer[NeiProc[iP]], NeiProc[iP] - 1, tag, MPI.COMM_WORLD)
  end

  stats = MPI.Waitall!(rreq)
  stats = MPI.Waitall!(sreq)

  MPI.Barrier(MPI.COMM_WORLD)
  #Receive
  @inbounds for iP in eachindex(NeiProc)
    U[IndRecvBuffer[NeiProc[iP]]] .+= RecvBuffer[NeiProc[iP]]
  end
end

function GlobalSum2D(U,CG)
  SumLoc = 0.0
  @inbounds for iG = 1 : CG.NumG
    SumLoc += sum(abs.(U[:,iG] .* CG.MasterSlave[iG]))
  end  
  Sum = MPI.Allreduce(SumLoc, +, MPI.COMM_WORLD)
end  

function GlobalSum3D(U,CG)
  SumLoc = 0.0
  @inbounds for iG = 1 : CG.NumG
    SumLoc += sum(abs.(U[:,iG,:] .* CG.MasterSlave[iG]))
  end  
  Sum = MPI.Allreduce(SumLoc, +, MPI.COMM_WORLD)
end  

function GlobalIntegral(c,CG,Global)
  SumLoc = sum(c.*CG.MMass)
  SumGlobal = MPI.Allreduce(SumLoc, +, MPI.COMM_WORLD)
end
