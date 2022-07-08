mutable struct ExchangeStruct
  SendBuffer::Dict
  RecvBuffer::Dict
  NeiProc::Array{Int, 1}
  Proc::Int
  ProcNumber::Int
  Parallel::Bool
# rreq::Vector{MPI_Request}
# sreq::MPI_Request
end

function InitExchange(SubGrid,OrdPoly,CellToProc,Proc,ProcNumber,Parallel)

  if Parallel
    NumInBoundEdges = 0
    InBoundEdges = zeros(Int,NumInBoundEdges)
    InBoundEdgesP = zeros(Int,NumInBoundEdges)
    NumInBoundEdges = 0
    for i = 1:SubGrid.NumEdges
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
    for iE = 1 : NumInBoundEdges
      DictE[SubGrid.Edges[InBoundEdges[iE]].EG] = (SubGrid.Edges[InBoundEdges[iE]].E,
        InBoundEdgesP[iE])
    end  

    GlobBuffer = Dict()
    SendBuffer = Dict()
    for iP in eachindex(NeiProcE)
      LocTemp=zeros(Int,0)  
      GlobTemp=zeros(Int,0)  
      for iEB = 1 : NumInBoundEdges
        if InBoundEdgesP[iEB] == NeiProcE[iP]
          iE = InBoundEdges[iEB] 
          push!(GlobTemp,SubGrid.Edges[iE].EG)
          for k = 1 : OrdPoly - 1
            push!(LocTemp,k + SubGrid.NumNodes + (iE - 1)*(OrdPoly - 1))
          end
        end  
      end
      SendBuffer[NeiProcE[iP]] = LocTemp
      GlobBuffer[NeiProcE[iP]] = GlobTemp
    end  

    GlobRecvBuffer = Dict()
    rreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProcE .- 1)]
    for iP in eachindex(NeiProcE)
      GlobRecvBuffer[NeiProcE[iP]] = similar(GlobBuffer[NeiProcE[iP]])
      tag = Proc + ProcNumber*NeiProcE[iP]
      rreq[iP] = MPI.Irecv!(GlobRecvBuffer[NeiProcE[iP]], NeiProcE[iP] - 1, tag, MPI.COMM_WORLD)
    end  
    sreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProcE .- 1)]
    for iP in eachindex(NeiProcE)
      tag = NeiProcE[iP] + ProcNumber*Proc
      sreq[iP] = MPI.Isend(GlobBuffer[NeiProcE[iP]], NeiProcE[iP] - 1, tag, MPI.COMM_WORLD)
    end  

    stats = MPI.Waitall!(rreq)
    stats = MPI.Waitall!(sreq)

    MPI.Barrier(MPI.COMM_WORLD)
    RecvBuffer = Dict()
    for iP in eachindex(NeiProcE)
      GlobInd = GlobRecvBuffer[NeiProcE[iP]]  
      LocTemp=zeros(Int,0)
      iEB1 = 0
      for iEB = 1 : NumInBoundEdges
        if InBoundEdgesP[iEB] == NeiProcE[iP]
          iEB1 += 1  
          (iE,) = DictE[GlobInd[iEB1]]  
          for k = 1 : OrdPoly - 1
            push!(LocTemp,k + SubGrid.NumNodes + (iE - 1)*(OrdPoly - 1))
          end
        end
      end
      RecvBuffer[NeiProcE[iP]] = LocTemp
    end

    NumInBoundNodes = 0
    InBoundNodes = zeros(Int,0)
    InBoundNodesP = zeros(Int,0)
    for i = 1 : SubGrid.NumNodes
      NeiProcLoc = zeros(Int,0)  
      for iF in eachindex(SubGrid.Nodes[i].FG)
        if CellToProc[SubGrid.Nodes[i].FG[iF]] != Proc
          push!(NeiProcLoc,CellToProc[SubGrid.Nodes[i].FG[iF]])  
        end
      end
      NeiProcLoc = unique(NeiProcLoc)
      for iC in eachindex(NeiProcLoc)
        push!(InBoundNodes,SubGrid.Nodes[i].N)  
        push!(InBoundNodesP,NeiProcLoc[iC])
        NumInBoundNodes += 1
      end  
    end   
    NeiProcN = unique(InBoundNodesP)
    NumNeiProcN = length(NeiProcN)
    DictN=Dict()
    for iN = 1 : NumInBoundNodes
      DictN[(SubGrid.Nodes[InBoundNodes[iN]].NG,InBoundNodesP[iN])] = (SubGrid.Nodes[InBoundNodes[iN]].N,
        InBoundNodesP[iN])
    end
    GlobBuffer = Dict()
    SendBufferN = Dict()
    for iP in eachindex(NeiProcN)
      LocTemp=zeros(Int,0)
      GlobTemp=zeros(Int,0)
      for iNB = 1 : NumInBoundNodes
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
    rreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProcN .- 1)]
    for iP in eachindex(NeiProcN)
      GlobRecvBuffer[NeiProcN[iP]] = similar(GlobBuffer[NeiProcN[iP]])
      tag = Proc + ProcNumber*NeiProcN[iP]
      rreq[iP] = MPI.Irecv!(GlobRecvBuffer[NeiProcN[iP]], NeiProcN[iP] - 1, tag, MPI.COMM_WORLD)
    end  
    sreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProcN .- 1)]
    for iP in eachindex(NeiProcN)
      tag = NeiProcN[iP] + ProcNumber*Proc
      sreq[iP] = MPI.Isend(GlobBuffer[NeiProcN[iP]], NeiProcN[iP] - 1, tag, MPI.COMM_WORLD)
    end  

    stats = MPI.Waitall!(rreq)
    stats = MPI.Waitall!(sreq)

    MPI.Barrier(MPI.COMM_WORLD)

    RecvBufferN = Dict()
    for iP in eachindex(NeiProcN)
      GlobInd = GlobRecvBuffer[NeiProcN[iP]]
      LocTemp=zeros(Int,0)
      iNB1 = 0
      for iNB = 1 : NumInBoundNodes
        if InBoundNodesP[iNB] == NeiProcN[iP]
          iNB1 += 1
          (iN,) = DictN[(GlobInd[iNB1],NeiProcN[iP])]
          push!(LocTemp,iN)
        end
      end
      RecvBufferN[NeiProcN[iP]] = LocTemp
    end
    for iP in eachindex(NeiProcE)
      RecvBufferN[NeiProcE[iP]] = [RecvBufferN[NeiProcE[iP]];RecvBuffer[NeiProcE[iP]]] 
      SendBufferN[NeiProcE[iP]] = [SendBufferN[NeiProcE[iP]];SendBuffer[NeiProcE[iP]]] 
    end  
  else
    SendBufferN=Dict()  
    RecvBufferN=Dict()  
    NeiProcN=zeros(Int,0)
  end  

  return ExchangeStruct(
    SendBufferN,
    RecvBufferN,
    NeiProcN,
    Proc,
    ProcNumber,
    Parallel,
#   rreq,
#   sreq,
   )   
end


function ExchangeData!(U::Array{Float64,3},Exchange)

  if Exchange.Parallel

    IndSendBuffer = Exchange.SendBuffer
    IndRecvBuffer = Exchange.RecvBuffer
    NeiProc = Exchange.NeiProc
    Proc = Exchange.Proc
    ProcNumber = Exchange.ProcNumber
    nz = size(U,1)
    nT = size(U,3)
    SendBuffer = Dict()
    for iP in eachindex(NeiProc)
      SendBuffer[NeiProc[iP]] = deepcopy(U[:,IndSendBuffer[NeiProc[iP]],:])
    end

    RecvBuffer = Dict()
    rreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProc .- 1)]
    for iP in eachindex(NeiProc)
      RecvBuffer[NeiProc[iP]] = zeros(nz,length(IndRecvBuffer[NeiProc[iP]]),nT)
      tag = Proc + ProcNumber*NeiProc[iP]
      rreq[iP] = MPI.Irecv!(RecvBuffer[NeiProc[iP]], NeiProc[iP] - 1, tag, MPI.COMM_WORLD)
    end  
    sreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProc .- 1)]
    for iP in eachindex(NeiProc)
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
end    

function ExchangeData3D!(U,Exchange)

  if Exchange.Parallel

    IndSendBuffer = Exchange.SendBuffer
    IndRecvBuffer = Exchange.RecvBuffer
    NeiProc = Exchange.NeiProc
    Proc = Exchange.Proc
    ProcNumber = Exchange.ProcNumber
    nz = size(U,1)
    nT = size(U,3)
    SendBuffer = Dict()
    for iP in eachindex(NeiProc)
      SendBuffer[NeiProc[iP]] = deepcopy(U[:,IndSendBuffer[NeiProc[iP]],:])
    end

    RecvBuffer = Dict()
    rreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProc .- 1)]
    for iP in eachindex(NeiProc)
      RecvBuffer[NeiProc[iP]] = zeros(nz,length(IndRecvBuffer[NeiProc[iP]]),nT)
      tag = Proc + ProcNumber*NeiProc[iP]
      rreq[iP] = MPI.Irecv!(RecvBuffer[NeiProc[iP]], NeiProc[iP] - 1, tag, MPI.COMM_WORLD)
    end  
    sreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProc .- 1)]
    for iP in eachindex(NeiProc)
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
end    

function ExchangeData!(U::Array{Float64,2},Exchange)

  if Exchange.Parallel
    nz = size(U,1)

    IndSendBuffer = Exchange.SendBuffer
    IndRecvBuffer = Exchange.RecvBuffer
    NeiProc = Exchange.NeiProc
    Proc = Exchange.Proc
    ProcNumber = Exchange.ProcNumber
    SendBuffer = Dict()
    for iP in eachindex(NeiProc)
      SendBuffer[NeiProc[iP]] = deepcopy(U[:,IndSendBuffer[NeiProc[iP]]])
    end

    RecvBuffer = Dict()
    rreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProc .- 1)]
    for iP in eachindex(NeiProc)
      RecvBuffer[NeiProc[iP]] = zeros(nz,length(IndRecvBuffer[NeiProc[iP]]))
      tag = Proc + ProcNumber*NeiProc[iP]
      rreq[iP] = MPI.Irecv!(RecvBuffer[NeiProc[iP]], NeiProc[iP] - 1, tag, MPI.COMM_WORLD)
    end  
    sreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProc .- 1)]
    for iP in eachindex(NeiProc)
      tag = NeiProc[iP] + ProcNumber*Proc
      sreq[iP] = MPI.Isend(SendBuffer[NeiProc[iP]], NeiProc[iP] - 1, tag, MPI.COMM_WORLD)
    end

    stats = MPI.Waitall!(rreq)
    stats = MPI.Waitall!(sreq)

    MPI.Barrier(MPI.COMM_WORLD)
    #Receive
    for iP in eachindex(NeiProc)
      U[:,IndRecvBuffer[NeiProc[iP]]] .+= RecvBuffer[NeiProc[iP]]
    end
  end  
end  
function ExchangeData!(U::Array{Float64,1},Exchange)

  if Exchange.Parallel

    IndSendBuffer = Exchange.SendBuffer
    IndRecvBuffer = Exchange.RecvBuffer
    NeiProc = Exchange.NeiProc
    Proc = Exchange.Proc
    ProcNumber = Exchange.ProcNumber
    SendBuffer = Dict()
    for iP in eachindex(NeiProc)
      SendBuffer[NeiProc[iP]] = deepcopy(U[IndSendBuffer[NeiProc[iP]]])
    end

    RecvBuffer = Dict()
    rreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProc .- 1)]
    for iP in eachindex(NeiProc)
      RecvBuffer[NeiProc[iP]] = zeros(length(IndRecvBuffer[NeiProc[iP]]))
      tag = Proc + ProcNumber*NeiProc[iP]
      rreq[iP] = MPI.Irecv!(RecvBuffer[NeiProc[iP]], NeiProc[iP] - 1, tag, MPI.COMM_WORLD)
    end  
    sreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProc .- 1)]
    for iP in eachindex(NeiProc)
      tag = NeiProc[iP] + ProcNumber*Proc
      sreq[iP] = MPI.Isend(SendBuffer[NeiProc[iP]], NeiProc[iP] - 1, tag, MPI.COMM_WORLD)
    end

    stats = MPI.Waitall!(rreq)
    stats = MPI.Waitall!(sreq)

    MPI.Barrier(MPI.COMM_WORLD)
    #Receive
    for iP in eachindex(NeiProc)
      U[IndRecvBuffer[NeiProc[iP]]] .+= RecvBuffer[NeiProc[iP]]
    end
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
