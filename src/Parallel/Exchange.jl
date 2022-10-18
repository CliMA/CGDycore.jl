mutable struct ExchangeStruct
  IndSendBuffer::Dict{Int,Array{Int,1}}
  IndRecvBuffer::Dict{Int,Array{Int,1}}
  NeiProc::Array{Int, 1}
  Proc::Int
  ProcNumber::Int
  Parallel::Bool
  InitSendBuffer::Bool
  SendBuffer::Dict
  SendBuffer3::Dict{Int,Array{Float64, 3}}
  InitRecvBuffer::Bool
  RecvBuffer::Dict
  RecvBuffer3::Dict{Int,Array{Float64, 3}}
# rreq::Array{MPI.Request, 1}
# sreq::Array{MPI.Request, 1}
# rreq::Vector{MPI_Request}
  sreq::MPI.RequestSet
  rreq::MPI.RequestSet
end

function InitExchange(SubGrid,OrdPoly,CellToProc,Proc,ProcNumber,Parallel)

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
            push!(LocTemp,k + SubGrid.NumNodes + (iE - 1)*(OrdPoly - 1))
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
            push!(LocTemp,k + SubGrid.NumNodes + (iE - 1)*(OrdPoly - 1))
          end
        end
      end
      RecvBufferE[NeiProcE[iP]] = LocTemp
    end

    NumInBoundNodes = 0
    InBoundNodes = zeros(Int,0)
    InBoundNodesP = zeros(Int,0)
    @inbounds for i = 1 : SubGrid.NumNodes
      NeiProcLoc = zeros(Int,0)  
      @inbounds for iF in eachindex(SubGrid.Nodes[i].FG)
        if CellToProc[SubGrid.Nodes[i].FG[iF]] != Proc
          push!(NeiProcLoc,CellToProc[SubGrid.Nodes[i].FG[iF]])  
        end
      end
      NeiProcLoc = unique(NeiProcLoc)
      @inbounds for iC in eachindex(NeiProcLoc)
        push!(InBoundNodes,SubGrid.Nodes[i].N)  
        push!(InBoundNodesP,NeiProcLoc[iC])
        NumInBoundNodes += 1
      end  
    end   
    NeiProcN = unique(InBoundNodesP)
    NumNeiProcN = length(NeiProcN)
    DictN=Dict()
    @inbounds for iN = 1 : NumInBoundNodes
      DictN[(SubGrid.Nodes[InBoundNodes[iN]].NG,InBoundNodesP[iN])] = (SubGrid.Nodes[InBoundNodes[iN]].N,
        InBoundNodesP[iN])
    end
    GlobBuffer = Dict()
    SendBufferN = Dict()
    @inbounds for iP in eachindex(NeiProcN)
      LocTemp=zeros(Int,0)
      GlobTemp=zeros(Int,0)
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
    rreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProcN .- 1)]
    @inbounds for iP in eachindex(NeiProcN)
      GlobRecvBuffer[NeiProcN[iP]] = similar(GlobBuffer[NeiProcN[iP]])
      tag = Proc + ProcNumber*NeiProcN[iP]
      rreq[iP] = MPI.Irecv!(GlobRecvBuffer[NeiProcN[iP]], NeiProcN[iP] - 1, tag, MPI.COMM_WORLD)
    end  
    sreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProcN .- 1)]
    @inbounds for iP in eachindex(NeiProcN)
      tag = NeiProcN[iP] + ProcNumber*Proc
      sreq[iP] = MPI.Isend(GlobBuffer[NeiProcN[iP]], NeiProcN[iP] - 1, tag, MPI.COMM_WORLD)
    end  

    stats = MPI.Waitall!(rreq)
    stats = MPI.Waitall!(sreq)

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
  else
    SendBufferN=Dict()  
    RecvBufferN=Dict()  
    NeiProcN=zeros(Int,0)
  end  

  InitSendBuffer = true
  SendBuffer = Dict()
  SendBuffer3 = Dict{Int,Array{Int,3}}()
  @inbounds for iP in eachindex(NeiProcN)
    SendBuffer[iP] = zeros(0)
    SendBuffer3[iP] = zeros(0,0,0)
  end
  InitRecvBuffer = true
  RecvBuffer = Dict()
  RecvBuffer3 = Dict{Int,Array{Int,3}}()
  @inbounds for iP in NeiProcN
    RecvBuffer[iP] = zeros(0)
    RecvBuffer3[iP] = zeros(0,0,0)
  end
  #rreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProcN .- 1)]
  #sreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProcN .- 1)]
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

  if Exchange.Parallel

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
end    

function ExchangeData3D!(U,Exchange)

  if Exchange.Parallel

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
      @views SendBuffer[iP] .= U[:,IndSendBuffer[iP],:]
    end

#   rreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProc .- 1)]
    i = 0
    @inbounds for iP in NeiProc
      tag = Proc + ProcNumber*iP
      i += 1
      rreq[i] = MPI.Irecv!(RecvBuffer[iP], iP - 1, tag, MPI.COMM_WORLD)
    end  
#   sreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProc .- 1)]
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
      i = 0
      @inbounds for Ind in IndRecvBuffer[iP]
        i += 1
        @views @. U[:,Ind,:] += RecvBuffer[iP][:,i,:]
      end
    end
  end  
end    

function ExchangeData3DSend(U,p,Exchange)

  if Exchange.Parallel

    IndSendBuffer = Exchange.IndSendBuffer
    IndRecvBuffer = Exchange.IndRecvBuffer
    NeiProc = Exchange.NeiProc
    Proc = Exchange.Proc
    ProcNumber = Exchange.ProcNumber
    nz = size(U,1)
    nT = size(U,3)
    if Exchange.InitRecvBuffer
      @inbounds for iP in NeiProc
        Exchange.RecvBuffer3[iP] = zeros(nz,length(IndRecvBuffer[iP]),nT + 1)
        Exchange.SendBuffer3[iP] = zeros(nz,length(IndRecvBuffer[iP]),nT + 1)
      end  
      RecvBuffer3 = Exchange.RecvBuffer3
      SendBuffer3 = Exchange.SendBuffer3
      Exchange.InitRecvBuffer = false
      Exchange.InitSendBuffer = false
#     Exchange.rreq = MPI.RequestSet[MPI.REQUEST_NULL for _ in (NeiProc .- 1)]
#     Exchange.sreq = MPI.RequestSet[MPI.REQUEST_NULL for _ in (NeiProc .- 1)]
      Exchange.rreq = MPI.RequestSet(MPI.Request[MPI.REQUEST_NULL for _ in (NeiProc .- 1)])
      Exchange.sreq = MPI.RequestSet(MPI.Request[MPI.REQUEST_NULL for _ in (NeiProc .- 1)])
      rreq = Exchange.rreq
      sreq = Exchange.sreq
    else
      RecvBuffer3 = Exchange.RecvBuffer3
      SendBuffer3 = Exchange.SendBuffer3
      rreq = Exchange.rreq
      sreq = Exchange.sreq
    end    

    @inbounds for iP in NeiProc
      i = 0
      @inbounds for Ind in IndSendBuffer[iP]
        i += 1
        @views @. SendBuffer3[iP][:,i,1:nT] = U[:,Ind,:]
        @views @. SendBuffer3[iP][:,i,nT + 1] = p[:,Ind]
      end
    end
#  Test  
#   @time rreq = [MPI.Irecv!(RecvBuffer3[iP], iP - 1, Proc + ProcNumber*iP, MPI.COMM_WORLD) for iP in NeiProc]
#  Test  
    i = 0
    @inbounds for iP in NeiProc
      tag = Proc + ProcNumber*iP
      i += 1
      rreq[i] = MPI.Irecv!(RecvBuffer3[iP], iP - 1, tag, MPI.COMM_WORLD)
    end  
    i = 0
    @inbounds for iP in NeiProc
      tag = iP + ProcNumber*Proc
      i += 1
      sreq[i] = MPI.Isend(SendBuffer3[iP], iP - 1, tag, MPI.COMM_WORLD)
#     count = length(SendBuffer3[iP][:,:,:])
#     datatype = Float64
#     dest = iP - 1
#     @views API.MPI_Isend(SendBuffer3[iP][:,:,:], count, datatype, dest, tag, MPI.COMM_WORLD, sreq[i])
    end
  end
end  

function ExchangeData3DRecv!(U,p,Exchange)

  if Exchange.Parallel

    nT = size(U,3)
    IndRecvBuffer = Exchange.IndRecvBuffer
    NeiProc = Exchange.NeiProc
    RecvBuffer3 = Exchange.RecvBuffer3
    rreq = Exchange.rreq
    sreq = Exchange.sreq

#   stats = MPI.Waitall!(rreq)
#   stats = MPI.Waitall!(sreq)
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
end    

function ExchangeData!(U::Array{Float64,2},Exchange)

  if Exchange.Parallel
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
end  
function ExchangeData!(U::Array{Float64,1},Exchange)

  if Exchange.Parallel

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
