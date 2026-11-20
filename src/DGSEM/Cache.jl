mutable struct CacheAuxStructDG{FT<:AbstractFloat,
                              AT3<:AbstractArray,
                              AT4<:AbstractArray,
                              IT1<:AbstractArray,
                              IT3<:AbstractArray}
  Aux::AT4
  KV::AT3
  # horizontal interface fluxes without atomics: edge nodes whose element node
  # gets contributions from several edges (corners) write into the slot
  # SlotE[side,I,IE] of FC, the slots of the shared node IndS[i] are
  # PtrS[i]:PtrS[i+1]-1 and are summed in this order
  FC::AT4
  SlotE::IT3
  IndS::IT1
  PtrS::IT1
end
  
function CacheAuxStruct(backend,FT,FE::FiniteElements.DGElement,M,nz,Model,Grid)
  NumG = FE.NumG
  NumI = FE.NumI
  Aux = KernelAbstractions.zeros(backend,FT,M,nz,NumG,Model.NumAux)
  KV = KernelAbstractions.zeros(backend,FT,M,nz,NumI)
  SlotECPU, IndSCPU, PtrSCPU = SharedEdgeNodes(Array(FE.GlobE),Array(Grid.EF),Grid.NumFaces,NumG)
  FC = KernelAbstractions.zeros(backend,FT,M,nz,max(length(IndSCPU) > 0 ? PtrSCPU[end] - 1 : 0,1),Model.NumV)
  SlotE = KernelAbstractions.zeros(backend,Int,size(SlotECPU))
  copyto!(SlotE,SlotECPU)
  IndS = KernelAbstractions.zeros(backend,Int,length(IndSCPU))
  copyto!(IndS,IndSCPU)
  PtrS = KernelAbstractions.zeros(backend,Int,length(PtrSCPU))
  copyto!(PtrS,PtrSCPU)
  return CacheAuxStructDG{FT,
                        typeof(KV),
                        typeof(Aux),
                        typeof(IndS),
                        typeof(SlotE)}(
    Aux,
    KV,
    FC,
    SlotE,
    IndS,
    PtrS,
  )
end

# Element nodes that are written by more than one edge node (only edges whose
# side lies in a face 1:NF are written), slots in the order edge, edge node, side
function SharedEdgeNodes(GlobE,EF,NF,NumG)
  DoFE = size(GlobE,2)
  NE = size(GlobE,3)
  Count = zeros(Int,NumG)
  for IE = 1 : NE, I = 1 : DoFE, s = 1 : 2
    if 0 < EF[s,IE] <= NF
      Count[GlobE[s,I,IE]] += 1
    end
  end
  IndS = findall(c -> c > 1, Count)
  PtrS = ones(Int,length(IndS)+1)
  Pos = zeros(Int,NumG)
  for i = 1 : length(IndS)
    PtrS[i+1] = PtrS[i] + Count[IndS[i]]
    Pos[IndS[i]] = PtrS[i]
  end
  SlotE = zeros(Int,2,DoFE,NE)
  for IE = 1 : NE, I = 1 : DoFE, s = 1 : 2
    ind = GlobE[s,I,IE]
    if 0 < EF[s,IE] <= NF && Count[ind] > 1
      SlotE[s,I,IE] = Pos[ind]
      Pos[ind] += 1
    end
  end
  return SlotE, IndS, PtrS
end

# Vertically implicit solver of the HEVI scheme
#   HDG           : hybridized DG, scalar tridiagonal system for the face traces
#   TwoSidedSplit : two-sided interface condensation with split-form gravity,
#                   block tridiagonal system with 4x4 blocks
#   OneSided      : one-sided interface condensation (arXiv:2610.15775, Sec. 4.2.2),
#                   block tridiagonal system with 3x3 blocks
function JacHEVI(backend,FT,M,nz,FE::FiniteElements.DGElement,HEVISolver)
  if HEVISolver == "HDG"
    return JacHDG(backend,FT,M,nz,FE)
  elseif HEVISolver == "TwoSidedSplit"
    return JacTwoSidedSplit(backend,FT,M,nz,FE)
  elseif HEVISolver == "OneSided"
    return JacOneSided(backend,FT,M,nz,FE)
  else
    error("Unknown HEVISolver $HEVISolver, use HDG, TwoSidedSplit or OneSided")
  end
end
