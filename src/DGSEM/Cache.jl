mutable struct CacheAuxStructDG{FT<:AbstractFloat,
                              AT3<:AbstractArray,
                              AT4<:AbstractArray}
  Aux::AT4
  KV::AT3
end
  
function CacheAuxStruct(backend,FT,FE::FiniteElements.DGElement,M,nz,Model,Grid)
  NumG = FE.NumG
  NumI = FE.NumI
  Aux = KernelAbstractions.zeros(backend,FT,M,nz,NumG,Model.NumAux)
  KV = KernelAbstractions.zeros(backend,FT,M,nz,NumI)
  return CacheAuxStructDG{FT,
                        typeof(KV),
                        typeof(Aux)}(
    Aux,
    KV,
  )
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
