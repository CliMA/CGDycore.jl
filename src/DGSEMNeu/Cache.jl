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
