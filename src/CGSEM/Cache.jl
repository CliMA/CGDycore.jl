mutable struct CacheAuxStructCG{FT<:AbstractFloat,
                              AT3<:AbstractArray,
                              AT4<:AbstractArray}
  Aux::AT4
  KV::AT3
  Temp1::AT4
  q::AT4
end

function CacheAuxStruct(backend,FT,FE::FiniteElements.CGElement,M,nz,Model,Grid)
  NumG = FE.NumG
  NumI = FE.NumI
  Aux = KernelAbstractions.zeros(backend,FT,M,nz,NumG,Model.NumAux)
  KV = KernelAbstractions.zeros(backend,FT,M,nz,NumI)
  lenghthTemp1 = 5 + 1 + 1 + 1 + Model.NumTr + Model.NDEDMF*(1 + 1 +  1 + Model.NumTr)
  Temp1 = KernelAbstractions.zeros(backend,FT,M,nz,NumG,lenghthTemp1)
  q = KernelAbstractions.zeros(backend,FT,M,nz,Grid.NumFaces+Grid.NumFacesG,2*(Model.NumTr+1))
  return CacheAuxStructCG{FT,
                        typeof(KV),
                        typeof(Aux)}(
    Aux,
    KV,
    Temp1,
    q,
  )
end
