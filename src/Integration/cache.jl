function CacheJac(backend,FT,M,nz,Model,FE::FiniteElements.DGElement)
  JCache = DGSEM.JacDGVert{FT}(backend,M,nz,FE.NumI)
end
function CacheJac(backend,FT,M,nz,Model,FE::FiniteElements.CGElement)
  JCache = CGSEM.JStruct{FT}(backend,FE.NumG,nz,Model.NumTr,Model.TkePos)
end
function CacheAuxStruct(backend,FT,FE::FiniteElements.CGElement,M,nz,Model,Grid)
  CGSEM.CacheAuxStruct(backend,FT,FE::FiniteElements.CGElement,M,nz,Model,Grid)
end
function CacheAuxStruct(backend,FT,FE::FiniteElements.DGElement,M,nz,Model,Grid)
  DGSEM.CacheAuxStruct(backend,FT,FE::FiniteElements.DGElement,M,nz,Model,Grid)
end
