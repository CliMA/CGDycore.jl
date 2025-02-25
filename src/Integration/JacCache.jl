mutable struct JStruct{FT<:AbstractFloat,
                           IT1<:AbstractArray, 
                           AT1<:AbstractArray,
                           AT3<:AbstractArray,
                           AT4<:AbstractArray}
    JRhoW::AT3
    JWRhoTh::AT3
    JWRho::AT3
    JWRhoV::AT3
    JRhoThW::AT3
    JTrW::AT4
    JWW::AT3
    tri::AT3
    JDiff::AT3
    JAdvC::AT3
    JAdvF::AT3
    ListTracer::IT1
    SN::AT1
    C::AT1
    CompTri::Bool
    CompJac::Bool
    CacheCol1::AT1
    CacheCol2::AT1
    CacheCol3::AT1
end

function JStruct{FT}(backend) where FT<:AbstractFloat
  JRhoW=KernelAbstractions.zeros(backend,FT,0,0,0)
  JWRhoTh=KernelAbstractions.zeros(backend,FT,0,0,0)
  JWRho=KernelAbstractions.zeros(backend,FT,0,0,0)
  JWRhoV=KernelAbstractions.zeros(backend,FT,0,0,0)
  JRhoThW=KernelAbstractions.zeros(backend,FT,0,0,0)
  JTrW=KernelAbstractions.zeros(backend,FT,0,0,0,0)
  JWW=KernelAbstractions.zeros(backend,FT,0,0,0)
  tri=KernelAbstractions.zeros(backend,FT,0,0,0)
  JDiff=KernelAbstractions.zeros(backend,FT,0,0,0)
  JAdvC=KernelAbstractions.zeros(backend,FT,0,0,0)
  JAdvF=KernelAbstractions.zeros(backend,FT,0,0,0)
  ListTracer = KernelAbstractions.zeros(backend,Int32,0)
  SN = KernelAbstractions.zeros(backend,FT,0)
  C = KernelAbstractions.zeros(backend,FT,0)
  CompTri=false
  CompJac=false
  CacheCol1=KernelAbstractions.zeros(backend,0)
  CacheCol2=KernelAbstractions.zeros(backend,0)
  CacheCol3=KernelAbstractions.zeros(backend,0)
  return JStruct{FT,
                 typeof(ListTracer),
                 typeof(CacheCol1),
                 typeof(JRhoW),
                 typeof(JTrW)}(
    JRhoW,
    JWRhoTh,
    JWRho,
    JWRhoV,
    JRhoThW,
    JTrW,
    JWW,
    tri,
    JDiff,
    JAdvC,
    JAdvF,
    ListTracer,
    SN,
    C,
    CompTri,
    CompJac,
    CacheCol1,
    CacheCol2,
    CacheCol3,
  )
end  

function JStruct{FT}(backend,NumG,nz,NumTr,TkePos) where FT<:AbstractFloat
  JRhoW=KernelAbstractions.zeros(backend,FT,2,nz-1,NumG)
  JWRhoTh=KernelAbstractions.zeros(backend,FT,2,nz-1,NumG)
  JWRho=KernelAbstractions.zeros(backend,FT,2,nz-1,NumG)
  JWRhoV=KernelAbstractions.zeros(backend,FT,2,nz-1,NumG)
  JRhoThW=KernelAbstractions.zeros(backend,FT,2,nz-1,NumG)
  JTrW=KernelAbstractions.zeros(backend,FT,2,nz-1,NumG,NumTr)
  JWW=KernelAbstractions.zeros(backend,FT,1,nz-1,NumG)
  tri=KernelAbstractions.zeros(backend,FT,3,nz-1,NumG)
  JDiff=KernelAbstractions.zeros(backend,FT,3,nz,NumG)
  JAdvC=KernelAbstractions.zeros(backend,FT,3,nz,NumG)
  JAdvF=KernelAbstractions.zeros(backend,FT,3,nz-1,NumG)
  if TkePos > 0
    ListTracer = KernelAbstractions.zeros(backend,Int32,NumTr+4)
    ListTracer[1] = 2
    ListTracer[2] = 3
    ListTracer[3] = 5
    ListTracer[4] = TkePos
    for iT = 1 : NumTr
      ListTracer[4+iT] = 6 + iT
    end  
    SN = KernelAbstractions.ones(backend,FT,NumTr+4)
    C = KernelAbstractions.zeros(backend,FT,NumTr+4)
  else    
    ListTracer = KernelAbstractions.zeros(backend,Int32,NumTr+3)
    ListTracer[1] = 2
    ListTracer[2] = 3
    ListTracer[3] = 5
    for iT = 1 : NumTr
      ListTracer[3+iT] = 5 + iT
    end  
    SN = KernelAbstractions.ones(backend,FT,NumTr+3)
    C = KernelAbstractions.zeros(backend,FT,NumTr+3)
  end  
  CompTri=false
  CompJac=false
  CacheCol1=KernelAbstractions.zeros(backend,FT,nz)
  CacheCol2=KernelAbstractions.zeros(backend,FT,nz)
  CacheCol3=KernelAbstractions.zeros(backend,FT,nz)
  return JStruct{FT,
                 typeof(ListTracer),
                 typeof(CacheCol1),
                 typeof(JRhoW),
                 typeof(JTrW)}(
    JRhoW,
    JWRhoTh,
    JWRho,
    JWRhoV,
    JRhoThW,
    JTrW,
    JWW,
    tri,
    JDiff,
    JAdvC,
    JAdvF,
    ListTracer,
    SN,
    C,
    CompTri,
    CompJac,
    CacheCol1,
    CacheCol2,
    CacheCol3,
  )
end
