mutable struct MetricStruct{FT<:AbstractFloat,
                            AT2<:AbstractArray,
                            AT3<:AbstractArray,
                            AT4<:AbstractArray,
                            AT5<:AbstractArray,
                            AT6<:AbstractArray}
  J::AT4
  X::AT5
  dXdxI::AT6
  dXdx::AT6
  Rotate::AT6
  nSS::AT2
  nS::AT3
  FS::AT2
  dz::AT2
  zP::AT2
  JC::AT3
  JCW::AT3
  xS::AT2
  VolSurfH::AT4
  NH::AT5
  VolSurfV::AT3
  NV::AT4
end
function MetricStruct{FT}(backend,nQuad,OPZ,NF,nz,NumG) where FT<:AbstractFloat
    J      = KernelAbstractions.zeros(backend,FT,nQuad,OPZ,nz,NF)
    X      = KernelAbstractions.zeros(backend,FT,nQuad,OPZ,3,nz,NF)
    dXdxI  = KernelAbstractions.zeros(backend,FT,3,3,OPZ,nQuad,nz,NF)
    dXdx   = KernelAbstractions.zeros(backend,FT,3,3,OPZ,nQuad,nz,NF)
    Rotate  = KernelAbstractions.zeros(backend,FT,3,3,OPZ,nQuad,nz,NF)
    nSS  = KernelAbstractions.zeros(backend,FT,3,NumG)
    nS = KernelAbstractions.zeros(backend,FT,nQuad,3,NF)
    FS = KernelAbstractions.zeros(backend,FT,nQuad,NF)
    dz = KernelAbstractions.zeros(backend,FT,0,0)
    zP = KernelAbstractions.zeros(backend,FT,0,0)
    JC     = KernelAbstractions.zeros(backend,FT,0,0,0)
    JCW    = KernelAbstractions.zeros(backend,FT,0,0,0)
    xS    = KernelAbstractions.zeros(backend,FT,2,NumG)
    VolSurfH = KernelAbstractions.zeros(backend,FT,0,0,0,0)
    NH = KernelAbstractions.zeros(backend,FT,0,0,0,0,0)
    VolSurfV = KernelAbstractions.zeros(backend,FT,0,0,0)
    NV = KernelAbstractions.zeros(backend,FT,0,0,0,0)
    return MetricStruct{FT,
                        typeof(zP),
                        typeof(nS),
                        typeof(J),
                        typeof(X),
                        typeof(dXdxI)}(
        J,
        X,
        dXdxI,
        dXdx,
        Rotate,
        nSS,
        nS, 
        FS, 
        dz,
        zP,
        JC,
        JCW,
        xS,
        VolSurfH,
        NH,
        VolSurfV,
        NV,
    )
end
