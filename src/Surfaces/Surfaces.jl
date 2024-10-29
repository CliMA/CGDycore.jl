module Surfaces

import ..Thermodynamics
import ..Grids

using StaticArrays
using KernelAbstractions

include("solver.jl")
include("wrfscheme.jl")
include("functions.jl")
include("surface.jl")


mutable struct LandUseData{FT<:AbstractFloat,
                           AT1<:AbstractArray,
                           IT1<:AbstractArray}
  z0M::AT1
  z0H::AT1
  LandClass::IT1
end

function LandUseData{FT}(backend,NumG) where FT<:AbstractFloat
  z0M = KernelAbstractions.zeros(backend,FT,NumG)
  z0H = KernelAbstractions.zeros(backend,FT,NumG)
  LandClass = KernelAbstractions.zeros(backend,Int,NumG)
  return LandUseData{FT,
                     typeof(z0M),
                     typeof(LandClass)}(
    z0M,
    z0H,
    LandClass,
  )
end  

mutable struct SurfaceData{FT<:AbstractFloat,
                           AT1<:AbstractArray}
  TS::AT1
  RhoVS::AT1
  uStar::AT1
  CM::AT1
  CT::AT1
  CH::AT1
  RiBS::AT1
  hBL::AT1
end

function SurfaceData{FT}(backend,NumG) where FT<:AbstractFloat
  TS = KernelAbstractions.zeros(backend,FT,NumG)
  RhoVS = KernelAbstractions.zeros(backend,FT,NumG)
  uStar = KernelAbstractions.zeros(backend,FT,NumG)
  CM = KernelAbstractions.zeros(backend,FT,NumG)
  CT = KernelAbstractions.zeros(backend,FT,NumG)
  CH = KernelAbstractions.zeros(backend,FT,NumG)
  RiBS = KernelAbstractions.zeros(backend,FT,NumG)
  hBL = KernelAbstractions.zeros(backend,FT,NumG)
  return SurfaceData{FT,
                     typeof(TS)}(
    TS,
    RhoVS,
    uStar,
    CM,
    CT,
    CH,
    RiBS,
    hBL,
  )
end  
end
