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

function LandUseData{FT}(backend,DoF,NF) where FT<:AbstractFloat
  z0M = KernelAbstractions.zeros(backend,FT,DoF,NF)
  z0H = KernelAbstractions.zeros(backend,FT,DoF,NF)
  LandClass = KernelAbstractions.zeros(backend,Int,DoF,NF)
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
  CT::AT1
  CH::AT1
end

function SurfaceData{FT}(backend,DoF,NF) where FT<:AbstractFloat
  TS = KernelAbstractions.zeros(backend,FT,DoF,NF)
  RhoVS = KernelAbstractions.zeros(backend,FT,DoF,NF)
  uStar = KernelAbstractions.zeros(backend,FT,DoF,NF)
  CT = KernelAbstractions.zeros(backend,FT,DoF,NF)
  CH = KernelAbstractions.zeros(backend,FT,DoF,NF)
  return SurfaceData{FT,
                     typeof(TS)}(
    TS,
    RhoVS,
    uStar,
    CT,
    CH,
  )
end  
end
