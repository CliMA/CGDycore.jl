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
                           AT2<:AbstractArray}
  Data::AT2                         
  TS::Int
  RhoVS::Int
  uStar::Int
  CM::Int
  CT::Int
  CH::Int
  RiBSurf::Int
  hBL::Int
end

function SurfaceData{FT}(backend,NumG) where FT<:AbstractFloat
  Data = KernelAbstractions.zeros(backend,FT,8,NumG)
  TS = 1
  RhoVS = 2
  uStar = 3
  CM = 4
  CT = 5
  CH = 6
  RiBSurf = 7
  hBL = 8
  return SurfaceData{FT,
                     typeof(Data)}(
    Data,                 
    TS,
    RhoVS,
    uStar,
    CM,
    CT,
    CH,
    RiBSurf,
    hBL,
  )
end  
end
