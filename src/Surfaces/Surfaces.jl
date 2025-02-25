module Surfaces

import ..Thermodynamics
import ..Grids

using StaticArrays
using KernelAbstractions

include("solver.jl")
include("wrfscheme.jl")
include("functions.jl")
include("surface.jl")
include("BoundaryLayer.jl")

global const TSurfPos = 1
global const RhoVSurfPos = 2
global const uStarPos = 3
global const CMPos = 4
global const CTPos = 5
global const CHPos = 6
global const RiBSurfPos = 7
global const hBLPos = 8
global const LenSurfaceData = 8


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
end

function SurfaceData{FT}(backend,LenSurfaceData,NumG) where FT<:AbstractFloat
  Data = KernelAbstractions.zeros(backend,FT,LenSurfaceData,NumG)
  return SurfaceData{FT,
                     typeof(Data)}(
    Data,                 
  )
end  
end
