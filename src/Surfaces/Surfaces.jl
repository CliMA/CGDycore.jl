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
global const zetaPos = 8
global const LenSurfaceData = 9

global const NumLandClasses = 2
global const SeaClass = 2



mutable struct LandUseData{FT<:AbstractFloat,
                           AT1<:AbstractArray,
                           IT1<:AbstractArray}
  z0M::AT1
  z0H::AT1
  LandClass::IT1
end


function LandUseData{FT}(backend,NumG) where FT<:AbstractFloat
  z0M = KernelAbstractions.zeros(backend,FT,NumLandClasses)
  z0H = KernelAbstractions.zeros(backend,FT,NumLandClasses)
  LandClass = KernelAbstractions.zeros(backend,Int,NumG)
  z0MCPU = zeros(NumLandClasses)
  z0HCPU = zeros(NumLandClasses)
  z0MCPU[1] = 0.01
  z0MCPU[2] = 0.01
  z0HCPU[1] = 0.01
  z0HCPU[2] = 0.01
  copyto!(z0M,z0MCPU)
  copyto!(z0H,z0HCPU)
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

function SurfaceData{FT}(backend,NumG) where FT<:AbstractFloat
  Data = KernelAbstractions.zeros(backend,FT,LenSurfaceData,NumG)
  if NumG > 0
    @. Data[zetaPos,:] = FT(1)
    @. Data[uStarPos,:] = FT(1)
  end
  return SurfaceData{FT,
                     typeof(Data)}(
    Data,                 
  )
end  
end
