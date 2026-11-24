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
global const zetaPos = 9
# wind speed |v| at the lowest level; uStarPos holds the friction velocity
global const UaPos = 10
global const LenSurfaceData = 10

global const NumLandClasses = 2
global const SeaClass = 2



mutable struct LandUseData{FT<:AbstractFloat,
                           AT1<:AbstractArray,
                           IT1<:AbstractArray}
  z0M::AT1
  z0H::AT1
  LandClass::IT1
  # factor of the surface momentum flux per column (1 no-slip, 0 free-slip)
  DragFac::AT1
end


function LandUseData{FT}(backend,NumG;z0M=0.01,z0H=0.01) where FT<:AbstractFloat
  z0MD = KernelAbstractions.zeros(backend,FT,NumLandClasses)
  z0HD = KernelAbstractions.zeros(backend,FT,NumLandClasses)
  LandClass = KernelAbstractions.zeros(backend,Int,NumG)
  DragFac = KernelAbstractions.ones(backend,FT,NumG)
  copyto!(z0MD,fill(FT(z0M),NumLandClasses))
  copyto!(z0HD,fill(FT(z0H),NumLandClasses))
  return LandUseData{FT,
                     typeof(z0MD),
                     typeof(LandClass)}(
    z0MD,
    z0HD,
    LandClass,
    DragFac,
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
