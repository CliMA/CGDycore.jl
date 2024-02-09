module Surfaces

import ..Thermodynamics
import ..Grids

include("solver.jl")
include("wrfscheme.jl")
include("functions.jl")
include("surface.jl")

mutable struct SurfaceData
  TS::Float64
  RhoVS::Float64
  z0M::Float64
  z0H::Float64
  uStar::Float64
  CT::Float64
  CH::Float64
  LandClass::Int
end
function SurfaceData()
  TS = 300.0
  RhoVS = 0.0
  z0M = 0.0
  z0H = 0.0
  uStar = 0.0
  CT = 0.0
  CH = 0.0
  LandClass = 0
  return SurfaceData(
    TS,
    RhoVS,
    z0M,
    z0H,
    uStar,
    CT,
    CH,
    LandClass,
  )
end  

end
