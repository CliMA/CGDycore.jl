module Surfaces

import ..Thermodynamics
import ..Grids

include("solver.jl")
include("wrfscheme.jl")
include("functions.jl")
include("surface.jl")

struct LandUse{FT<:AbstractFloat}
  thetaS::FT
  qvS::FT
  z0M::FT
  z0H::FT
  LandClass::Int
end

end
