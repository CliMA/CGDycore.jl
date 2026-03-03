module IMEXRosenbrock

import ..Integration as INT

using Plots
using Optim
using LinearAlgebra

#include("RungeKuttaMethod.jl")
#include("RosenbrockMethod.jl")
include("IMEXDirkMethod.jl")
include("SSPMethod.jl")
include("StabilityRegion.jl")
include("OrderConditions.jl")
include("residual.jl")
include("FindRosenbrockMethod.jl")
include("RosenbrockToIMEXDirk.jl")
include("IMEXDirkToRosenbrock.jl")
end
