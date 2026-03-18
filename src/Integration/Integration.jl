module Integration

import ..Grids
import ..FiniteElements
import ..Outputs
import ..Parallels
import ..Statistics
import ..Models
import ..CGSEM
import ..DGSEM

using KernelAbstractions
using KernelAbstractions: @atomic, @atomicswap, @atomicreplace
using KernelAbstractions.Extras: @unroll
using StaticArrays
using NVTX

abstract type IntegrationMethod end
mutable struct NoMethod <: IntegrationMethod end

include("IntegrationMethods.jl")
include("LinIMEXMethod.jl")
include("IMEXDirkMethod.jl")
include("SSPRungeKuttaMethod.jl")
include("MISMethod.jl")
include("RosenbrockMethod.jl")
include("RungeKuttaMethod.jl")
include("RungeKuttaIMEXMethod.jl")
include("LSRungeKuttaMethod.jl")
include("Rosenbrock.jl")
include("LinIMEX.jl")
include("IMEXDirk.jl")
include("RungeKuttaExplicit.jl")
include("SSPRungeKutta.jl")
include("LSRungeKutta.jl")
include("MIS.jl")
include("MISSemi.jl")
include("MISLin.jl")
include("Solve.jl")
include("TimeStepper.jl")
include("cache.jl")
include("FastLinAlg.jl")

end
