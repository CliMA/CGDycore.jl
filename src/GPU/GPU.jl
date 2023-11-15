module GPU

import ..Parallels
import ..Models
import ..Grids
using StaticArrays
using KernelAbstractions
using KernelAbstractions: @atomic, @atomicswap, @atomicreplace

include("OperatorKernel.jl")
include("FcnGPU.jl")
include("DiagnosticKernel.jl")
include("InitialConditions.jl")
include("InitialKernel.jl")
include("surface.jl")
include("HorLimiterKernel.jl")

end
