module GPU

import ..Parallels
import ..Models
import ..Grids
import ..Thermodynamics
using StaticArrays
using KernelAbstractions
using KernelAbstractions: @atomic, @atomicswap, @atomicreplace
using KernelAbstractions.Extras

include("OperatorKernel.jl")
include("FcnGPU.jl")
include("DiagnosticKernel.jl")
include("InitialConditions.jl")
include("InitialKernel.jl")
include("surface.jl")
include("HorLimiterKernel.jl")
include("dampingGPU.jl")

end
