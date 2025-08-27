module GPUS

import ..Parallels
import ..Models
import ..Grids
import ..FiniteElements
import ..Thermodynamics
import ..Surfaces

using StaticArrays
using KernelAbstractions
using KernelAbstractions: @atomic, @atomicswap, @atomicreplace
using KernelAbstractions.Extras
using NVTX

include("OperatorKernel.jl")
include("FcnGPU.jl")
include("DiagnosticKernel.jl")
include("InitialConditions.jl")
include("InitialKernel.jl")
include("HorLimiterKernel.jl")
include("dampingGPU.jl")
include("MomentumKernel.jl")
include("ViscKernel.jl")
include("GradKineticKernel.jl")
include("GradPressureKernel.jl")
include("Coriolis.jl")

end
