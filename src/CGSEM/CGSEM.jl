module CGSEM

import ..Parameters as P
import ..Parallels
import ..Examples
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
include("Fcn.jl")
include("Jac.jl")
include("MassCG.jl")
include("DiagnosticKernel.jl")
include("InitialConditions.jl")
include("InitialKernel.jl")
include("HorLimiterKernel.jl")
include("MomentumKernel.jl")
include("ViscKernel.jl")
include("GradKineticKernel.jl")
include("GradPressureKernel.jl")
include("DiscretizationCG.jl")

end
