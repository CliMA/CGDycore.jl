module DGSEM

import ..DG: DerivativeX!, DerivativeY!
import ..DG
import ..Parallels
import ..Models
import ..Surfaces
import ..Grids
import ..Outputs
import ..Integration
import ..FiniteElements

using StaticArrays
using KernelAbstractions
using KernelAbstractions: @atomic, @atomicswap, @atomicreplace
using KernelAbstractions.Extras: @unroll


include("RiemannNonLin.jl")
include("FluxVolumeNonLin.jl")
include("PresSh.jl")
include("Source.jl")
include("Fcn.jl")

end
