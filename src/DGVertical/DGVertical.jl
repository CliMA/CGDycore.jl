module DGVertical

import ..DG
import ..Parallels
import ..Models
import ..Surfaces
import ..Grids
import ..Outputs
import ..Integration
import ..FiniteElements

using StaticArrays
using LinearAlgebra
using SparseArrays
using AMD
using KernelAbstractions
using KernelAbstractions: @atomic, @atomicswap, @atomicreplace
using KernelAbstractions.Extras: @unroll


include("LocalFlux.jl")
#include("Rotation.jl")
include("RiemannNonLin.jl")
include("FluxVolumeNonLin.jl")
#include("PresSh.jl")
#include("Source.jl")
include("Fcn.jl")
include("Jac.jl")

end
