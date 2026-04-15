module DGSEM

import ..DG
import ..Parallels
import ..Examples
import ..Models
import ..Surfaces
import ..Grids
import ..Sources
import ..Outputs
import ..FiniteElements

using StaticArrays
using SparseArrays
using LinearAlgebra
using SparseArrays
using BandedMatrices
using FillArrays
using KernelAbstractions
using KernelAbstractions: @atomic, @atomicswap, @atomicreplace
using KernelAbstractions.Extras: @unroll



include("LocalFlux.jl")
include("Rotation.jl")
include("RiemannNonLin.jl")
include("FluxVolumeNonLin.jl")
include("Fcn.jl")
include("Jac.jl")
include("FcnLin.jl")
include("Orography.jl")
include("GeoPot.jl")
include("Cache.jl")
include("Permutation.jl")

end
