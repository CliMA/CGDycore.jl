module DGSEM

import ..DG
import ..Parallels
import ..Examples
import ..Models
import ..Surfaces
import ..Grids
import ..Sources
import ..Outputs
import ..Integration
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
include("PresSh.jl")
include("Source.jl")
include("Fcn.jl")
include("DiscretizationDG.jl")
include("TimeStepper.jl")
include("Jac.jl")
include("MIS.jl")
include("MISMethod.jl")
include("FcnLin.jl")
include("Orography.jl")

end
