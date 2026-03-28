module FiniteElements

import ..DG
import ..Grids

using LinearAlgebra
using FastGaussQuadrature
using KernelAbstractions
using DynamicPolynomials
using KernelAbstractions
using KernelAbstractions: @atomic, @atomicswap, @atomicreplace
using KernelAbstractions.Extras: @unroll

include("ConstructDG.jl")
include("FiniteElement.jl")	
include("NumberingFem.jl")
include("NumberingFemDG.jl")
include("Metric.jl")

end
