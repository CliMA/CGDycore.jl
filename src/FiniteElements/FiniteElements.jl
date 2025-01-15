module FiniteElements

import ..DG
import ..Grids

using LinearAlgebra
using FastGaussQuadrature
using KernelAbstractions

include("FiniteElement.jl")	
include("NumberingFem.jl")
include("NumberingFemDG.jl")

end
