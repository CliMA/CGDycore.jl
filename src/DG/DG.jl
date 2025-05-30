module DG

using MuladdMacro
using LinearAlgebra
using FastGaussQuadrature
using DynamicPolynomials

include("DLagrange.jl")
include("DerivativeMatrixSingle.jl")
include("Lagrange.jl")
include("Tools.jl")

end
