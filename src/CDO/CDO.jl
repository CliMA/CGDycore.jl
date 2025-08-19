module CDO

import ..Models
import ..Grids
import ..Outputs
import ..Grids
import ..FEM

using MPI
using LinearAlgebra
using SparseArrays
using StaticArrays
using KernelAbstractions
using KernelAbstractions: @atomic, @atomicswap, @atomicreplace
using DynamicPolynomials
using FastGaussQuadrature
using LinearAlgebra


mutable struct GridCDO{FT<:AbstractFloat,
                       AS<:AbstractArray} 
  IncEN::AS
  IncFE::AS
  IncCF::AS
end

mutable struct MetricCDO{FT<:AbstractFloat,
                        AT1<:AbstractArray, 
                        AT2<:AbstractArray} 
  gV::AT1
  gE::AT2
  gF::AT2
  gC::AT1
end
#include("MetricCDO.jl")
include("GridCDO.jl")
end
