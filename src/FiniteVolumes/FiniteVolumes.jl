module FiniteVolumes

import ..Models
import ..Grids
import ..Outputs
import ..Grids
import ..FEMSei

using MPI
using LinearAlgebra
using SparseArrays
using StaticArrays
using KernelAbstractions
using KernelAbstractions: @atomic, @atomicswap, @atomicreplace
using DynamicPolynomials


mutable struct MetricFiniteVolume{FT<:AbstractFloat,
                        AT1<:AbstractArray} 
  PrimalVolume::AT1
  DualVolume::AT1
  PrimalEdge::AT1
  DualEdge::AT1
  DualEdgeVolume::AT1                      
end

include("Divergence.jl")
include("Gradient.jl")
include("Project.jl")
include("MetricFV.jl")
include("MPFA.jl")
include("TangentialRec.jl")

end
