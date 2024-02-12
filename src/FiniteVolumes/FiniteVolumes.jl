module FiniteVolumes

import ..Models
import ..Grids
import ..Outputs
import ..FEMSei

using MPI
using LinearAlgebra
using SparseArrays
using StaticArrays
using KernelAbstractions
using KernelAbstractions: @atomic, @atomicswap, @atomicreplace


mutable struct MetricFiniteVolume{FT<:AbstractFloat,
                        AT1<:AbstractArray} 
  PrimalVolume::AT1
  DualVolume::AT1
  PrimalEdge::AT1
  DualEdge::AT1
  DualEdgeVolume::AT1                      
end

include("Divergence.jl")
include("Project.jl")
include("MetricFV.jl")

end
