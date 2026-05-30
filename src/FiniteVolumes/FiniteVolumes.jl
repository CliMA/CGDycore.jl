module FiniteVolumes

import ..DG
import ..Models
import ..Grids
import ..Outputs
import ..Grids
import ..FiniteElements
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


mutable struct MetricFiniteVolume3{FT<:AbstractFloat,
                        AT2<:AbstractArray,
                        AT3<:AbstractArray}
  PrimalVolume::AT2
  PrimalMidPoints::AT3
  PrimalSideFaces::AT2
  PrimalTopFaces::AT2
  PrimalSideNormals::AT3
  PrimalTopNormals::AT3
  PrimalSideMidPoints::AT3
  PrimalTopMidPoints::AT3
  PrimalPoints::AT3
end

mutable struct MetricFiniteVolume{FT<:AbstractFloat,
                        AT1<:AbstractArray, 
                        AT2<:AbstractArray, 
                        AT3<:AbstractArray} 
  PrimalVolume::AT1
  DualVolume::AT1
  PrimalEdge::AT1
  DualEdge::AT1
  DualEdgeVolume::AT3                      
  PrimalNormal::AT2
end

mutable struct CacheFV
  CurlUu::Array{Float64, 1} 
  TangUu::Array{Float64, 1}
  hE::Array{Float64, 1}
  K::Array{Float64, 1}
  Grad
  Inter
  Div
  Curl
  Tang
# TangV::Array{Float64, 2}
end

include("Divergence.jl")
include("Gradient.jl")
include("Project.jl")
include("MetricFV.jl")
include("MetricFV3.jl")
include("MPFA.jl")
include("TangentialRec.jl")
include("Curl.jl")
include("Advection.jl")
include("FcnFV.jl")
include("AdvectionUpwind.jl")
include("ComputeVolume.jl")

end
