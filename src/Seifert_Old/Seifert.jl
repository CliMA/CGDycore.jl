module Seifert

import ..DG
import ..Parallels
import ..Models
import ..Grids
import ..Outputs
import ..Integration

using MPI
using LinearAlgebra
using SparseArrays
using ArgParse
using UnPack
using StaticArrays
using Statistics
using StrideArraysCore: @gc_preserve, StrideArray, StaticInt
using KernelAbstractions
using KernelAbstractions: @atomic, @atomicswap, @atomicreplace
using DynamicPolynomials

export parse_commandline 

abstract type FiniteElement end
abstract type ScalarElement <: FiniteElement end
abstract type VectorElement <: FiniteElement end
abstract type HDivElement <: VectorElement end
abstract type TensorElement <: VectorElement end

include("Quadratur.jl")
include("RT0.jl")
include("DG0.jl")
include("MassMatrix.jl")
include("Jacobi.jl")



end
