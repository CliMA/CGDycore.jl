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
using Polynomials
using SpecialPolynomials

export parse_commandline 

include("Quadratur.jl")
include("MassMatrixVec.jl")
include("RT0.jl")

end