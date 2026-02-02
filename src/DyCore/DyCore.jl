module DyCore

import ..DG
import ..Parallels
import ..Models
import ..Surfaces
import ..Grids
import ..Outputs
import ..Integration
import ..FiniteElements
import ..CGSEM
import ..DGSEM

using MPI
using LinearAlgebra
using SparseArrays
using ArgParse
using UnPack
using StaticArrays
using NetCDF
using NCDatasets
using Statistics
using StrideArraysCore: @gc_preserve, StrideArray, StaticInt
using KernelAbstractions
using KernelAbstractions: @atomic, @atomicswap, @atomicreplace

export parse_commandline 

include("Average.jl")
include("AverageFB.jl")
include("Damping.jl")
include("DiscretizationDG2.jl")
include("Project.jl")
include("ProjectW.jl")
include("ProjectVec.jl")
include("Source.jl")
include("parse_commandline.jl")
include("GlobalVariables.jl")
include("InitDriver.jl")

end
