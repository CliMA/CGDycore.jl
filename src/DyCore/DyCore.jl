module DyCore

import ..DG: DerivativeX!, DerivativeY!
import ..DG
import ..Parallels
import ..Models
import ..Surfaces
import ..Grids
import ..Outputs
import ..Integration
import ..FiniteElements

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

export parse_commandline 

include("Average.jl")
include("AverageFB.jl")
include("Damping.jl")
include("DiscretizationCG.jl")
include("DiscretizationDG.jl")
include("Jac.jl")
include("MassCG.jl")
include("Project.jl")
include("ProjectW.jl")
include("ProjectVec.jl")
include("Source.jl")
include("parse_commandline.jl")
include("GlobalVariables.jl")
include("InitDriver.jl")

end
