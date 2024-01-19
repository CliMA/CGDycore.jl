module CGDycore

using LinearAlgebra
using SparseArrays
using UnPack
using StructArrays
using StaticArrays
using WriteVTK
using ArgParse
using NetCDF
using Statistics
using StrideArraysCore: @gc_preserve, StrideArray, StaticInt
using RootSolvers
using CUDA
using KernelAbstractions
using KernelAbstractions: @atomic, @atomicswap, @atomicreplace
using PairedLinkedLists


include("DG/DG.jl")
include("Parallels/Parallels.jl")
include("Grids/Grids.jl")
include("Thermodynamics/Thermodynamics.jl")
include("Examples/Examples.jl")
include("Models/Models.jl")
include("Statistics/Statistics.jl")
include("Outputs/Outputs.jl")
include("Integration/Integration.jl")
include("DyCore/DyCore.jl")
include("GPU/GPU.jl")
include("FiniteElements/FiniteElements.jl")
include("Surfaces/Surfaces.jl")

OOP = 5

end
