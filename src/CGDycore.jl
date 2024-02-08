module CGDycore

using LinearAlgebra
using SparseArrays
using UnPack
using StructArrays
using StaticArrays
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
include("Thermodynamics/Thermodynamics.jl")
include("Models/Models.jl")
include("Grids/Grids.jl")
include("Surfaces/Surfaces.jl")
include("Examples/Examples.jl")
include("Statistics/Statistics.jl")
include("Outputs/Outputs.jl")
include("Integration/Integration.jl")
include("DyCore/DyCore.jl")
include("GPU/GPU.jl")
include("FiniteElements/FiniteElements.jl")
include("FEMSei/FEMSei.jl")

OOP = 5

end
