module FEMSei

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
using FastGaussQuadrature

export parse_commandline 

abstract type FiniteElement end
abstract type ScalarElement <: FiniteElement end
abstract type ScalarKitePElement <: ScalarElement end
abstract type VectorElement <: FiniteElement end
abstract type VectorKiteElement <: FiniteElement end
abstract type HDivElement <: VectorElement end
abstract type HDivKiteDElement <: HDivElement end
abstract type HDivConfElement <: HDivElement end
abstract type HCurlElement <: VectorElement end
abstract type HCurlKiteDElement <: HCurlElement end
abstract type HCurlConfElement <: HCurlElement end
abstract type TensorElement <: VectorElement end

include("Quadratur.jl")
include("RT.jl")
include("ND.jl")
include("DG.jl")
include("CG.jl")
include("VecDG.jl")
include("BDM.jl")
include("CGKite.jl")
include("MassMatrix.jl")
include("Jacobi.jl")
include("StiffMatrix.jl")
include("Project.jl")
include("ConvertVelocity.jl")
include("ModelFEM.jl")
include("Fcn.jl")
include("TimestepperFEM.jl")
include("Interpolate.jl")



end
