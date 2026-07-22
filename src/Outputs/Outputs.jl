module Outputs

import ..DG
import ..Grids
import ..FiniteElements
import ..Thermodynamics

#using WriteVTK
using StaticArrays
using KernelAbstractions
using KernelAbstractions: @atomic, @atomicswap, @atomicreplace
using KernelAbstractions.Extras: @unroll
using LinearAlgebra
using WriteVTK
using HDF5


include("vtkCG.jl")	
include("vtkCGGrid.jl")	
include("vtkSphere.jl")
include("vtkOutputKernel.jl")
include("HDF5vtk.jl")

end
