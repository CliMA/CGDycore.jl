module Outputs

import ..DG
import ..Grids
import ..FiniteElements
import ..Thermodynamics

#using WriteVTK
using KernelAbstractions
using KernelAbstractions: @atomic, @atomicswap, @atomicreplace
using WriteVTK


include("vtkCG.jl")	
include("vtkCGGrid.jl")	
include("vtkSphere.jl")
include("vtkOutputKernel.jl")

end
