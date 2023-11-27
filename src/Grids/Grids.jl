module Grids

import ..DG
import ..Parallels

using LinearAlgebra
using RootSolvers
using Dierckx
using NCDatasets
using PairedLinkedLists
using KernelAbstractions
using KernelAbstractions: @atomic, @atomicswap, @atomicreplace
using CUDA
using StaticArrays


include("Geometry.jl")
include("Node.jl")
include("Edge.jl")
include("Face.jl")
include("GridStruct.jl")
include("CubedGrid.jl")
include("AddVerticalGrid.jl")
include("CartGrid.jl")
include("FacesInNodes.jl")
include("JacobiDG3.jl")
include("JacobiDG3GPU.jl")
include("JacobiSphere3.jl")
include("JacobiSphere3GPU.jl")
include("OrientFaceCart.jl")
include("OrientFaceSphere.jl")
include("Orientation.jl")
include("Renumbering.jl")
include("Topo.jl")
include("TopoNeu.jl")
include("Trans.jl")
include("NearestNeighbour.jl")
include("vtkWriteHex.jl")
include("Connectivity.jl")
include("SubGrid.jl")
include("Triangular.jl")
include("InputGrid.jl")

end
