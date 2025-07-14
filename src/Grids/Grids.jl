module Grids

import ..DG
import ..Parallels
import ..Models

using LinearAlgebra
using RootSolvers
using Dierckx
using NetCDF
using NCDatasets
using PairedLinkedLists
using KernelAbstractions
using KernelAbstractions: @atomic, @atomicswap, @atomicreplace
using CUDA
using StaticArrays    
using NearestNeighbors
using Distances
using LinearAlgebra
using SparseArrays
using Interpolations
using MPI
using FastGaussQuadrature
using DynamicPolynomials

EPS::Float64 = eps(Float64)
EPS2 = EPS * EPS
MAXSIZE::Int = 10

abstract type ElementType end
struct Tri <: ElementType end
struct TriPlanar <: ElementType end
struct Quad <: ElementType end
struct QuadPrimal <: ElementType end
struct QuadDual <: ElementType end
struct Line <: ElementType end
struct Polygonal <: ElementType end


include("point.jl")
include("Node.jl")
include("Edge.jl")
include("Face.jl")
include("geometry_circle.jl")
include("GridStruct.jl")
include("CubedGrid.jl")
include("AddVerticalGrid.jl")
include("CartGrid.jl")
include("CartGridTri.jl")
include("FacesInNodes.jl")
include("Adapt.jl")
include("JacobiDG3.jl")
include("JacobiDG3GPU.jl")
include("JacobiSphere3.jl")
include("JacobiSphere3GPU.jl")
include("JacobiDG2GPU.jl")
include("JacobiSphereDG3GPU.jl")
include("OrientFaceCart.jl")
include("OrientFaceSphere.jl")
include("Orientation.jl")
include("Renumbering.jl")
include("Topo.jl")
include("TopographySmoothing.jl")
include("TopoNeu.jl")
include("NearestNeighbour.jl")
include("vtkWriteHex.jl")
include("Connectivity.jl")
include("SubGrid.jl")
include("Triangular.jl")
include("InputGrid.jl")
include("polygon.jl")
include("intersect.jl")
include("interpolate.jl")
include("InitGrid.jl")
include("SphericalGrid.jl")
include("Grid2KiteGrid.jl")
include("GridLength.jl")
include("QuadGrid.jl")
include("HealpixGrid.jl")
include("TestGrid.jl")
include("OrientTriangle.jl")
include("RefinePoints.jl")

end
