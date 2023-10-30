# CGDycore.jl

GDycore is an expermental Julia code [Julia](https://julialang.org/) to study numerical methods and implementational details of numerical dycores for weather prediction. The code is parallelized by MPI and the time integration part of the code is already ported to accelerator by using the Julia package KernelAbstractions.jl. within this programming paradigm the code runs on CPUs and different GPU backends (CUDA, Metal). At the moment a spectral continuous Galerkin method following the HOMME philosphy is implemented as first numerical dycore.

## Code structure and usage
* The code has a general data structure for unstructured grids on the sphere. There exist grid generators for cubed sphere grids, icosahedral triangular and hexagonal grids. In addition there are options to read in other types of unstructured quad grids n the sphere.
* Three dimensional grids (extension to the vertical) are constructed by extruding the two-dimensional horizontal grids.
