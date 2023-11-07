# CGDycore.jl

GDycore is an expermental Julia code [Julia](https://julialang.org/) to study numerical methods and implementational details of numerical dycores for weather prediction. The code is parallelized by MPI and the time integration part of the code is already ported to accelerator by using the Julia package KernelAbstractions.jl. within this programming paradigm the code runs on CPUs and different GPU backends (CUDA, Metal). At the moment a spectral continuous Galerkin method following the HOMME philosphy is implemented as first numerical dycore.

## Code structure and usage
* The code has a general data structure for unstructured grids on the sphere. There exist grid generators for cubed sphere grids, icosahedral triangular and hexagonal grids. In addition there are options to read in other types of unstructured quad grids n the sphere.
* Three dimensional grids (extension to the vertical) are constructed by extruding the two-dimensional horizontal grids.
* The code is parallelized by two-dimensional domain decomposition of the spherical grid. The decomposition is based either on space filling curves or an equa-area decomposition of the sphere in two polar caps and further spherical bands. 
* AbstractionKernels.jl is used to run the code solely on CPU's or also partly on GPU's of different vendors. 
* For time integration different explicit and partially implicit methods are available. Our favoured scheme is the  Rosenbrock-W method.
