<div align="center">
<img src="https://github.com/user-attachments/assets/bf53eeb8-d610-4114-98a8-de3bc9067e69" width="250"/>
</div>

GDycore is an experimental [Julia code](https://julialang.org/) for studying numerical methods and implementing the details of numerical dycores for weather prediction. The code is parallelised by MPI and the time integration part of the code is already ported to an accelerator by using the Julia package KernelAbstractions.jl. Within this programming paradigm, the code runs on CPUs and different GPU backends (CUDA, Metal). Currently, a spectral continuous Galerkin method, following the HOMME philosophy, is implemented as the first numerical dycore.

## Code structure and usage
* The code has a general data structure for unstructured grids on the sphere. There exist grid generators for cubed sphere grids, icosahedral triangular and hexagonal grids. In addition, there are options to read in other types of unstructured quad grids on the sphere.
* Three-dimensional grids (extension to the vertical) are constructed by extruding the two-dimensional horizontal grids.
* The code is parallelised by two-dimensional domain decomposition of the spherical grid. The decomposition is based either on space-filling curves or an equi-area decomposition of the sphere in two polar caps and further spherical bands. 
* AbstractionKernels.jl is used to run the code solely on CPU's or also partly on GPU's of different vendors. 
* For time integration, different explicit and partially implicit methods are available. Our favoured scheme is the  Rosenbrock-W method.

## CGDycore.jl/FEM.jl

In addition, you will find an experimental environment based on a finite element method for calculating test examples such as Galewsky, Haurwitz, Modon Collison on the sphere. The Bickley jet can also be computed on a planar grid by appropriately adjusting the manifold’s metric terms. An overview of all these options can be found in the table below. This Code is currently not parallelised for CPU or GPU. 

### ✅ Grid Type and Possible Test Cases

| **Nonlinear SWE**       | **Output/Grid** | **Galewsky** (Sphere) | **Haurwitz** (Sphere) | **Modon** (Sphere) | **Galewsky** (Flat) | **Haurwitz** (Flat) | **Modon** (Flat) | **Bickley** (Planar) |
|------------------------|------------------|------------------------|------------------------|---------------------|----------------------|-----------------------|---------------------|---------------------------|
| **Conservative**       | Tri              | ✓                      |   ✓                    |  ✓                  |      ✓               |    ✓                  |    ✓                 |    ✓                       |✓
|                        | Quad             | ✓                      | ✓                      | ✓                 |     ✓                | ✓                     |    ✓                  |      ✓                     |✓
| **Vector Invariant**   | Tri              | ✓                      | ✓                      |  ✓                  |   ✓                  |   ✓                   | ✓                 |   ✓                        |✓
|                        | Quad             | ✓                      |      ✓                 | ✓                |  ✓                   |  ✓                    | ✓                 |   ✓                        |✓

Our main setting is calculating with finite elements order 1, grid-output refining level 1. 
### 📌 Legend:
- **✓** = fully supported and tested 
- **(✓)** = partially tested  
- Terms like "Sphere", "Flat", and "Cartesian" indicate the geometric configuration. Flat output here refers to reprinting the sphere back onto the plane.

### 🚀 Installation Guide

This guide explains how to install [Julia](https://julialang.org/) and set up the `CGDycore.jl` repository.

#### 🧰 Prerequisites

- Git
- Julia ≥ 1.9
- Linux, macOS, or Windows

---

#### 🟣 1. Install Julia

Download and install Julia from the official site:

https://julialang.org/downloads/

After installation, verify it works:
```
julia --version
```
#### 📂 2. Clone the CGDycore.jl repository from GitHub:
```
git clone https://github.com/CliMA/CGDycore.jl.git
cd CGDycore.jl
```
#### 📦 3. Set Up the Julia Environment
Start Julia in the project directory:
```
julia
```
Then activate the environment and install all dependencies:
```
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```
Exit Julia by pressing Ctrl+D or:
```
exit()
```
#### ▶️ Usage
To run a simulation, use one of the example scripts provided in the Jobs/ folder. For example:
```
./Jobs/FEM/GalewskyFEMVecIQuad
```
You can modify or create new configurations based on these examples.

#### :cyclone: Example Simulations
<img src="https://github.com/user-attachments/assets/e56bc457-d1a3-4255-addd-fc00465707d7" width="500"/>

<sub>Figure: Galewsky barotropic instability test after 6 days of simulation on a triangular grid with refinement level 6.</sub>


## 📂 Project Structure
- src/ – Core implementation of the dynamical core
- BatchScripts/ – Example configurations and runs for the /Jobs
- Jobs/ - Examples with all parameters 
- Project.toml – Dependency specification
- Manifest.toml – Exact versions of dependencies
 
## 🤝 Contributing

Contributions are welcome! Feel free to:
- Open issues for bugs or questions
- Suggest improvements
- Submit pull requests
For major changes, please open an issue first to discuss your ideas.
