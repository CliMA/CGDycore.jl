# CGDycore.jl

Welcome to the documentation of **CGDycore.jl** – a flexible, high-performance finite element dynamical core for geophysical fluid dynamics, developed by the CliMA project.

## Overview

CGDycore.jl provides a modular framework for the simulation of shallow water and atmospheric flows on various grids (e.g., cubed sphere, traingular, planar, icosahedral) using continuous and discontinuous Galerkin finite element methods. The code is designed for research and teaching in numerical weather prediction, climate modeling, and computational fluid dynamics.

**Key features:**
- High-order finite element discretizations (CG/DG)
- Support for various grid types (cubed- and triangular sphere, planar, etc.)
- Matrix-free and efficient sparse matrix operators
- Modular structure for easy extension and experimentation
- Parallelization via MPI and GPU support

## Getting Started

### Installation

Clone the repository and install dependencies:
```bash
git clone https://github.com/CliMA/CGDycore.jl.git
cd CGDycore.jl
julia --project -e 'using Pkg; Pkg.instantiate()'
```

Running an Example
You can run a standard test case (e.g., the Bickley Jet) directly from the command line. For example:
```bash
julia --project Jobs/FEM/BickleyJetFEMConsQuad
```
or for the Galewsky test case:
```bash
julia --project Jobs/FEM/GalewskyFEMConsQuad
```
All job scripts in the Jobs/FEM/ directory can be made executable with:
```bash
chmod +x Jobs/FEM/*
```
and then started directly:
```bash
./Jobs/FEM/BickleyJetFEMConsQuad
```

Customizing Simulations
You can adjust simulation parameters (e.g., grid type, order, simulation time) by editing the job scripts or passing command line arguments. See the Examples section for more details and available test cases.

## Documentation Structure

- [Mass Matrix](MassMatrix.md): Theory and implementation of the mass matrix operator.
- [Stiffness Matrix](StiffMatrix.md): Details on the stiffness matrix and related operators.
- [Examples](Examples.md): Overview of available test cases and how to run them.

## Further Information

- [CGDycore.jl on GitHub](https://github.com/CliMA/CGDycore.jl)
- [CliMA Project](https://clima.caltech.edu/)

---

If you are new to Julia, see the [Julia Language documentation](https://docs.julialang.org/) for installation and
