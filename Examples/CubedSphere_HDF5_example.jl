# Example: Write CubedSphere grid to HDF5 VTK format

using Pkg
Pkg.activate(".")

using CGDycore
using CGDycore.FEM.Outputs.HDF5vtk
using MPI
using KernelAbstractions

# Initialize MPI (required for grid creation)
MPI.Init()

# Create a simple CubedSphere grid
backend = CPU()
FTB = Float64
OrdPoly = 1
nz = 1  # Single layer for simplicity
nPanel = 6  # Standard CubedSphere
RefineLevel = 0
ns = 4  # Grid resolution
nLon = 0
nLat = 0
LatB = 0.0
GridType = "CubedSphere"
Decomp = "EqualArea"
RadEarth = 6.371e6

# Create minimal model and parallel structs
Model = CGDycore.DyCore.ModelStruct{FTB}()
ParallelCom = CGDycore.DyCore.ParallelComStruct()
ParallelCom.Proc = 1
ParallelCom.ProcNumber = 1

# Create the grid
Grid, CellToProc = CGDycore.Grids.InitGridSphere(backend, FTB, OrdPoly, nz, nPanel, RefineLevel, ns, nLon, nLat, LatB,
    GridType, Decomp, RadEarth, Model, ParallelCom; Discretization="DG")

println("Grid created with $(Grid.NumFaces) faces and $(Grid.NumNodes) nodes")

# Create some dummy cell data (one value per face)
cell_data = Dict("face_id" => collect(1:Grid.NumFaces))

# Write to HDF5 VTK format
HDF5vtk.write_vtk_hdf5("cubed_sphere_grid", Grid, cell_data)

println("CubedSphere grid written to cubed_sphere_grid.hdf")

# Clean up
MPI.Finalize()
