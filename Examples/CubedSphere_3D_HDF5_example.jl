# Example: Write 3D extruded CubedSphere grid to HDF5 VTK format

using Pkg
Pkg.activate(".")

using CGDycore
using CGDycore.FEM.Outputs.HDF5vtk
using MPI
using KernelAbstractions

# Initialize MPI (required for grid creation)
MPI.Init()

# Create a simple CubedSphere grid with vertical layers
backend = CPU()
FTB = Float64
OrdPoly = 1
nz = 4  # 1 vertical layer with dz[1] = 0.5 (from z[1]=0 to z[2]=0.5)
nPanel = 10  # Full CubedSphere (6 panels)
RefineLevel = 0
ns = 4  # Grid resolution
nLon = 0
nLat = 0
LatB = 0.0
GridType = "CubedSphere"
Decomp = "EqualArea"
RadEarth = 1.0  # Small radius for easier visualization
H = 0.5  # Atmosphere height (50% of radius)

# Create minimal model and parallel structs
Model = CGDycore.DyCore.ModelStruct{FTB}()
ParallelCom = CGDycore.DyCore.ParallelComStruct()
ParallelCom.Proc = 1
ParallelCom.ProcNumber = 1

# Create the grid
Grid, CellToProc = CGDycore.Grids.InitGridSphere(backend, FTB, OrdPoly, nz, nPanel, RefineLevel, ns, nLon, nLat, LatB,
    GridType, Decomp, RadEarth, Model, ParallelCom; Discretization="DG")

# Add vertical grid
CGDycore.Grids.AddVerticalGrid!(Grid, nz, H)

println("Grid created with $(Grid.NumFaces) faces, $(Grid.NumNodes) nodes, nz=$(Grid.nz)")
println("Vertical levels: z = $(Grid.z)")
println("Layer thicknesses: dzeta = $(Grid.dzeta)")

# Create some dummy cell data (one value per face)
cell_data = Dict("face_id" => collect(1:nz*Grid.NumFaces),"face_id1" => collect(1:nz*Grid.NumFaces))

# Write to HDF5 VTK format
HDF5vtk.write_vtk_hdf5("cubed_sphere_3d_grid", Grid, cell_data)

println("3D CubedSphere grid written to cubed_sphere_3d_grid.hdf")

# Clean up
MPI.Finalize()
