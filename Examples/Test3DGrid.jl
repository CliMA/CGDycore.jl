# Test 3D extruded grid HDF5 output

using Pkg
Pkg.activate(".")

using CGDycore
using CGDycore.FEM.Outputs.HDF5vtk
using MPI
using KernelAbstractions

# Initialize MPI
MPI.Init()

function test_3d_grid()
    backend = CPU()
    FTB = Float64
    OrdPoly = 1
    nz = 4  # 4 vertical layers
    nPanel = 1  # Small grid for testing
    RefineLevel = 0
    GridType = "CubedSphere"
    Decomp = "EqualArea"
    RadEarth = 6.371e6
    H = 10000.0  # 10km atmosphere

    Model = CGDycore.DyCore.ModelStruct{FTB}()
    ParallelCom = CGDycore.DyCore.ParallelComStruct()
    ParallelCom.Proc = 1
    ParallelCom.ProcNumber = 1

    println("Testing 3D grid with nz=$nz, nPanel=$nPanel")
    
    Grid, CellToProc = CGDycore.Grids.InitGridSphere(backend, FTB, OrdPoly, nz, nPanel, RefineLevel, 4, 0, 0, 0.0,
        GridType, Decomp, RadEarth, Model, ParallelCom; Discretization="DG")

    # Add vertical grid
    CGDycore.Grids.AddVerticalGrid!(Grid, nz, H)
    
    println("  NumNodes: $(Grid.NumNodes), NumFaces: $(Grid.NumFaces), nz: $(Grid.nz)")
    println("  Vertical levels: z = $(Grid.z)")
    println("  Layer thicknesses: dzeta = $(Grid.dzeta)")
    
    # Create cell data and write
    cell_data = Dict("face_id" => collect(1:Grid.NumFaces))
    filename = "test_3d_grid"
    
    HDF5vtk.write_vtk_hdf5(filename, Grid, cell_data)
    println("  Written to $filename.hdf")
end

test_3d_grid()