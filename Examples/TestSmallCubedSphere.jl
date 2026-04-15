# Test small CubedSphere grids with different nPanel values

using Pkg
Pkg.activate(".")

using CGDycore
using CGDycore.FEM.Outputs.HDF5vtk
using MPI
using KernelAbstractions

# Initialize MPI
MPI.Init()

function test_grid(nPanel, ns, filename_prefix)
    backend = CPU()
    FTB = Float64
    OrdPoly = 1
    nz = 1
    RefineLevel = 0
    GridType = "CubedSphere"
    Decomp = "EqualArea"
    RadEarth = 6.371e6

    Model = CGDycore.DyCore.ModelStruct{FTB}()
    ParallelCom = CGDycore.DyCore.ParallelComStruct()
    ParallelCom.Proc = 1
    ParallelCom.ProcNumber = 1

    println("\nTesting nPanel=$nPanel, ns=$ns")
    
    Grid, CellToProc = CGDycore.Grids.InitGridSphere(backend, FTB, OrdPoly, nz, nPanel, RefineLevel, ns, 0, 0, 0.0,
        GridType, Decomp, RadEarth, Model, ParallelCom; Discretization="DG")

    println("  NumNodes: $(Grid.NumNodes), NumFaces: $(Grid.NumFaces)")
    
    # Check node indices in a few faces
    println("  First 3 faces node indices:")
    for i in 1:min(3, Grid.NumFaces)
        nodes = Grid.Faces[i].N
        valid = all(n -> n >= 1 && n <= Grid.NumNodes, nodes)
        println("    Face[$i].N = $nodes (valid=$valid)")
    end
    
    # Create cell data and write
    cell_data = Dict("face_id" => collect(1:Grid.NumFaces))
    filename = "$(filename_prefix)_nPanel$(nPanel)_ns$(ns)"
    
    HDF5vtk.write_vtk_hdf5(filename, Grid, cell_data)
    println("  Written to $filename.hdf")
end

# Test with different configurations
test_grid(1, 1, "test_grid")
test_grid(2, 1, "test_grid")
test_grid(2, 2, "test_grid")

println("\nAll tests completed")
