# Example usage of HDF5vtk module

using Pkg
Pkg.activate(".")

using CGDycore.FEM.Outputs.HDF5vtk

# Example data for a simple triangle mesh
points = [0.0 1.0 0.5;  # x coordinates
          0.0 0.0 1.0;  # y coordinates
          0.0 0.0 0.0]  # z coordinates

connectivity = [1, 2, 3]  # single triangle
types = [5]  # VTK_TRIANGLE = 5

# Compute offsets automatically
cell_sizes = [3]  # one triangle with 3 nodes
offsets = HDF5vtk.compute_offsets(cell_sizes)
println("Computed offsets: ", offsets)

cell_data = Dict("pressure" => [1.0])

println("About to write HDF5 file...")
# Write to HDF5 VTK format
HDF5vtk.write_vtk_hdf5("example_triangle", points, connectivity, types, offsets, cell_data)

println("HDF5 VTK file written: example_triangle.hdf")