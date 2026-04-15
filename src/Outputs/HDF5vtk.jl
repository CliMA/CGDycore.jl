module HDF5vtk

using HDF5
using Printf
using LinearAlgebra  # For norm and cross
using ..Grids  # Import Grids module to access GridStruct

export write_vtk_hdf5, compute_offsets

"""
    compute_offsets(cell_sizes)

Compute the offsets array for VTK HDF5 format from cell sizes.

# Arguments
- `cell_sizes`: Array of integers representing the number of nodes per cell

# Returns
- Array of offsets starting with 0
"""
function compute_offsets(cell_sizes::AbstractArray{<:Integer})
    offsets = zeros(Int, length(cell_sizes) + 1)
    for i in 1:length(cell_sizes)
        offsets[i+1] = offsets[i] + cell_sizes[i]
    end
    return offsets
end

"""
    write_vtk_hdf5(filename, points, connectivity, types, offsets, cell_data; point_data=nothing)

Write VTK UnstructuredGrid data to HDF5 format.

# Arguments
- `filename`: Output filename (without extension, .hdf will be added)
- `points`: 3xNpoints array of point coordinates
- `connectivity`: Connectivity array (flattened)
- `types`: Cell types array
- `offsets`: Offset array for connectivity
- `cell_data`: Dictionary of cell data arrays
- `point_data`: Dictionary of point data arrays (optional)
"""
function write_vtk_hdf5(filename::String, points::AbstractArray, connectivity::AbstractArray, 
                       types::AbstractArray, offsets::AbstractArray, cell_data::Dict; 
                       point_data::Union{Dict,Nothing}=nothing)
    
    # Validate inputs
    ncells = length(types)
    if length(offsets) != ncells + 1
        error("Offsets array must have length $(ncells + 1), got $(length(offsets))")
    end
    if offsets[1] != 0
        error("Offsets array must start with 0, got $(offsets[1])")
    end
    if size(points, 1) != 3
        error("Points array must have 3 rows (x,y,z), got $(size(points, 1))")
    end
    
    h5open("$filename.hdf", "w") do file
        VTKHDF = create_group(file, "VTKHDF")
        HDF5.attributes(VTKHDF)["Version"] = [2, 0]
        
        # Set type attribute
        HDF5.attributes(VTKHDF)["Type"] = "UnstructuredGrid"
        
        # Write metadata
        VTKHDF["NumberOfConnectivityIds"] = [length(connectivity)]
        VTKHDF["NumberOfPoints"] = [size(points, 2)]
        VTKHDF["NumberOfCells"] = [length(types)]
        
        # Write geometry
        # Reshape points from (3, N) to (N, 3) for VTK HDF5 format
        npts = size(points, 2)
        points_reshaped = reshape(points, 3, npts)  # Ensure it's (3, npts)
        points_vtk = zeros(npts, 3)
        for i in 1:npts
            points_vtk[i, 1] = points_reshaped[1, i]
            points_vtk[i, 2] = points_reshaped[2, i]
            points_vtk[i, 3] = points_reshaped[3, i]
        end
        VTKHDF["Points"] = Matrix(transpose(points_vtk))
        
        VTKHDF["Types"] = UInt8.(types)  # Convert to unsigned char for VTK compatibility
        VTKHDF["Connectivity"] = connectivity  # Already converted to 0-based indexing above
        VTKHDF["Offsets"] = offsets
        
        # Write cell data
        if !isempty(cell_data)
            CellData = create_group(VTKHDF, "CellData")
            for (name, data) in cell_data
                CellData[name] = data
            end
        end
        
        # Write point data
        if point_data !== nothing && !isempty(point_data)
            PointData = create_group(VTKHDF, "PointData")
            for (name, data) in point_data
                PointData[name] = data
            end
        end
    end
end

"""
    write_vtk_hdf5(filename::String, grid::Grids.GridStruct, cell_data::Dict; point_data::Union{Dict,Nothing}=nothing)

Write a CGDycore GridStruct to HDF5 VTK format.

# Arguments
- `filename`: Output filename (without extension, .hdf will be added)
- `grid`: CGDycore GridStruct containing nodes and faces
- `cell_data`: Dictionary of cell data arrays (one value per face)
- `point_data`: Dictionary of point data arrays (optional)
"""
function write_vtk_hdf5(filename::String, grid::Grids.GridStruct, cell_data::Dict; 
                       point_data::Union{Dict,Nothing}=nothing)
    
    if grid.nz == 0
        # 2D surface grid (no vertical layers)
        # Convert to VTK format: (npoints, 3) with contiguous memory
        npoints = grid.NumNodes
        points_vtk = zeros(3, npoints)
        for i in 1:npoints
            points_vtk[1, i] = grid.Nodes[i].P.x
            points_vtk[2, i] = grid.Nodes[i].P.y  
            points_vtk[3, i] = grid.Nodes[i].P.z
        end
        
        # Build connectivity and types
        ncells = grid.NumFaces
        connectivity = Int[]
        types = UInt8[]
        offsets = zeros(Int, ncells + 1)
        offsets[1] = 0
        
        for iF in 1:ncells
            face = grid.Faces[iF]
            nverts = length(face.N)
            
            # Get vertex coordinates
            verts = zeros(3, nverts)
            for i in 1:nverts
                node_idx = face.N[i]
                verts[1, i] = grid.Nodes[node_idx].P.x
                verts[2, i] = grid.Nodes[node_idx].P.y
                verts[3, i] = grid.Nodes[node_idx].P.z
            end
            
            # Compute face center
            face_center = zeros(3)
            for i in 1:nverts
                face_center[1] += verts[1, i]
                face_center[2] += verts[2, i]
                face_center[3] += verts[3, i]
            end
            face_center ./= nverts
            face_center ./= norm(face_center)  # Normalize for spherical grid
            
            # Compute face normal using cross products (right-hand rule)
            normal = zeros(3)
            for i in 1:nverts-1
                v1 = verts[:, i]
                v2 = verts[:, i+1]
                normal += cross(v1, v2)
            end
            v1 = verts[:, nverts]
            v2 = verts[:, 1]
            normal += cross(v1, v2)
            normal ./= norm(normal)
            
            # Check orientation: normal should point in same direction as face_center
            orientation_correct = dot(normal, face_center) > 0
            
            # Add node indices in correct order
            verts_to_add = if orientation_correct
                face.N
            else
                # Reverse order for correct VTK orientation
                reverse(face.N)
            end
            
            append!(connectivity, verts_to_add)
            
            # Update offset
            offsets[iF+1] = offsets[iF] + nverts
            
            # Determine cell type based on number of vertices
            if nverts == 3
                push!(types, UInt8(5))  # VTK_TRIANGLE
            elseif nverts == 4
                push!(types, UInt8(9))  # VTK_QUAD
            else
                push!(types, UInt8(7))  # VTK_POLYGON
            end
        end
        
        # Convert connectivity to 0-based
        connectivity = connectivity .- 1
        
        # Write the file
        write_vtk_hdf5(filename, points_vtk, connectivity, types, offsets, cell_data; point_data=point_data)
    else
        # 3D extruded grid (nz >= 1)
        nz = grid.nz
        npoints_2d = grid.NumNodes
        npoints_3d = npoints_2d * (nz + 1)  # interfaces
        
        # Create 3D points by extruding radially
        points = zeros(3, npoints_3d)
        point_idx = 1
        
        for iz in 1:(nz + 1)
            height = grid.z[iz]  # z[1] = 0 (surface), z[nz+1] = H (top)
            radius = grid.Rad + height
            
            for i in 1:npoints_2d
                # Get surface point (normalized direction)
                surf_x = grid.Nodes[i].P.x
                surf_y = grid.Nodes[i].P.y
                surf_z = grid.Nodes[i].P.z
                
                # Normalize to get direction vector
                r_surf = sqrt(surf_x^2 + surf_y^2 + surf_z^2)
                dir_x = surf_x / r_surf
                dir_y = surf_y / r_surf
                dir_z = surf_z / r_surf
                
                # Extrude radially outward
                points[1, point_idx] = dir_x * radius
                points[2, point_idx] = dir_y * radius
                points[3, point_idx] = dir_z * radius
                point_idx += 1
            end
        end
        
        # Build 3D connectivity and types
        ncells_2d = grid.NumFaces
        ncells_3d = ncells_2d * nz  # nz 3D cell per face (not per layer)
        connectivity = Int[]
        types = UInt8[]
        offsets = zeros(Int, ncells_3d + 1)
        offsets[1] = 0
        
        cell_idx = 1
        # Create only one layer of 3D cells (from bottom to top interface)
        for iz in 1:nz
        for iF in 1:ncells_2d
            face = grid.Faces[iF]
            nverts = length(face.N)
            
            # Bottom face nodes (surface level, iz=1)
            bottom_nodes = face.N .+ (iz - 1) * npoints_2d
            # Top face nodes (top level, iz=nz+1)  
            top_nodes = face.N .+ iz * npoints_2d
            
            # For radial extrusion on CubedSphere, use consistent orientation
            # VTK_HEXAHEDRON connectivity: bottom face CCW, top face CCW (same order)
            for i in 1:nverts
                push!(connectivity, bottom_nodes[i] - 1)  # 0-based
            end
            # Top face (same order as bottom for radial extrusion)
            for i in 1:nverts
                push!(connectivity, top_nodes[i] - 1)  # 0-based
            end
                
                # Update offset
                offsets[cell_idx + 1] = offsets[cell_idx] + 2 * nverts
                
                # Determine 3D cell type
                if nverts == 3
                    push!(types, UInt8(13))  # VTK_WEDGE
                elseif nverts == 4
                    push!(types, UInt8(12))  # VTK_HEXAHEDRON
                else
                    # For polygons, we'd need VTK_POLYHEDRON but that's more complex
                    push!(types, UInt8(42))  # VTK_POLYHEDRON (placeholder)
                end
                
                cell_idx += 1
            end
        end
        end
        
        # Write the 3D file
        write_vtk_hdf5(filename, points, connectivity, types, offsets, cell_data; point_data=point_data)
    end
end
