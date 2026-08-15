using LinearAlgebra

# Simple container to hold both topology and geometry
struct CubedSphereMesh
    plex::JuliaDMPlex
    coordinates::Vector{Vector{Float64}} # Vertex ID -> [x, y, z] on sphere
end

function create_cubed_sphere_mesh(N::Int)
    # 1. Calculate global total counts for sizing arrays
    num_vertices = 6 * N^2 + 2
    num_edges = 12 * N^2 + 12 * N
    num_cells = 6 * N^2
    total_points = num_vertices + num_edges + num_cells
    
    # 2. Offset blocks for layout
    v_offset = 0
    e_offset = num_vertices
    c_offset = num_vertices + num_edges

    # Prepare DMPlex array sizes
    # Vertices have 0 downward elements. Edges have 2. Quad Cells have 4.
    cone_counts = zeros(Int, total_points)
    for e in (v_offset + num_vertices + 1):(v_offset + num_vertices + num_edges)
        cone_counts[e] = 2 # Every edge connects 2 vertices
    end
    for c in (c_offset + 1):total_points
        cone_counts[c] = 4 # Every cell is a quad with 4 edges
    end

    # Generate cumulative offsets
    cone_offsets = zeros(Int, total_points + 1)
    cone_offsets[1] = 1
    for i in 1:total_points
        cone_offsets[i+1] = cone_offsets[i] + cone_counts[i]
    end

    # Prepare the target data holder
    cone_data = Vector{Int}(undef, cone_offsets[end] - 1)
    
    # Track spatial position coordinates for our geometry layout
    coords = Vector{Vector{Float64}}(undef, num_vertices)

    # 3. Stitching Engine (Conceptual loop)
    # In practice, you map out an adjacency table for the 6 faces:
    # Face 1 (+X), Face 2 (-X), Face 3 (+Y), Face 4 (-Y), Face 5 (+Z), Face 6 (-Z)
    
    # Helper index trackers
    # You would implement helper functions that match internal mesh locations:
    # get_global_vertex_id(face, i, j, N) -> returns stable 1..num_vertices ID
    # get_global_edge_id(face, i, j, direction, N) -> returns stable unique edge ID
    
    # Example filling loop for cells to edges:
    cell_idx = 1
    for face in 1:6
        for j in 1:N, i in 1:N
            global_cell = c_offset + cell_idx
            
            # Extract stitched unique global edge tracking numbers for this quad
            e_bottom = e_offset + get_global_edge_id(face, i, j, :bottom, N)
            e_right  = e_offset + get_global_edge_id(face, i, j, :right, N)
            e_top    = e_offset + get_global_edge_id(face, i, j, :top, N)
            e_left   = e_offset + get_global_edge_id(face, i, j, :left, N)
            
            # Place into the flat Cone array matching the DMPlex definition
            start_pos = cone_offsets[global_cell]
            cone_data[start_pos]   = e_bottom
            cone_data[start_pos+1] = e_right
            cone_data[start_pos+2] = e_top
            cone_data[start_pos+3] = e_left
            
            cell_idx += 1
        end
    end

    # 4. Generate Geometric Node Projections
    # For every unique vertex ID mapped out during stitching:
    for v_id in 1:num_vertices
        # 1. Grab raw un-normalized coordinates matching its cube position (-1.0 to 1.0)
        x_c, y_c, z_c = get_cube_coords(v_id, N)
        
        # 2. Project outward radially to the surface of the sphere
        mag = sqrt(x_c^2 + y_c^2 + z_c^2)
        coords[v_id] = [x_c / mag, y_c / mag, z_c / mag]
    end

    # 5. Automatically reverse the relations to get Upward Supports
    support_offsets, support_data = compute_support!(total_points, cone_offsets, cone_data)

    plex = JuliaDMPlex(total_points, cone_offsets, cone_data, support_offsets, support_data)
    return CubedSphereMesh(plex, coords)
end

