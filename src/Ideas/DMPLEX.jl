using LinearAlgebra
struct JuliaDMPlex
    # Total range of points (Vertices, Edges, Faces, Cells lumped together)
    num_points::Int
    
    # Cone (Downward: Point -> sub-elements)
    cone_offsets::Vector{Int}
    cone_data::Vector{Int}
    
    # Support (Upward: Point -> parent elements)
    support_offsets::Vector{Int}
    support_data::Vector{Int}
end


function compute_support!(mesh_size::Int, cone_offsets::Vector{Int}, cone_data::Vector{Int})
    # Step 1: Count how many times each point appears as a sub-element (degree count)
    support_counts = zeros(Int, mesh_size)
    for p in 1:mesh_size
        start_idx = cone_offsets[p]
        end_idx = cone_offsets[p+1] - 1
        for sub_p in @view cone_data[start_idx:end_idx]
            support_counts[sub_p] += 1
        end
    end

    # Step 2: Build support_offsets using a prefix sum (cumulative sum)
    support_offsets = Vector{Int}(undef, mesh_size + 1)
    support_offsets[1] = 1
    for i in 1:mesh_size
        support_offsets[i+1] = support_offsets[i] + support_counts[i]
    end

    # Step 3: Populate support_data using a working cursor array
    support_data = Vector{Int}(undef, length(cone_data))
    current_positions = copy(support_offsets) # Tracks where to insert next item

    for p in 1:mesh_size
        start_idx = cone_offsets[p]
        end_idx = cone_offsets[p+1] - 1
        for sub_p in @view cone_data[start_idx:end_idx]
            # 'p' is a parent/support of 'sub_p'
            dest_idx = current_positions[sub_p]
            support_data[dest_idx] = p
            current_positions[sub_p] += 1
        end
    end

    return support_offsets, support_data
end

function create_two_triangles()
    num_points = 11 # 4 vertices + 5 edges + 2 cells

    # 1-indexed array of offsets (length = num_points + 1)
    cone_offsets = zeros(Int, num_points + 1)

    # Vertices (1-4) have no cones (base layer) -> 0 entries
    # Edges (5-9) have 2 vertices each -> 2 entries
    # Cells (10-11) have 3 edges each -> 3 entries
    counts = [0, 0, 0, 0,  2, 2, 2, 2, 2,  3, 3]

    cone_offsets[1] = 1
    for i in 1:num_points
        cone_offsets[i+1] = cone_offsets[i] + counts[i]
    end

    # Raw topology data
    cone_data = [
        1, 2,  # Edge 5
        2, 3,  # Edge 6
        3, 1,  # Edge 7 (Shared edge!)
        3, 4,  # Edge 8
        4, 1,  # Edge 9
        5, 6, 7,   # Cell 10 (Triangle 1)
        7, 8, 9    # Cell 11 (Triangle 2)
    ]

    # Compute the inverse relation automatically
    support_offsets, support_data = compute_support!(num_points, cone_offsets, cone_data)

    return JuliaDMPlex(num_points, cone_offsets, cone_data, support_offsets, support_data)
end


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


# Helper to query the graph
get_support(m::JuliaDMPlex, p::Int) = @view m.support_data[m.support_offsets[p]:m.support_offsets[p+1]-1]

# Execution
mesh = create_two_triangles()

# Let's inspect the shared edge (Edge 7)
println("Cells sharing Edge 7: ", get_support(mesh, 7))
# Output will be: [10, 11]

# Let's inspect Vertex 1
println("Edges meeting at Vertex 1: ", get_support(mesh, 1))
# Output will show the connected edge IDs

