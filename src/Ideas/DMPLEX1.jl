struct StandardDMPlex
    total_points::Int
    # Standard top-down PETSc ranges
    c_range::UnitRange{Int}
    f_range::UnitRange{Int}
    e_range::UnitRange{Int}
    v_range::UnitRange{Int}

    cone_offsets::Vector{Int}
    cone_data::Vector{Int}
    support_offsets::Vector{Int}
    support_data::Vector{Int}
end

"""
Constructs a standard PETSc-ordered DMPlex from a bottom-up user specification.
- user_edge_to_v: Vector of pairs/vectors mapping User Edge ID -> User Vertex IDs
- user_face_to_e: Vector mapping User Face ID -> User Edge IDs
- user_cell_to_f: Vector mapping User Cell ID -> User Face IDs
"""
function create_standard_plex_from_bottom_up(
    num_vertices::Int,
    user_edge_to_v::Vector{Vector{Int}},
    user_face_to_e::Vector{Vector{Int}},
    user_cell_to_f::Vector{Vector{Int}}
)
    num_edges = length(user_edge_to_v)
    num_faces = length(user_face_to_e)
    num_cells = length(user_cell_to_f)
    total_points = num_vertices + num_edges + num_faces + num_cells

    # 1. Define standard PETSc top-down ranges
    c_range = 1:num_cells
    f_range = (num_cells + 1):(num_cells + num_faces)
    e_range = (num_cells + num_faces + 1):(num_cells + num_faces + num_edges)
    v_range = (num_cells + num_faces + num_edges + 1):total_points

    # 2. Setup ID Translation Helpers (User ID -> Standard DMPlex ID)
    # User Vertices are 1:num_vertices
    transform_v(u_v) = v_range.start - 1 + u_v
    # User Edges are 1:num_edges
    transform_e(u_e) = e_range.start - 1 + u_e
    # User Faces are 1:num_faces
    transform_f(u_f) = f_range.start - 1 + u_f

    # 3. Calculate Cone Sizes for the Standard Top-Down Layout
    cone_counts = zeros(Int, total_points)
    
    for (i, c_id) in enumerate(c_range)
        cone_counts[c_id] = length(user_cell_to_f[i]) # Cells point to Faces
    end
    for (i, f_id) in enumerate(f_range)
        cone_counts[f_id] = length(user_face_to_e[i]) # Faces point to Edges
    end
    for (i, e_id) in enumerate(e_range)
        cone_counts[e_id] = length(user_edge_to_v[i]) # Edges point to Vertices
    end
    # Vertices have cone size 0 (bottom of the standard DAG)

    # 4. Generate Standard Cone Offsets
    cone_offsets = zeros(Int, total_points + 1)
    cone_offsets[1] = 1
    for i in 1:total_points
        cone_offsets[i+1] = cone_offsets[i] + cone_counts[i]
    end

    # 5. Populate Standard Cone Data with Translated IDs
    cone_data = Vector{Int}(undef, cone_offsets[end] - 1)

    # Fill Cells -> Faces
    for (i, c_id) in enumerate(c_range)
        idx = cone_offsets[c_id]
        for (j, u_f) in enumerate(user_cell_to_f[i])
            cone_data[idx + j - 1] = transform_f(u_f)
        end
    end

    # Fill Faces -> Edges
    for (i, f_id) in enumerate(f_range)
        idx = cone_offsets[f_id]
        for (j, u_e) in enumerate(user_face_to_e[i])
            cone_data[idx + j - 1] = transform_e(u_e)
        end
    end

    # Fill Edges -> Vertices
    for (i, e_id) in enumerate(e_range)
        idx = cone_offsets[e_id]
        for (j, u_v) in enumerate(user_edge_to_v[i])
            cone_data[idx + j - 1] = transform_v(u_v)
        end
    end

    # 6. Compute Standard Support Arrays (Upward Relations)
    support_offsets, support_data = compute_support!(total_points, cone_offsets, cone_data)

    return StandardDMPlex(
        total_points, c_range, f_range, e_range, v_range,
        cone_offsets, cone_data, support_offsets, support_data
    )
end

# Verification pipeline for Constructive-to-Standard DMPlex transformation

# Let's define a minimal test mesh:
# A single square face partitioned into a 2x2 grid of Quads.
# User input data is defined strictly Bottom-Up:
function run_dmplex_verification()
    println("=== STEP 1: Defining Constructive Input Data ===")
    num_vertices = 9  # 3x3 grid of vertices

    # User Edges defined by two user Vertex IDs (1-indexed)
    user_edge_to_v = [[1,2], [2, 3], [4, 5], [5, 6], [7, 8], [8, 9], [1, 4], [4, 7], [2, 5], [5, 8], [3, 6], [6, 9]]

    # User Faces defined by four user Edge IDs
    user_face_to_e = [[1, 9, 3, 7], [2, 10, 4, 
        [4, 12, 6, 10]  # Face 4 (Top-Right Quad)
    ]

    # User Cells defined by their constituent Face IDs
    # Since this is a 2D manifold shell, each Cell corresponds to exactly 1 Face
    user_cell_to_f = [[1], [2], [3], [4]]

    println("Input counts -> Vertices: $num_vertices, Edges: $(length(user_edge_to_v)), Faces: $(length(user_face_to_e)), Cells: $(length(user_cell_to_f))")

    println("\n=== STEP 2: Running Transformation Engine ===")
    plex = create_standard_plex_from_bottom_up(
        num_vertices, user_edge_to_v, user_face_to_e, user_cell_to_f
    )

    println("Transformation successful! Internal standard DMPlex ranges:")
    println("  Cells (Dim 2) : ", plex.c_range)
    println("  Faces (Dim 2) : ", plex.f_range)
    println("  Edges (Dim 1) : ", plex.e_range)
    println("  Vertices (Dim 0): ", plex.v_range)

    println("\n=== STEP 3: Verifying Top-Down Topological DAG ===")

    # Query Helper
    get_cone(p) = @view plex.cone_data[plex.cone_offsets[p]:plex.cone_offsets[p+1]-1]
    get_support(p) = @view plex.support_data[plex.support_offsets[p]:plex.support_offsets[p+1]-1]

    # Test 1: Pick Standard Cell 1. Does it point to a Face, which points to Edges, which point to Vertices?
    test_cell = 1
    faces_of_cell = get_cone(test_cell)
    println("Standard Cell $test_cell points down to Face(s): $faces_of_cell")

    target_face = faces_of_cell[1]
    edges_of_face = get_cone(target_face)
    println("  Standard Face $target_face points down to Edges: $edges_of_face")

    target_edge = edges_of_face[1]
    vertices_of_edge = get_cone(target_edge)
    println("    Standard Edge $target_edge points down to Vertices: $vertices_of_edge")

    # Test 2: Check standard Upward Support connectivity.
    # Let's find an internal edge. Standard edge 13 corresponds to user horizontal edge 5.
    # It should be shared by Face 3 and Face 4 internally.
    test_edge = 13
    println("\nUpward support verification:")
    println("  Standard Edge $test_edge is bounded upward by Faces: $(get_support(test_edge))")

    return plex
end

# To execute, ensure compute_support! and create_standard_plex_from_bottom_up are in scope:
plex_verified = run_dmplex_verification()

