using Printf
using NCDatasets
using WriteVTK

struct Node
    id::Int
    x::Float64
    y::Float64
    is_boundary::Bool
end

struct Edge
    id::Int
    node1::Int
    node2::Int
    is_boundary::Bool
end

struct Cell
    id::Int
    ring::Int
    node_ids::Vector{Int}
    edge_ids::Vector{Int}
    is_boundary::Bool
end

struct HexGrid
    nodes::Vector{Node}
    edges::Vector{Edge}
    cells::Vector{Cell}
end

"""
Generate a HexGrid with exact concentric rings around (0,0).
"""
function generate_hex_grid(num_rings::Int; r::Float64=1.0)
    axial_dirs = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]

    # 1. Generate Cell Centers in axial coordinates (q, r)
    cell_centers = Tuple{Int, Int, Int}[]
    push!(cell_centers, (0, 0, 0)) # Ring 0

    for ring in 1:num_rings
        curr_q = axial_dirs[5][1] * ring
        curr_r = axial_dirs[5][2] * ring
        for d in 1:6
            dq, dr = axial_dirs[d]
            for _ in 1:ring
                push!(cell_centers, (curr_q, curr_r, ring))
                curr_q += dq
                curr_r += dr
            end
        end
    end

    # Exact vertex integer offsets relative to hexagon center (pointy-topped)
    # Scaled by 1/sqrt(3) to form standard radius 1
    vertex_offsets = [
        ( 1/sqrt(3),  0.0),
        ( 0.5/sqrt(3), 0.5),
        (-0.5/sqrt(3), 0.5),
        (-1/sqrt(3),  0.0),
        (-0.5/sqrt(3), -0.5),
        ( 0.5/sqrt(3), -0.5)
    ]

    node_map = Dict{Tuple{Int64, Int64}, Int}()
    nodes = Node[]
    edge_map = Dict{Tuple{Int, Int}, Int}()
    edges = Edge[]
    cells = Cell[]

    # Helper: Deduplicate nodes using EXACT integer grid coordinates
    # Grid coordinates scaled by 1e6 to avoid float inaccuracies
    function get_node_id(x::Float64, y::Float64)
        ix = round(Int64, x * 1e8)
        iy = round(Int64, y * 1e8)
        key = (ix, iy)
        if haskey(node_map, key)
            return node_map[key]
        else
            id = length(nodes) + 1
            node_map[key] = id
            push!(nodes, Node(id, x, y, false))
            return id
        end
    end

    function get_edge_id(n1::Int, n2::Int)
        key = n1 < n2 ? (n1, n2) : (n2, n1)
        if haskey(edge_map, key)
            return edge_map[key]
        else
            id = length(edges) + 1
            edge_map[key] = id
            push!(edges, Edge(id, key[1], key[2], false))
            return id
        end
    end

    # 2. Build Cells & Vertices
    for (cell_idx, (q, r_coord, ring)) in enumerate(cell_centers)
        cx = r * sqrt(3) * (q + r_coord / 2.0)
        cy = r * (1.5 * r_coord)

        raw_node_ids = Int[]
        for i in 0:5
            angle = π/6 + i * (π/3)
            vx = cx + r * cos(angle)
            vy = cy + r * sin(angle)
            push!(raw_node_ids, get_node_id(vx, vy))
        end

        # Counter-clockwise sort relative to center
        c_node_ids = sort(raw_node_ids, by = nid -> begin
            n = nodes[nid]
            atan(n.y - cy, n.x - cx)
        end)

        c_edge_ids = Int[]
        for i in 1:6
            n1 = c_node_ids[i]
            n2 = c_node_ids[mod1(i + 1, 6)]
            push!(c_edge_ids, get_edge_id(n1, n2))
        end

        is_cell_boundary = (ring == num_rings)
        push!(cells, Cell(cell_idx, ring, c_node_ids, c_edge_ids, is_cell_boundary))
    end

    # 3. Boundary Flag Updates
    edge_cell_count = zeros(Int, length(edges))
    for c in cells
        for e_id in c.edge_ids
            edge_cell_count[e_id] += 1
        end
    end

    updated_edges = Edge[]
    for e in edges
        is_bnd = (edge_cell_count[e.id] == 1)
        push!(updated_edges, Edge(e.id, e.node1, e.node2, is_bnd))
    end

    # Count how many cells share each node
    node_cell_count = zeros(Int, length(nodes))
    for c in cells
        for n_id in c.node_ids
            node_cell_count[n_id] += 1
        end
    end

    # Nodes shared by 3 cells are fully interior (is_boundary = false)
    updated_nodes = [
        Node(n.id, n.x, n.y, node_cell_count[n.id] < 3)
        for n in nodes
    ]

    return HexGrid(updated_nodes, updated_edges, cells)
end

function export_to_netcdf(grid::HexGrid, filename::String="hex_mesh.nc")
    NCDataset(filename, "c") do ds
        # 1. Define Dimensions
        defDim(ds, "num_nodes", length(grid.nodes))
        defDim(ds, "num_edges", length(grid.edges))
        defDim(ds, "num_cells", length(grid.cells))
        defDim(ds, "nodes_per_cell", 6)

        # 2. Define Node Variables
        x_var = defVar(ds, "node_x", Float64, ("num_nodes",))
        y_var = defVar(ds, "node_y", Float64, ("num_nodes",))
        nbnd_var = defVar(ds, "node_is_boundary", Int32, ("num_nodes",))

        x_var[:] = [n.x for n in grid.nodes]
        y_var[:] = [n.y for n in grid.nodes]
        nbnd_var[:] = [n.is_boundary ? 1 : 0 for n in grid.nodes]

        # 3. Define Edge Variables
        e1_var = defVar(ds, "edge_node1", Int32, ("num_edges",))
        e2_var = defVar(ds, "edge_node2", Int32, ("num_edges",))
        ebnd_var = defVar(ds, "edge_is_boundary", Int32, ("num_edges",))

        e1_var[:] = [e.node1 for e in grid.edges]
        e2_var[:] = [e.node2 for e in grid.edges]
        ebnd_var[:] = [e.is_boundary ? 1 : 0 for e in grid.edges]

        # 4. Define Cell Variables
        conn_var = defVar(ds, "cell_node_ids", Int32, ("nodes_per_cell", "num_cells"))
        ring_var = defVar(ds, "cell_ring", Int32, ("num_cells",))
        cbnd_var = defVar(ds, "cell_is_boundary", Int32, ("num_cells",))

        # Fill connectivity matrix (6 x num_cells)
        conn_mat = hcat([c.node_ids for c in grid.cells]...)
        conn_var[:, :] = conn_mat
        ring_var[:] = [c.ring for c in grid.cells]
        cbnd_var[:] = [c.is_boundary ? 1 : 0 for c in grid.cells]
        
        # Attributes
        ds.attrib["title"] = "Concentric Hexagonal Mesh Export"
    end
    println("Saved NetCDF file to: $filename")
end

"""
Export grid to ParaView VTK (.vtu) file
"""
function export_to_paraview(grid::HexGrid, filename_prefix::String="hex_mesh")
    # 1. Prepare node coordinates array (3 x num_nodes for VTK)
    num_nodes = length(grid.nodes)
    points = zeros(Float64, 3, num_nodes)
    for (i, n) in enumerate(grid.nodes)
        points[1, i] = n.x
        points[2, i] = n.y
        points[3, i] = 0.0  # 2D plane
    end

    # 2. Build VTK Polygon Cells (6-vertex hexagons)
    cells = [MeshCell(VTKCellTypes.VTK_POLYGON, c.node_ids) for c in grid.cells]

    # 3. Create VTK Unstructured Grid and write point/cell scalar data
    vtk_grid(filename_prefix, points, cells) do vtk
        # Point/Node attributes
        vtk["node_is_boundary"] = [n.is_boundary ? 1 : 0 for n in grid.nodes]

        # Cell attributes
        vtk["cell_ring"] = [c.ring for c in grid.cells]
        vtk["cell_is_boundary"] = [c.is_boundary ? 1 : 0 for c in grid.cells]
    end
    println("Saved ParaView VTK file to: $(filename_prefix).vtu")
end

# --- Run Exports ---
grid = generate_hex_grid(10) # 2 rings = 19 hexagons
export_to_netcdf(grid, "hex_mesh.nc")
export_to_paraview(grid, "hex_mesh")
