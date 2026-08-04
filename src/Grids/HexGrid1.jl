using Printf

# Struct definitions representing the mesh topology
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
    node_ids::Vector{Int}  # 6 nodes forming the hexagon (counter-clockwise)
    edge_ids::Vector{Int}  # 6 edges forming the hexagon
    is_boundary::Bool
end

struct HexGrid
    nodes::Vector{Node}
    edges::Vector{Edge}
    cells::Vector{Cell}
end

"""
Generate a HexGrid with `num_rings` around a central cell.
`r` is the side length (radius) of each hexagon.
"""
function generate_hex_grid(num_rings::Int; r::Float64=1.0)
    # Axial directions for hex grid centers (q, r)
    axial_dirs = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
    
    # 1. Gather all cell centers in axial coordinates
    cell_centers = Tuple{Int, Int, Int}[] # (q, r, ring)
    push!(cell_centers, (0, 0, 0)) # Center cell
    
    for ring in 1:num_rings
        q, r_coord = -ring, ring
        for d in 1:6
            dq, dr = axial_dirs[d]
            for step in 1:ring
                push!(cell_centers, (q, r_coord, ring))
                q += dq
                r_coord += dr
            end
        end
    end

    # Helper maps for unique node and edge deduplication
    node_map = Dict{Tuple{Float64, Float64}, Int}()
    nodes = Node[]

    edge_map = Dict{Tuple{Int, Int}, Int}()
    edges = Edge[]

    cells = Cell[]

    # Helper: get or create node ID based on rounded coordinates
    function get_node_id(x::Float64, y::Float64)
        key = (round(x, digits=6), round(y, digits=6))
        if haskey(node_map, key)
            return node_map[key]
        else
            id = length(nodes) + 1
            node_map[key] = id
            # Initial boundary status is false; we will detect actual boundaries after
            push!(nodes, Node(id, key[1], key[2], false))
            return id
        end
    end

    # Helper: get or create edge ID given two node IDs
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

    # 2. Build Cells, Nodes, and Edges
    for (cell_idx, (q, r_coord, ring)) in enumerate(cell_centers)
        # Center in Cartesian space (Pointy-topped orientation)
        cx = r * sqrt(3) * (q + r_coord / 2.0)
        cy = r * (3/2 * r_coord)

        c_node_ids = Int[]
        c_edge_ids = Int[]

        # 6 vertices of the hexagon
        for i in 0:5
            angle = π/6 + i * (π/3)
            vx = cx + r * cos(angle)
            vy = cy + r * sin(angle)
            push!(c_node_ids, get_node_id(vx, vy))
        end

        # 6 edges connecting adjacent vertices
        for i in 1:6
            n1 = c_node_ids[i]
            n2 = c_node_ids[mod1(i + 1, 6)]
            push!(c_edge_ids, get_edge_id(n1, n2))
        end

        # Cell sits on grid boundary if it belongs to the outer ring
        is_cell_boundary = (ring == num_rings)
        push!(cells, Cell(cell_idx, ring, c_node_ids, c_edge_ids, is_cell_boundary))
    end

    # 3. Mark Boundary Edges & Nodes
    # Count how many cells share each edge
    edge_cell_count = zeros(Int, length(edges))
    for c in cells
        for e_id in c.edge_ids
            edge_cell_count[e_id] += 1
        end
    end

    # An edge is a boundary edge if it is only shared by 1 cell
    boundary_nodes_set = Set{Int}()
    updated_edges = Edge[]
    for e in edges
        is_bnd = (edge_cell_count[e.id] == 1)
        push!(updated_edges, Edge(e.id, e.node1, e.node2, is_bnd))
        if is_bnd
            push!(boundary_nodes_set, e.node1)
            push!(boundary_nodes_set, e.node2)
        end
    end

    # Update node boundary statuses
    updated_nodes = [
        Node(n.id, n.x, n.y, n.id in boundary_nodes_set)
        for n in nodes
    ]

    return HexGrid(updated_nodes, updated_edges, cells)
end

# Printable inspection function
function print_hex_grid(grid::HexGrid)
    println("==================================================")
    println("                 CELLS SUMMARY                   ")
    println("==================================================")
    @printf("%-6s | %-6s | %-12s | %-24s\n", "ID", "Ring", "Boundary?", "Node IDs (1..6)")
    println("-"^58)
    for c in grid.cells
        @printf("%-6d | %-6d | %-12s | %s\n", c.id, c.ring, c.is_boundary ? "Yes" : "No", string(c.node_ids))
    end

    println("\n==================================================")
    println("                 EDGES SUMMARY                   ")
    println("==================================================")
    @printf("%-6s | %-14s | %-12s\n", "ID", "Node Pair", "Boundary?")
    println("-"^40)
    for e in grid.edges
        @printf("%-6d | (%-3d -> %-3d)   | %-12s\n", e.id, e.node1, e.node2, e.is_boundary ? "Yes" : "No")
    end

    println("\n==================================================")
    println("                 NODES SUMMARY                   ")
    println("==================================================")
    @printf("%-6s | %-20s | %-12s\n", "ID", "Coordinates (x, y)", "Boundary?")
    println("-"^46)
    for n in grid.nodes
        coord_str = @sprintf("(%.2f, %.2f)", n.x, n.y)
        @printf("%-6d | %-20s | %-12s\n", n.id, coord_str, n.is_boundary ? "Yes" : "No")
    end
end

# Example Execution: Grid with 1 Ring (7 total hexagons)
grid = generate_hex_grid(5)
print_hex_grid(grid)
