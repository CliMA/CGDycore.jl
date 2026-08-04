struct NodeHex
    id::Int
    x::Float64
    y::Float64
    is_boundary::Bool
end

struct EdgeHex
    id::Int
    node1::Int
    node2::Int
    is_boundary::Bool
end

struct CellHex
    id::Int
    ring::Int
    node_ids::Vector{Int}
    edge_ids::Vector{Int}
    is_boundary::Bool
end


"""
Generate a HexGrid with exact concentric rings around (0,0).
"""
function generate_hex_grid(backend, FT, num_rings::Int, OrientFace, nz; r::Float64=50.0)
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
    nodes = NodeHex[]
    edge_map = Dict{Tuple{Int, Int}, Int}()
    edges = EdgeHex[]
    cells = CellHex[]

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
            push!(nodes, NodeHex(id, x, y, false))
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
            push!(edges, EdgeHex(id, key[1], key[2], false))
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
        push!(cells, CellHex(cell_idx, ring, c_node_ids, c_edge_ids, is_cell_boundary))
    end

    # 3. Boundary Flag Updates
    edge_cell_count = zeros(Int, length(edges))
    for c in cells
        for e_id in c.edge_ids
            edge_cell_count[e_id] += 1
        end
    end

    updated_edges = EdgeHex[]
    for e in edges
        is_bnd = (edge_cell_count[e.id] == 1)
        push!(updated_edges, EdgeHex(e.id, e.node1, e.node2, is_bnd))
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
        NodeHex(n.id, n.x, n.y, node_cell_count[n.id] < 3)
        for n in nodes
    ]

    nBar=[ 0  1   0   1
        -1  0  -1   0]
    nBar3=[ 0  1   0   1
             -1  0  -1   0
             0   0  0    0]
    Dim = 3
    Type = Quad()
    Form = CartesianGrid()

    NumNodes = length(updated_nodes)
    Nodes = map(1:NumNodes) do i
    Node()
    end

    NodeNumber = 1
    for iN = 1 : NumNodes
        if updated_nodes[iN].is_boundary
            TypeN = 'B'
        else
            TypeN = ' '
        end    
        x = updated_nodes[iN].x
        y = updated_nodes[iN].y
        Nodes[NodeNumber]=Node(Point([x,y,0.0]),NodeNumber,TypeN)
        NodeNumber += 1
    end    

    NumEdges = length(updated_edges)
    Edges = map(1:NumEdges) do i
    Edge()
    end

    EdgeNumber = 1
    for iE = 1 : NumEdges
        N1 = updated_edges[iE].node1 
        N2 = updated_edges[iE].node2 
        if Nodes[N1].Type == 'B' && Nodes[N2].Type == 'B'
            TypeE = "B"
        else
            TypeE = ""
        end    
        Edges[EdgeNumber]=Edge([N1,N2],Nodes,EdgeNumber,EdgeNumber,TypeE,0)
        EdgeNumber += 1
    end    

    NumFaces = length(cells)
    Faces=map(1:NumFaces) do i
       Face()
    end

    FaceNumber = 1
    TypeF = "o"
    for iF = 1 : NumFaces
        E1 = cells[iF].edge_ids[1]
        E2 = cells[iF].edge_ids[2]
        E3 = cells[iF].edge_ids[3]
        E4 = cells[iF].edge_ids[4]
        E5 = cells[iF].edge_ids[5]
        E6 = cells[iF].edge_ids[6]
        (Faces[FaceNumber],Edges)=Face([E1,E2,E3,E4,E5,E6],Nodes,Edges,FaceNumber,TypeF,OrientFace;P=zeros(Float64,0,0))
        FaceNumber += 1 
    end    
        
  FacesInNodes!(Nodes,Faces)

  zP=zeros(nz)
  z=KernelAbstractions.zeros(backend,FT,nz+1)
  dzeta=zeros(nz)
  H=0.0
  NumFacesB = 0
  NumFacesG = 0
  NumEdgesB = 0
  NumEdgesG = 0
  NumNodesB = 0
  NumNodesG = 0
  nBar3 = zeros(0,0)
  AdaptGrid = ""
  Rad = 1.0
  EF=KernelAbstractions.zeros(backend,Int,0,0)
  FE=KernelAbstractions.zeros(backend,Int,0,0)
  return GridStruct{FT,
                    typeof(EF),
                    typeof(z)}(
    nz,
    zP,
    z,
    dzeta,
    H,
    NumFaces,
    NumFacesB,
    NumFacesG,
    Faces,
    NumEdges,
    NumEdgesB,
    NumEdgesG,
    Edges,
    NumNodes,
    NumNodesB,
    NumNodesG,
    Nodes,
    Form,
    Type,
    Dim,
    Rad,
    nBar3,
    nBar,
    AdaptGrid,
    EF,
    FE,
    )

end

