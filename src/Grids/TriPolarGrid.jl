"""
    generate_tripolar_components(Ni, Nj)
    Builds the grid by explicitly defining edges first, then assembling quads from them.
"""
function generate_tripolar_components(Ni, Nj, R)
  # 1. Generate Vertices (Cartesian)
  FT = eltype(R)
  vertices = []
  for j in 1:Nj, i in 1:Ni
    lat = -90.0 + (j - 1) * (180.0 / (Nj - 1))
    lon = (i - 1) * (360.0 / Ni)
    ϕ, λ = deg2rad(lat), deg2rad(lon)
    push!(vertices, (x=R*cos(ϕ)*cos(λ), y=R*cos(ϕ)*sin(λ), z=R*sin(ϕ)))
  end
#   grid = Oceananigans.Grids.TripolarGrid(arch = Oceananigans.CPU(), FT;
  grid = Oceananigans.Grids.TripolarGrid(;
                      size = (Ni, Nj, 1),
                      southernmost_latitude = -80,
                      halo = (1, 1, 1),
                      radius = R,
                      z = (0, 1),
                      north_poles_latitude = 55,
                      first_pole_longitude = 70, # second pole is at longitude `first_pole_longitude + 180ᵒ`
                      fold_topology = Oceananigans.Grids.RightFaceFolded)
  vertices = []
  for j in 1:Nj, i in 1:Ni
    lon = grid.λᶠᶠᵃ[i, j]
    lat = grid.φᶠᶠᵃ[i, j]
    ϕ, λ = deg2rad(lat), deg2rad(lon)
    push!(vertices, (x=R*cos(ϕ)*cos(λ), y=R*cos(ϕ)*sin(λ), z=R*sin(ϕ)))
  end

  get_v_idx(i, j) = (j - 1) * Ni + ((i - 1) % Ni + 1)

  # 2. Construct Edges
  # We store them in a Dict to easily retrieve them by their logical "direction"
  h_edges = Dict() # Horizontal (East-West)
  v_edges = Dict() # Vertical (North-South)
  all_edges = []
  edge_counter = 1

  for j in 1:Nj
    for i in 1:Ni
      v1 = get_v_idx(i, j)
      # Horizontal Edge
      v2 = get_v_idx(i + 1, j)
      h_edges[(i, j)] = edge_counter
      push!(all_edges, (v1, v2))
      edge_counter += 1

      # Vertical Edge (skip the very top row for standard quads)
      if j < Nj
        v3 = get_v_idx(i, j + 1)
        v_edges[(i, j)] = edge_counter
        push!(all_edges, (v1, v3))
        edge_counter += 1
      end
    end
  end

  # 3. Construct Quads FROM the Edge IDs
  quads_from_edges = []
  for j in 1:(Nj - 1)
    for i in 1:Ni
      # A quad at (i,j) is composed of:
      # Bottom edge: h_edges(i, j)
      # Right edge:  v_edges(i+1, j)  -- with periodic wrapping for i
      # Top edge:    h_edges(i, j+1)
      # Left edge:   v_edges(i, j)
          
      i_next = (i % Ni) + 1
        
      e_bottom = h_edges[(i, j)]
      e_right  = v_edges[(i_next, j)]
      e_top    = h_edges[(i, j + 1)]
      e_left   = v_edges[(i, j)]

      push!(quads_from_edges, (e_bottom, e_right, e_top, e_left))
    end
  end

  return vertices, all_edges, quads_from_edges
end

function TriPolarGrid(backend,FT,Ni,Nj,Rad,nz;order=true)

  nBar=[ 0  1   0   1
        -1  0  -1   0];
  Dim=3
  Type=Quad()
  Rad=Rad
  Form=SphericalGrid()

  verts, edgs, quads = generate_tripolar_components(Ni, Nj, Rad)

  NumNodes = length(verts)
  Nodes = map(1:NumNodes) do i
    Node()
  end
  NodeNumber = 1
  @inbounds for p = 1 : length(verts)
    x = verts[p].x
    y = verts[p].y
    z = verts[p].z
    Nodes[NodeNumber]=Node(Point(x,y,z),NodeNumber,' ');
    NodeNumber=NodeNumber+1;
  end

  NumEdges = length(edgs)
  Edges = map(1:NumEdges) do i
    Edge()
  end
  EdgeNumber = 1
  @inbounds for p = 1 : length(edgs)
    n1 = edgs[p][1]
    n2 = edgs[p][2]
    e = [n1,n2]
    Edges[EdgeNumber]=Edge(e,Nodes,EdgeNumber,EdgeNumber,"",EdgeNumber,Form=Form,Rad=Rad)
    EdgeNumber += 1
  end

  NumFaces = length(quads)
  Faces = map(1:NumFaces) do i
    Face()
  end
  @inbounds for p = 1 : length(quads)
    e1 = quads[p][1]
    e2 = quads[p][2]
    e3 = quads[p][3]
    e4 = quads[p][4]
    f = [e1,e2,e3,e4]
    (Faces[p], Edges) = Face(f,Nodes,Edges,p,"Sphere",OrientFaceSphere;
      Form=Form,Rad=Rad,P=zeros(Float64,0,0))
  end

  if order
    Orientation!(Edges,Faces);
    Renumbering!(Edges,Faces);
  end
  FacesInNodes!(Nodes,Faces)
  SortFacesInNodes!(Nodes,Faces)

  zP=zeros(nz)
  z=KernelAbstractions.zeros(backend,FT,nz+1)
  dzeta=zeros(nz)
  H=0.0
  nBar3 = zeros(0,0)
  NumNodesB = 0
  NumNodesG = 0
  NumEdgesB = 0
  NumEdgesG = 0
  NumFacesB = 0
  NumFacesG = 0
  AdaptGrid = ""
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

