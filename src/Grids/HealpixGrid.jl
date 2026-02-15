function HealpixGrid(backend,FT,ns,Rad,nz;order=true)
  """
  healpix.jl: Generator of HEALPix grid

  Healpix points are plotted on the sphere or in the lon-lat plane.
  Rhomboids are created with their centroids in healpix points.
  On the sphere, the edges of the rhomboids are great arcs.
  """
  nBar=[ 0  1   0   1
        -1  0  -1   0];
  Dim=3
  Type=Quad()
  Rad=Rad
  Form=SphericalGrid()


  npix = 12 * ns * ns
  nrings = 4 * ns - 1
  rad2deg = 180 / pi
  deg2rad = 1 / rad2deg
  ncap_ids = zeros(Int, 2)

  # HEALPix dual points
  dualp = zeros(npix, 2)

  # North cap
  for p in 1:npix
    ph = 0.5 * p
    i = floor(Int,sqrt(ph - sqrt(floor(ph)))) + 1
    if i > ns
      ncap_ids[1] = p - 1  # Adjust for Julia's 1-based index
      break
    end
    j = p - 2 * i * (i - 1)
    dualp[p, 1] = acos(1 - i^2 / (3 * ns^2)) * rad2deg
    dualp[p, 2] = (0.25 * pi * (2 * j - 1) / i) * rad2deg
  end

  # North belt
  for p in ncap_ids[1] + 1:npix
    ph = p - 1 - 2 * ns * (ns - 1)
    i = floor(Int,0.25 * ph / ns) + ns
    if i > 2 * ns
      ncap_ids[2] = p - 1  # Adjust for Julia's 1-based index
      break
    end
    j = 1 + ph % (4 * ns)
    s = (i - ns + 1) % 2
    dualp[p, 1] = acos((4 * ns - 2 * i) / (3 * ns)) * rad2deg
    dualp[p, 2] = (0.25 * pi * (2 * j - s) / ns) * rad2deg
  end

  # South belt & south cap
  for p in ncap_ids[2] + 1:npix
    dualp[p, 1] = 180 - dualp[npix - p + 1, 1]
    dualp[p, 2] = dualp[npix - p + 1, 2]
  end

  # HEALPix vertex points
  vertp = zeros(npix + 2, 2)
  vertpC = zeros(npix + 2, 3)

  vertp[1, 1] = 0
  vertp[1, 2] = 180
  vertp[npix + 2, 1] = 180
  vertp[npix + 2, 2] = 180

  # North cap
  for p in 1:ncap_ids[1]
    ph = 0.5 * p
    i = floor(Int,sqrt(ph - sqrt(floor(ph)))) + 1
    last_id = 2 * i * (i - 1)
    j = p - last_id
    vertp[p + 1, 1] = dualp[p, 1]
    pp = floor(Int,j > 1 ? p - 2 : last_id + 4 * i - 1) + 1
    vertp[p + 1, 2] = dualp[pp, 2] + dualp[p, 2]
    vertp[p + 1, 2] = p > pp ? 0.5 * vertp[p + 1, 2] : 0.5 * (vertp[p + 1, 2] + 360)

  end

  # North belt
  for p in ncap_ids[1] + 1:ncap_ids[2]
    ph = p - 2 * ns * (ns - 1) - 1
    i = floor(Int,0.25 * ph / ns) + ns
    j = 1 + ph % (4 * ns)
    vertp[p + 1, 1] = dualp[p, 1]
    pp = j > 1 ? p - 2 : ncap_ids[1] + 4 * ns * (i - ns) - 1
    pp += 1
    vertp[p + 1, 2] = dualp[pp, 2] + dualp[p, 2]
    vertp[p + 1, 2] = p > pp ? 0.5 * vertp[p + 1, 2] : 0.5 * (vertp[p + 1, 2] + 360)
  end

  # South belt & south cap
  for p in ncap_ids[2] + 1:npix
    vertp[p + 1, 1] = 180 - vertp[npix - p + 2, 1]
    vertp[p + 1, 2] = vertp[npix - p + 2, 2]
  end


  # HEALPix cells
  el_vid = zeros(Int, npix, 4)

  # North cap for cells
  for el in 1:ncap_ids[1]
    ph = 0.5 * el
    i = floor(sqrt(ph - sqrt(floor(ph)))) + 1
    j = el - 2 * i * (i - 1)
    if i > 1
      if j < 4 * i
        el_vid[el,1] = j-floor(Int,((j-1)/i))
      else
        el_vid[el,1] = 1
      end
    else
      el_vid[el,1] = 0
    end
    if i > 2
      el_vid[el,1] += 2*(i-1)*(i-2)
    end
    el_vid[el, 2] = j + 2 * i * (i - 1)
    el_vid[el, 4] = j < 4 * i ? j + 1 + 2 * i * (i - 1) : 1 + 2 * i * (i - 1)
#   el_vid[el, 3] = el < 2 * (ns - 1) * ns ? 2 * (i + 1) * i + j + 1 : 2 * (i + 1) * i + j
    if el - 1 < 2*(ns-1)*ns 
        el_vid[el,3] = 2*(i+1)*i + j + 1 + floor(Int,(j-1)/i)
    else
        el_vid[el,3] = 2*(i+1)*i + j
    end    
  end

  # North belt for cells
  for el in ncap_ids[1] + 1:6 * ns * ns + 2 * ns
    ph = el - 1 - 2 * ns * (ns - 1)
    i = floor(0.25 * ph / ns) + ns
    j = 1 + ph % (4 * ns)
    s = (i - ns + 1) % 2
    el_vid[el, 1] = 2 * ns * (ns - 1) + 4 * ns * (i - ns - 1) + 1
    if j == 4 * ns && s == 0
      el_vid[el, 1] -= s  
    else
      el_vid[el, 1] += j - s  
    end  
    el_vid[el, 2] = el - 1 + 1
    el_vid[el, 4] = el - 1 + 2 - (j == 4 * ns ? 4 * ns : 0)
    if el - 1 < 6*ns*ns-2*ns 
      el_vid[el,3] = 1 - s 
      if j == 4 * ns && s == 0
        el_vid[el,3] += 2*ns*(ns-1)+4*ns*(i-ns+1) 
      else
         el_vid[el,3] += el - 1+4*ns+1 
      end   
    end  
  end

# Equator for cells
  for el in 6 * ns * ns - 2 * ns + 1:6 * ns * ns + 2 * ns
    el_vid[el, 3] = npix + 1 - el_vid[el, 1]
  end

# South belt & south cap for cells
  for el in 6 * ns * ns + 2 * ns + 1:npix
    el_vid[el, 1] = npix + 1 - el_vid[npix - el + 1, 1]
    el_vid[el, 4] = npix + 1 - el_vid[npix - el + 1, 2]
    el_vid[el, 2] = npix + 1 - el_vid[npix - el + 1, 4]
#   el_vid[el, 3] = el_vid[npix - el + 1, 3] > 6 * ns * ns - 2 * ns && el_vid[npix - el + 1, 3] <= 
#     6 * ns * ns + 2 * ns ? el_vid[npix - el + 1, 3] : npix + 1 - el_vid[npix - el + 1, 3]
    if el_vid[npix-el+1,3] > 6*ns*ns-2*ns && el_vid[npix-el+1,3] <= 6*ns*ns+2*ns 
        el_vid[el,3] = el_vid[npix-el+1,3]
    else
        el_vid[el,3] = npix + 1 - el_vid[npix-el+1,3]  
    end    
  end

  @. el_vid += 1

  Node2Edge = Dict()
  E1 = zeros(Int,0)
  E2 = zeros(Int,0)
  iE = 1
  for el = 1 : size(el_vid,1)
    # Edge 1
    n1 = el_vid[el,1]
    n2 = el_vid[el,2]
    if n1 < n2
      Node2Edge[[n1,n2]] = iE
      iE += 1
      push!(E1,n1)
      push!(E2,n2)
    end  
    # Edge 2
    n1 = el_vid[el,2]
    n2 = el_vid[el,3]
    if n1 < n2
      Node2Edge[[n1,n2]] = iE
      iE += 1
    end  
    # Edge 3
    n1 = el_vid[el,3]
    n2 = el_vid[el,4]
    if n1 < n2
      Node2Edge[[n1,n2]] = iE
      iE += 1
    end  
    # Edge 4
    n1 = el_vid[el,4]
    n2 = el_vid[el,1]
    if n1 < n2
      Node2Edge[[n1,n2]] = iE
      iE += 1
    end  
  end   
  NumEdges = iE - 1

  F = zeros(Int,size(el_vid))
  for el = 1 : size(el_vid,1)
    # Edge 1
    n1 = el_vid[el,1]
    n2 = el_vid[el,2]
    if n1 < n2
      e = [n1,n2]  
    else
      e = [n2,n1]  
    end  
    F[el,1] = Node2Edge[e]
    # Edge 2
    n1 = el_vid[el,2]
    n2 = el_vid[el,3]
    if n1 < n2
      e = [n1,n2]  
    else
      e = [n2,n1]  
    end  
    F[el,2] = Node2Edge[e]
    # Edge 3
    n1 = el_vid[el,3]
    n2 = el_vid[el,4]
    if n1 < n2
      e = [n1,n2]  
    else
      e = [n2,n1]  
    end  
    F[el,3] = Node2Edge[e]
    # Edge 4
    n1 = el_vid[el,4]
    n2 = el_vid[el,1]
    if n1 < n2
      e = [n1,n2]  
    else
      e = [n2,n1]  
    end  
    F[el,4] = Node2Edge[e]
  end

  NumNodes = size(vertp,1)
  Nodes = map(1:NumNodes) do i
    Node()
  end
  NodeNumber = 1
  @inbounds for p = 1 : size(vertp,1)
    lat = vertp[p,1] - 90.0
    lon = vertp[p,2] 
    x = sind(lon)*cosd(lat)
    y = cosd(lon)*cosd(lat)
    z = sind(lat)
    Nodes[NodeNumber]=Node(Point(x,y,z),NodeNumber,' ');
    NodeNumber=NodeNumber+1;
  end

  Edges = map(1:NumEdges) do i
    Edge()
  end
  EdgeNumber = 1

  for el = 1 : size(el_vid,1)
    # Edge 1
    n1 = el_vid[el,1]
    n2 = el_vid[el,2]
    if n1 < n2
      e = [n1,n2]  
      Edges[EdgeNumber]=Edge(e,Nodes,EdgeNumber,EdgeNumber,"",EdgeNumber,Form=Form,Rad=Rad)
      EdgeNumber += 1
    end  
    # Edge 2
    n1 = el_vid[el,2]
    n2 = el_vid[el,3]
    if n1 < n2
      e = [n1,n2]  
      Edges[EdgeNumber]=Edge(e,Nodes,EdgeNumber,EdgeNumber,"",EdgeNumber,Form=Form,Rad=Rad)
      EdgeNumber += 1
    end  
    # Edge 3
    n1 = el_vid[el,3]
    n2 = el_vid[el,4]
    if n1 < n2
      e = [n1,n2]  
      Edges[EdgeNumber]=Edge(e,Nodes,EdgeNumber,EdgeNumber,"",EdgeNumber,Form=Form,Rad=Rad)
      EdgeNumber += 1
    end  
    # Edge 4
    n1 = el_vid[el,4]
    n2 = el_vid[el,1]
    if n1 < n2
      e = [n1,n2]  
      Edges[EdgeNumber]=Edge(e,Nodes,EdgeNumber,EdgeNumber,"",EdgeNumber,Form=Form,Rad=Rad)
      EdgeNumber += 1
    end  
  end

  NumFaces = size(F,1)
  Faces = map(1:NumFaces) do i
    Face()
  end

    
  for i = 1 : NumFaces
    (Faces[i], Edges) = Face(F[i,:],Nodes,Edges,i,"Sphere",OrientFaceSphere;
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
