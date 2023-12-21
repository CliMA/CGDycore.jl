mutable struct Triangle
  centroid::Point
  area::Float64
end

mutable struct Polygon
  P::Array{Point,1}
  size::Int
  valid::Bool
end

function Polygon()
  P = Array{Point,1}(undef,MAXSIZE)
  size = 0
  valid = true
  return Polygon(
    P,
    size,
    valid,
  )
end  

function Polygon(Face::Face,Nodes::Array{Node,1})
  P = Array{Point}(undef,MAXSIZE)
  i = 0
  for N in Face.N
    i += 1
    P[i] = Nodes[N].P
  end
  return Polygon(
    P,
    i,
    true
    )
end

function area(triangles::Array{Triangle,1})
  area = 0.0
  for i = 1 : length(triangles)
    area += triangles[i].area
  end
  return area
end  

function distance2(P1::Point, P2::Point) 
  dx = P2.x - P1.x
  dy = P2.y - P1.y
  dz = P2.z - P1.z
  dx * dx + dy * dy + dz * dz
end

function distance(P1::Point, P2::Point) 
  dx = P2.x - P1.x
  dy = P2.y - P1.y
  dz = P2.z - P1.z
  sqrt(dx * dx + dy * dy + dz * dz)
end

function approx_eq(v1::Float64, v2::Float64) 
  abs(v1 - v2) <= EPS
end  

function approx_eq(P1::Point, P2::Point) 
  distance2(P1, P2) <= EPS2
end  

function lonlat2cart(P)
  x = sin(P[1])*sin(P[2])
  y = cos(P[1])*sin(P[2])
  z = cos(P[2])
  return Point(x, y, z)
end

function cart2lonlat(P)
  r = sqrt(P.x^2 + P.y^2 + P.z^2)
  lat = asin(P.z / r)
  lon = atan(P.y,P.x)
  if P.y < 0.0
    lon = lon + 2.0 * pi
  end
  return SVector{2}(lon, lat)
end

function lonlat2Polygon(Points::Array{Float64,2})

  lenP = size(Points,2)
  P = Arrayr{Point}(undef,MAXSIZE)
  @. P[1] = lonlat2cart(Points[:,1])
  isp = 2
  for i = 2 : lenP - 1 
    @. P[isp] = lonlat2cart(Points[:,i])
    if ~approx_eq(P[isp], P[isp - 1])
      isp += 1  
    end  
  end    
  @. P[isp] = lonlat2cart(Points[lenp])
  if ~(approx_eq(P[isp], P[1]) || approx_eq(P[isp], P[isp - 1]))
    isp += 1  
  end
  valid = isp > 2
  return Polygon(
    P,
    isp,
    valid,
  )
end  

function cart2Polygon(Points::Array{Float64,2})

  lenP = size(Points,2)
  P = Arrayr{Point}(undef,MAXSIZE)
  P[1] = Point(Points[:,1])
  isp = 2
  for i = 2 : lenP - 1 
    @. P[isp] = Point(Points[:,i])
    if ~approx_eq(P[isp], P[isp - 1])
      isp += 1  
    end  
  end    
  P[isp] = Points[lenp]
  if ~(approx_eq(P[isp], P[1]) || approx_eq(P[isp], P[isp - 1]))
    isp += 1  
  end
  valid = isp > 2
  return Polygon(
    P,
    isp,
    valid,
  )
end  

function triangulate(Polygon)
  triangles = Array{Triangle}(undef,0)
  P1 = Polygon.P[1]
  for i = 2 : Polygon.size - 1
    P2 = Polygon.P[i]  
    P3 = Polygon.P[i+1]  
    centroid = normalize(P1 + P2 + P3)
    A = AreaSphericalTriangle(P1,P2,P3)
    push!(triangles,Triangle(centroid,A))  
  end  
  return triangles
end
