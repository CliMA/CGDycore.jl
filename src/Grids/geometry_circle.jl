mutable struct GreatCircle
  P1::Point
  P2::Point
end

function sphere2cart(lam,phi,r)
x=cos(lam)*cos(phi)*r
y=sin(lam)*cos(phi)*r
z=sin(phi)*r
return [x;y;z];
end

function sphereDeg2cart(lam,phi,r)
x=cosd(lam)*cosd(phi)*r
y=sind(lam)*cosd(phi)*r
z=sind(phi)*r
return [x;y;z];
end


@inline function cart2sphere(x,y,z)
  FT = eltype(x)
  r = sqrt(x^2 + y^2 + z^2)
  phi = asin(z / r)
  lam = atan(y,x) 
  if y < FT(0)
    lam = lam + FT(2*pi)
  end  

# ϕ = atan(z, hypot(y, x))
# if abs(ϕ) == 90
#     λ = zero(ϕ)
# else
#     λ = atan(y, x)
# end

# lam = FT(0)
# if abs(abs(phi) - FT(pi/2))>FT(1.e-14)
#   lam = atan(y,x) # TODO: check translation with Oswald
#   if lam < FT(0.0)
#     lam = lam + FT(2*pi)
#   end
# end
  return lam,phi,r
end

function cart2sphereDeg(x,y,z)
r=sqrt(x^2+y^2+z^2)
phi=asind(z/r)

# ϕ = atan(z, hypot(y, x))
# if abs(ϕ) == 90
#     λ = zero(ϕ)
# else
#     λ = atan(y, x)
# end

lam=-180
if abs(abs(phi)-90.0)>1.e-14
  lam=atand(y,x) # TODO: check translation with Oswald
end
return (float(lam),float(phi),r)
end

function SizeGreatCircle(C::GreatCircle) 
  return acos(dot(C.P1,C.P2)/(norm(C.P1)*norm(C.P2)))
end  

function SizeGreatCircle(Lon1,Lat1,Lon2,Lat2) 
  return acos(sin(Lat1) * sin(Lat2) +
        cos(Lat1) * cos(Lat2) * cos(Lon2-Lon1))
end  

function SizeGreatCircle(Edge::Grids.Edge,Nodes)
  P1 = Nodes[Edge.N[1]].P
  P2 = Nodes[Edge.N[2]].P
  return acos(dot(P1,P2)/(norm(P1)*norm(P2)))
end  

function SizeGreatCircle(P1::Point,P2::Point)
  return acos(dot(P1,P2)/(norm(P1)*norm(P2)))
end  


function AreaSphericalTriangle(P1,P2,P3)
  P1P2P3 = dot(P1,P2) + dot(P2,P3) + dot(P3,P1)
  P1_P2P3 = dot(P1,cross(P2,P3))
  area = 2.0 * atan(abs(P1_P2P3)  / (1.0 + P1P2P3))
end  

function AreaFace(Face,Nodes)
  P1 = Nodes[Face.N[1]].P
  Area = 0.0
  for i = 2 : length(Face.N) -1
    P2 = Nodes[Face.N[i]].P  
    P3 = Nodes[Face.N[i+1]].P  
    Area += AreaSphericalTriangle(P1,P2,P3)
  end  
  return Area
end

function RadiusFace(Face,Nodes)
  Radius = 0
  for N in Face.N
    Radius = max(Radius, SizeGreatCircle(GreatCircle(Face.Mid,Nodes[N].P)))
  end
  return Radius
end  

function inLeftHemisphere(P,C,offset)
  dot(cross(C.P1,C.P2),P) >= offset
end  


function CircumCenter(P1,P2,P3)
  C = cross(P2 - P1,P3 - P1)
  C = (1 / norm(C)) * C
end
