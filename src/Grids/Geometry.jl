import Base: -
import Base: /
import Base: +
import Base: *
import LinearAlgebra: norm
import LinearAlgebra: dot

struct Point
  x::Float64
  y::Float64
  z::Float64
end

mutable struct GreatCircle
  P1::Point
  P2::Point
end


function Point()
  return Point(
    0.0,
    0.0,
    0.0,
  )
end

function Point(P::Array{Float64, 1})
  return Point(
    P[1],
    P[2],
    P[3],
  )
end
-(P::Point)=Point([-P.x,-P.y,-P.z])
-(P1::Point,P2::Point)=Point([P1.x-P2.x,P1.y-P2.y,P1.z-P2.z])
+(P1::Point,P2::Point)=Point([P1.x+P2.x,P1.y+P2.y,P1.z+P2.z])
/(P::Point,s::Float64)=Point([P.x/s,P.y/s,P.z/s])
*(s::Float64,P::Point)=Point([s*P.x,s*P.y,s*P.z])
*(P::Point,s::Float64)=Point([s*P.x,s*P.y,s*P.z])
norm(P::Point)=sqrt(P.x*P.x + P.y*P.y + P.z*P.z)
dot(P1::Point,P2::Point)=P1.x*P2.x + P1.y*P2.y + P1.z*P2.z
cross(P1::Point,P2::Point)=Point([P1.y*P2.z-P1.z*P2.y,P1.z*P2.x-P1.x*P2.z,P1.x*P2.y-P1.y*P2.x])

function normalize(P)
  P = P / norm(P)
end  


function sphere2cart(lam,phi,r)
x=cos(lam)*cos(phi)*r
y=sin(lam)*cos(phi)*r
z=sin(phi)*r
return [xy;z];
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

function SizeGreatCircle(C) 
  return acos(dot(C.P1,C.P2)/(norm(P1)*norm(P2)))
end  

function GreatCircle(Lon1,Lat1,Lon2,Lat2) 
  return acos(sin(Lat1) * sin(Lat2) +
        cos(Lat1) * cos(Lat2) * cos(Lon2-Lon1))
end  

function AreaSphericalTriangle(P1,P2,P3)
  P1Loc = P1 / norm(P1)
  P2Loc = P2 / norm(P2)
  P3Loc = P3 / norm(P3)
  P1P2P3 = dot(P1Loc,P2Loc) + dot(P2Loc,P3Loc) + dot(P3Loc,P1Loc)
  P1_P2P3 = dot(P1Loc,cross(P2Loc,P3Loc))
  area = 2.0 * atan(abs(P1_P2P3)  / (1.0 + P1P2P3))
end  

function AreaFace(Face,Nodes)
  P1 = Nodes[Face.N[1]].P
  Area = 0.0
  for i = 2 : length(Face.N) -1
    P2 = Nodes[Face.N[i]].P  
    P3 = Nodes[Face.N[i+1]].P  
    Area += AreaSphericalTriangle(P1,P2,P3)
    @show Area
  end  
  stop
  return Area
end
function RadiusFace(Face,Nodes)
  Radius = 0
  for N in Face.N
    Radius = max(Radius, GreatCircle(Face.Mid,Nodes[N].P))
  end
  return Radius
end  


