import Base: -
import Base: /
import Base: +
import Base: *
import LinearAlgebra: norm
import LinearAlgebra: dot

struct Point
# P::Array{Float64,1}
  x::Float64
  y::Float64
  z::Float64
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


function sphere2cart(lam,phi,r)
x=cos(lam)*cos(phi)*r;
y=sin(lam)*cos(phi)*r;
z=sin(phi)*r;
return [x;y;z];
end

function sphereDeg2cart(lam,phi,r)
x=cosd(lam)*cosd(phi)*r;
y=sind(lam)*cosd(phi)*r;
z=sind(phi)*r;
return [x;y;z];
end


function cart2sphere(x,y,z)
r=sqrt(x^2+y^2+z^2);
phi=asin(z/r);

# ϕ = atan(z, hypot(y, x))
# if abs(ϕ) == 90
#     λ = zero(ϕ)
# else
#     λ = atan(y, x)
# end

lam=0;
if abs(abs(phi)-pi/2)>1.e-14
  lam=atan(y,x); # TODO: check translation with Oswald
  if lam<0.0
    lam=lam+2*pi;
  end
end
#return (float(lam),float(phi),r)
return lam,phi,r
end

function cart2sphereDeg(x,y,z)
r=sqrt(x^2+y^2+z^2);
phi=asind(z/r);

# ϕ = atan(z, hypot(y, x))
# if abs(ϕ) == 90
#     λ = zero(ϕ)
# else
#     λ = atan(y, x)
# end

lam=-180;
if abs(abs(phi)-90.0)>1.e-14
  lam=atand(y,x); # TODO: check translation with Oswald
end
return (float(lam),float(phi),r)
end


function GreatCircle(Lon1,Lat1,Lon2,Lat2) 
  return acos(sin(Lat1) * sin(Lat2) +
        cos(Lat1) * cos(Lat2) * cos(Lon2-Lon1))
end  

