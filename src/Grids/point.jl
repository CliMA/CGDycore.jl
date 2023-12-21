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
