#=
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
=#

mutable struct Node
    P::Point
    N::Int
    NG::Int
    E::Array{Int, 1}
    F::Array{Int, 1}
    FG::Array{Int, 1}
    FP::Array{Int, 1}
    MasterSlave::Int
end

function Node()
  P=Point()
  N=0
  NG=0
  E=zeros(Int,0)
  F=zeros(Int,0)
  FG=zeros(Int,0)
  FP=zeros(Int,0)
  MasterSlave = 0
  return Node(
    P,
    N,
    NG,
    E,
    F,
    FG,
    FP,
    MasterSlave,
  )
end  
  
function Node(Point::Point, Pos::Int)
  N=Node()
  N.P=Point
  N.N=Pos
  return N
end
