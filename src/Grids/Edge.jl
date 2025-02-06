mutable struct Edge
  N::Array{Int, 1}
  E::Int
  EG::Int
  EI::Int
  ET::Int
  NumF::Int
  F::Array{Int, 1}
  FG::Array{Int, 1}
  FP::Array{Int, 1}
  t::Point
  n::Point
  a::Float64
  Mid::Point
  FE::Array{Int,1}
  Type::String
  MasterSlave::Int
end

"""
  Edge()

This is my documentation
"""
function Edge()
  N=zeros(Int,2)
  E=0
  EG=0
  EI=0
  ET=0
  NumF=0
  F=zeros(Int,2)
  FG=zeros(Int,2)
  FP=zeros(Int,2)
  t=Point()
  n=Point()
  a=0.0
  Mid=Point()
  FE=zeros(Int,0)
  Type=""
  MasterSlave = 0
  return Edge(
    N,
    E,
    EG,
    EI,
    ET,
    NumF,
    F,
    FG,
    FP,
    t,
    n,
    a,
    Mid,
    FE,
    Type,
    MasterSlave,
  )
end  

function Edge(NodesE,Nodes,PosG,PosI,Type,PosT=nothing;Form="Cart",Rad=1.0)
  E = Edge()
  E.E=PosG;
  E.EI=PosI;
  if PosT ≠ nothing
    E.ET=PosT;
  else
    E.ET=0;
  end
  E.N=NodesE;
  E.t=Nodes[E.N[2]].P-Nodes[E.N[1]].P;
  if Form == "Sphere"
    E.a = SizeGreatCircle(Nodes[E.N[2]].P,Nodes[E.N[1]].P) * Rad  
  else    
    E.a=norm(E.t);
  end  
  E.t=E.t/norm(E.t)
  E.Mid=0.5*(Nodes[E.N[1]].P+Nodes[E.N[2]].P);
  if Form == "Sphere"
    E.Mid = E.Mid * (Rad / norm(E.Mid))  
  end  
  k = zeros(3)
  t = zeros(3)
  n = zeros(3)
  k[1] = E.Mid.x
  k[2] = E.Mid.y
  k[3] = E.Mid.z
  k = k / norm(k)
  t[1] = E.t.x
  t[2] = E.t.y
  t[3] = E.t.z
  n = LinearAlgebra.cross(t,k)
  E.n.x = n[1]
  E.n.y = n[2]
  E.n.z = n[3]
  E.Type=Type;
  return E
end

