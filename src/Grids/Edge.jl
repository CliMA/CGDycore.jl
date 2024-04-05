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
  a::Float64
  Mid::Point
  FE::Array{Int,1}
  Type::String
  MasterSlave::Int
end

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
  E.t=E.t/E.a;
  E.Mid=0.5*(Nodes[E.N[1]].P+Nodes[E.N[2]].P);
  if Form == "Sphere"
    E.Mid = E.Mid * (Rad / norm(E.Mid))  
  end  
  E.Type=Type;
  return E
end

