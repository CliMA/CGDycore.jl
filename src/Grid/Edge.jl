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

function Edge(Nodes,Grid,PosG,PosI,Type,PosT=nothing)
  E = Edge()
  E.E=PosG;
  E.EI=PosI;
  if PosT ≠ nothing
    E.ET=PosT;
  else
    E.ET=0;
  end
  E.N=Nodes;
  E.t=Grid.Nodes[E.N[2]].P-Grid.Nodes[E.N[1]].P;
  E.a=norm(E.t);
  E.t=E.t/E.a;
  E.Mid=0.5*(Grid.Nodes[E.N[1]].P+Grid.Nodes[E.N[2]].P);
  E.Type=Type;
  return E
end

