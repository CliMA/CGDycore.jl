Base.@kwdef mutable struct Edge
  N = nothing
  E = nothing
  EI = nothing
  ET = nothing
  F = nothing
  C = nothing
  t = nothing
  a = nothing
  Mid = nothing
  FE = nothing
  Type = nothing
  DGPoints = nothing
end

function Edge(Nodes,Grid,PosG,PosI,Type,PosT=nothing)
  E = Edge(;)
  E.E=PosG;
  E.EI=PosI;
  if PosT ≠ nothing
    E.ET=PosT;
  else
    E.ET=0;
  end
  E.N=Nodes;
  E.F=zeros(Int, 1,0);
  E.t=Grid.Nodes[E.N[2]].P-Grid.Nodes[E.N[1]].P;
  E.a=norm(E.t);
  E.t=E.t/E.a;
  E.Mid=0.5*(Grid.Nodes[E.N[1]].P+Grid.Nodes[E.N[2]].P);
  E.Type=Type;
  return E
end

