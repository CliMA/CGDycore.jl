classdef Edge
  
  properties
    N
    E
    EI
    ET
    F
    C
    t
    a
    Mid
    FE
    Type
    DGPoints
  end
  
  methods
    function E=Edge(Nodes,Grid,PosG,PosI,Type,PosT)
      E.E=PosG;
      E.EI=PosI;
      if nargin>5
        E.ET=PosT;
      else
        E.ET=0;
      end
      E.N=Nodes;
      E.F=zeros(1,0);
      E.t=Grid.Nodes(E.N(2)).P-Grid.Nodes(E.N(1)).P;
      E.a=norm(E.t);
      E.t=E.t/E.a;
      E.Mid=0.5*(Grid.Nodes(E.N(1)).P+Grid.Nodes(E.N(2)).P);
      E.Type=Type;
    end
  end
end

