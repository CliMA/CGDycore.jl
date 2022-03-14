classdef Node
  
  properties
    P
    N
    E
    F
  end
  
  methods
    function N=Node(Point,Pos)
      N.N=Pos;
      N.P=Point';
    end
  end
end

