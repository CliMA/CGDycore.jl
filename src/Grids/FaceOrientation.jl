function FaceOrientation(Face,Grid)
  NumN = length(Face.N)
  P1 = Grid.Nodes[Face.N[NumN]].P
  P2 = Grid.Nodes[Face.N[1]].P
  n = cross(P1,P2)
  @inbounds for i = 1 : NumN - 1
    P1 = Grid.Nodes[Face.N[i]].P
    P2 = Grid.Nodes[Face.N[i+1]].P
    n = n + cross(P1,P2)
  end
  n = n / norm(n)
  if dot(n,Face.Mid) < 0.0
    Orientation = -1
  else
    Orientation = 1  
  end  
  return Orientation
end  
