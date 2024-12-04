function OrientTriangle(Grid)
  for iF = 1 : Grid.NumFaces
    F = Grid.Faces[iF]
    N2 = [F.N;F.N]
    E2 = [F.E;F.E]
    OrientE2 = [F.OrientE;F.OrientE]
    
    #edges
    NE1 = Grid.Edges[F.E[1]].N[1]
    NE2 = Grid.Edges[F.E[2]].N[1]
    NE3 = Grid.Edges[F.E[3]].N[1]

    if NE1 == NE2
      NS = NE1
    elseif NE1 == NE3
      NS = NE1
    else
      NS = NE2
    end 

    if NS == F.N[1]
      s = 0
    elseif NS == F.N[2]
      s = 1
    else 
      s = 2
    end

    for i = 1 : 3
      F.N[i] = N2[s+i]
      F.E[i] = E2[s+i]
      F.OrientE[i] = OrientE2[s+i]
    end
  end
end