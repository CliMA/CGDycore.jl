function TestGrid(Grid)

# Orientation
  Orientation1 = Grid.Faces[1].Orientation
  for iF = 2 : length(Grid.Faces)
    if Grid.Faces[iF].Orientation == Orientation1
    else
      @show "Changed orienation",iF
    end  
  end    
# Arrangement of faces in node structure   
  for iN = 1 : Grid.NumNodes
    Np = zeros(Int,length(Grid.Nodes[iN].F))  
    Nm = zeros(Int,length(Grid.Nodes[iN].F))  
    for i = 1 : length(Grid.Nodes[iN].F)
      iF = Grid.Nodes[iN].F[i]  
      # Determine position of iN faces iF  
      for j = 1 : length(Grid.Faces[iF].N)
        if iN == Grid.Faces[iF].N[j]
          if j < length(Grid.Faces[iF].N)  
            Np[i] = Grid.Faces[iF].N[j+1]
          else
            Np[i] = Grid.Faces[iF].N[1]  
          end
          if j > 1
            Nm[i] = Grid.Faces[iF].N[j-1]
          else
            Nm[i] = Grid.Faces[iF].N[end]  
          end
          exit
        end
      end  
    end     
    Np1 = Np[1]
    @. Np[1:end-1] = Np[2:end]
    Np[end] = Np1
    if Np == Nm
    else
      @show iN  
      @show Nm  
      @show Np  
      stop
    end  
  end     
# Test nodes and edges in faces  
  for iF = 1 : Grid.NumFaces
    for i = 1 : length(Grid.Faces[iF].N)
      iN = Grid.Faces[iF].N[i]  
      iE = Grid.Faces[iF].E[i]
      if iN in Grid.Edges[iE].N
      else
        @show iN
        @show Grid.Edges[iE].N
        stop
      end  
    end  
  end  
  @show "All test successful"
end
