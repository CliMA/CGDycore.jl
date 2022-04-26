function ConnectivityGraph(Grid)
  I=[]
  J=[]
  V=[]
  for iF=1:Grid.NumFaces
    for iE in Grid.Faces[iF].E
      for iFF in Grid.Edges[iE].F
        if iFF != iF
          push!(I,iF)
          push!(J,iFF)
          push!(V,1)
        end
      end
    end
    for iN in Grid.Faces[iF].N
      for iFF in Grid.Nodes[iN].F
        if iFF != iF
          push!(I,iF)
          push!(J,iFF)
          push!(V,1)
        end
      end
    end
  end
  return sparse(I,J,V)
end  

function Coloring1(Grid)
  colors = Array{Array{Int, 1}, 1}(undef, 0)
  color = 1:Grid.NumFaces
  push!(colors,color)
  return colors
end

function Coloring(Grid)
  IndFace=ones(Grid.NumFaces,1);
  for k=1:10
    for iF=1:Grid.NumFaces
      if IndFace[iF]==k
        for iN in Grid.Faces[iF].N
          for iFF in Grid.Nodes[iN].F
            if iFF!=iF && IndFace[iFF]>=k
              IndFace[iFF]=k+1;
            end
          end
        end
      end
    end
  end
  iMax=maximum(IndFace);
  colors = Array{Array{Int, 1}, 1}(undef, 0)
  for i=1:iMax
    color=[]  
    for iF = 1:Grid.NumFaces  
      if i == IndFace[iF]
        push!(color,iF)
      end
    end
    push!(colors,color)
  end  
  return colors
end
