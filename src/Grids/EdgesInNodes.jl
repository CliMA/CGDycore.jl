function EdgesInNodes!(Nodes,Edges,Faces)

  NumNodes = size(Nodes,1)
  NumFaces = size(Faces,1)

  for iN = 1 : NumNodes
    EdgesInNode = Int64[]
    NumF = length(Nodes[iN].F) 
    for i = 1 : NumF
      iF = Nodes[iN].F[i]
      # Find Node in Face
      if iF > length(Faces)
        @show iN,iF,length(Faces)
      end  
      for j = 1 : length(Faces[iF].N)
        if iN == Faces[iF].N[j]
          push!(EdgesInNode,Faces[iF].E[j])
          exit
        end  
      end
    end  
    Nodes[iN].E = EdgesInNode
  end  
end
