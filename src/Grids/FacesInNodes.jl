function FacesInNodes!(Nodes,Faces)
  NumNodes = size(Nodes,1)
  NumFaces = size(Faces,1)
  NumFacesPerNode = zeros(Int,NumNodes,1)
  for iF = 1 : NumFaces
    Face = Faces[iF]
    for iN = 1 : size(Face.N,1)
      NumFacesPerNode[Face.N[iN]]=NumFacesPerNode[Face.N[iN]]+1
    end
  end
  FacesPerNode = zeros(Int,NumNodes,maximum(NumFacesPerNode))
  NumFacesPerNode = zeros(Int,NumNodes,1)
  for iF = 1 : NumFaces
    Face = Faces[iF]
    for iN = 1 : size(Face.N,1)
      NumFacesPerNode[Face.N[iN]] = NumFacesPerNode[Face.N[iN]]+1
      FacesPerNode[Face.N[iN],NumFacesPerNode[Face.N[iN]]] = iF
    end
  end
  for iN = 1 : NumNodes
    Nodes[iN].F = FacesPerNode[iN,1:NumFacesPerNode[iN]]
  end
end

function SortFacesInNodes!(Nodes,Faces)
  NumNodes = size(Nodes,1)

  for iN = 1 : NumNodes
    F = deepcopy(Nodes[iN].F)
    lenF = length(F)
    NPr = zeros(Int,lenF)
    NSu = zeros(Int,lenF)
    for i = 1 : lenF
      iF = F[i]
      for j = 1 : length(Faces[iF].N)
        jN = Faces[iF].N[j]
        if jN == iN
          if j > 1
            NPr[i] = Faces[iF].N[j-1]  
          else
            NPr[i] = Faces[iF].N[end]  
          end
          if j < length(Faces[iF].N)
            NSu[i] = Faces[iF].N[j+1]  
          else
            NSu[i] = Faces[iF].N[1]  
          end  
        end
      end
    end   
    if Nodes[iN].Type == 'B' || Nodes[iN].Type == 'P'  
      NSuStart = setdiff(NSu,NPr)
      for i = 1 : lenF
        if NSu[i] in NSuStart  
          NextSu = NPr[i]
          Nodes[iN].F[1] = F[i]
        end  
      end  
      for i = 1 : lenF
        for j = 1 : lenF
          if NextSu == NSu[j]  
            Nodes[iN].F[i+1] = F[j]  
            NextSu = NPr[j]
            break
          end
        end
      end  
    else    
      NextSu = NPr[end]
      for i = 1 : lenF
        for j = 1 : lenF
          if NextSu == NSu[j]  
            Nodes[iN].F[i] = F[j]  
            NextSu = NPr[j]
            break
          end
        end
      end  
    end  
  end
end
