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
