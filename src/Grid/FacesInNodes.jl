function FacesInNodes(Grid)
NumFacesPerNode=zeros(Int,Grid.NumNodes,1);
for iF=1:Grid.NumFaces
  Face=Grid.Faces[iF];
  for iN=1:size(Face.N,1)
    NumFacesPerNode[Face.N[iN]]=NumFacesPerNode[Face.N[iN]]+1;
  end
end
FacesPerNode=zeros(Int,Grid.NumNodes,maximum(NumFacesPerNode));
NumFacesPerNode=zeros(Int,Grid.NumNodes,1);
for iF=1:Grid.NumFaces
  Face=Grid.Faces[iF];
  for iN=1:size(Face.N,1)
    NumFacesPerNode[Face.N[iN]]=NumFacesPerNode[Face.N[iN]]+1;
    FacesPerNode[Face.N[iN],NumFacesPerNode[Face.N[iN]]]=iF;
#   if iN == 1
#     Ind1[Face.N[iN],NumFacesPerNode[Face.N[iN]]]=1
#     Ind2[Face.N[iN],NumFacesPerNode[Face.N[iN]]]=1
#   elseif iN == 2
#     Ind1[Face.N[iN],NumFacesPerNode[Face.N[iN]]]=OrdPoly+1
#     Ind2[Face.N[iN],NumFacesPerNode[Face.N[iN]]]=1
#   elseif iN == 3
#     Ind1[Face.N[iN],NumFacesPerNode[Face.N[iN]]]=OrdPoly+1
#     Ind2[Face.N[iN],NumFacesPerNode[Face.N[iN]]]=OrdPoly+1
#   elseif iN == 4
#     Ind1[Face.N[iN],NumFacesPerNode[Face.N[iN]]]=1
#     Ind2[Face.N[iN],NumFacesPerNode[Face.N[iN]]]=OrdPoly+1
#   end
  end
end
for iN=1:Grid.NumNodes
  Grid.Nodes[iN].F=FacesPerNode[iN,1:NumFacesPerNode[iN]];
end
return Grid
end
