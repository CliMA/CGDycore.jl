function FacesInNodes(Grid)
NumFacesPerNode=zeros(Int,Grid.NumNodes,1);
for iF=1:Grid.NumFaces
  Face=Grid.Faces[iF];
  for iN=1:size(Face.N,2)
    NumFacesPerNode[Face.N[iN]]=NumFacesPerNode[Face.N[iN]]+1;
  end
end
FacesPerNode=zeros(Int,Grid.NumNodes,maximum(NumFacesPerNode));
NumFacesPerNode=zeros(Int,Grid.NumNodes,1);
for iF=1:Grid.NumFaces
  Face=Grid.Faces[iF];
  for iN=1:size(Face.N,2)
    NumFacesPerNode[Face.N[iN]]=NumFacesPerNode[Face.N[iN]]+1;
    FacesPerNode[Face.N[iN],NumFacesPerNode[Face.N[iN]]]=iF;
  end
end
for iN=1:Grid.NumNodes
  Grid.Nodes[iN].F=FacesPerNode[iN,1:NumFacesPerNode[iN]];
end
return Grid
end
