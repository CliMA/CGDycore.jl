function [Grid] = FacesInEdges(Grid)
NumFacesPerEdge=zeros(Grid.NumEdges,1);
for iF=1:Grid.NumFaces
  Face=Grid.Faces(iF);
  for iE=1:size(Face.E,2)
    NumFacesPerEdge(Face.E(iE))=NumFacesPerEdge(Face.E(iE))+1;
  end
end
FacesPerEdge=zeros(Grid.NumEdges,max(NumFacesPerEdge));
NumFacesPerEdge=zeros(Grid.NumEdges,1);
for iF=1:Grid.NumFaces
  Face=Grid.Faces(iF);
  for iE=1:size(Face.E,2)
    NumFacesPerEdge(Face.E(iE))=NumFacesPerEdge(Face.E(iE))+1;
    FacesPerEdge(Face.E(iE),NumFacesPerEdge(Face.E(iE)))=iF;
  end
end
for iE=1:Grid.NumEdges
  Grid.Edges(iE).F=FacesPerEdge(iE,1:NumFacesPerEdge(iE));
end
end
