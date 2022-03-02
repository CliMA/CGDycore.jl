function Grid=SphereGrid(OrientFace,Param)
Grid.nBar=[ 0  1   0   1 
           -1  0  -1   0];
  Grid.Dim=3;
  Grid.Type='Quad';
  Grid.Rad=Param.RadEarth;
  Grid.Form='Sphere';
  lat1=0;
  lat2=lat1+pi/1000;
  lon1=0;
  lon2=lon1+pi/1000;
  
  NumNodes=4;
  Nodes(1:NumNodes)=Node([0 0 0],0);
  NodeNumber=1;
  Nodes(NodeNumber)=Node(sphere2cart(lat1,lon1,Grid.Rad),NodeNumber);
  NodeNumber=2;
  Nodes(NodeNumber)=Node(sphere2cart(lat2,lon1,Grid.Rad),NodeNumber);
  NodeNumber=3;
  Nodes(NodeNumber)=Node(sphere2cart(lat2,lon2,Grid.Rad),NodeNumber);
  NodeNumber=4;
  Nodes(NodeNumber)=Node(sphere2cart(lat1,lon2,Grid.Rad),NodeNumber);
  Grid.Nodes=Nodes;
  
  NumEdges=4;
  NumEdgesI=4;
  Edges(1:NumEdges)=Edge([1 2],Grid,0,0,'');
  EdgeNumber=1;
  Edges(EdgeNumber)=Edge([1 2],Grid,EdgeNumber,EdgeNumber,'');
  EdgeNumber=2;
  Edges(EdgeNumber)=Edge([2 3],Grid,EdgeNumber,EdgeNumber,'');
  EdgeNumber=3;
  Edges(EdgeNumber)=Edge([3 4],Grid,EdgeNumber,EdgeNumber,'');
  EdgeNumber=4;
  Edges(EdgeNumber)=Edge([4 1],Grid,EdgeNumber,EdgeNumber,'');
  Grid.Edges=Edges;
  
  NumFaces=1;
  FaceNumber=1;
  [Faces(FaceNumber),Grid]=Face([1 2 3 4],Grid,FaceNumber,'',OrientFace);
  Grid.Faces=Faces;
  
  Grid.NumNodes=size(Grid.Nodes,2);
  Grid.NumEdges=size(Grid.Edges,2);
  Grid.NumEdgesI=size(Grid.Edges,2);
  Grid.NumEdgesB=0;
  Grid.NumFaces=size(Grid.Faces,2);
  Grid.Dim=3;
  Grid=Orientation(Grid);
  Grid=Renumbering(Grid);
  Grid=FacesInNodes(Grid);
  

end
