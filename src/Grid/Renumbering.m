function [Grid] = Renumbering(Grid)
for iF=1:Grid.NumFaces
  Grid.Faces(iF)=RenumberingFace4(Grid.Faces(iF),Grid);
end
for iE=1:Grid.NumEdges
  Grid.Edges(iE)=PosEdgeInFace(Grid.Edges(iE),Grid);
end
end

function Edge=PosEdgeInFace(Edge,Grid)
Edge.FE=zeros(1,2);
for iF=1:size(Edge.F,2)
  F=Edge.F(iF);
  for iE=1:4
    if Edge.EI==Grid.Edges(Grid.Faces(F).E(iE)).EI
      Edge.FE(iF)=iE;
      break
    end
  end
end
if size(Edge.F,2)>1
  if Edge.FE(1) > Edge.FE(2)
    iTemp=Edge.FE(1);
    Edge.FE(1)=Edge.FE(2);
    Edge.FE(2)=iTemp;
    iTemp=Edge.F(1);
    Edge.F(1)=Edge.F(2);
    Edge.F(2)=iTemp;
  end
end
end
function Face=RenumberingFace4(Face,Grid)
for iN=1:4
  N=Face.N(iN);
  num=0;
  for iE=1:4
    if N==Grid.Edges(Face.E(iE)).N(1)
      num=num+1;
    end
    if num==2
      break
    end
  end
  if num==2
    break
  end
end
if iN>1
  NTemp=[Face.N Face.N];
  ETemp=[Face.E Face.E];
  PTemp=[Face.P Face.P];
  for i=1:4
    Face.N(i)=NTemp(i+iN-1);
    Face.E(i)=ETemp(i+iN-1);
    Face.P(:,i)=PTemp(:,i+iN-1);
  end
end
OrientL=zeros(4,1);
OrientL(1)=-1;
OrientL(2)=1;
OrientL(3)=1;
OrientL(4)=-1;
for i=1:4
  if Grid.Edges(Face.E(i)).N(1)==Face.N(i)
    Face.OrientE(i)=1*OrientL(i);
  else
    Face.OrientE(i)=-1*OrientL(i);
  end
end
end

