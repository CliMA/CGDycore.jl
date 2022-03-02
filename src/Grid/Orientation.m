function [Grid] = Orientation(Grid)
NumEdges=size(Grid.Edges,2);
VisitE=zeros(1,NumEdges);
OrientationE=zeros(1,NumEdges);
o=1;
for iE=1:NumEdges
  [VisitE,OrientationE,o,Grid]=Orient(iE,o,VisitE,OrientationE,Grid);
end
end

function [VisitE,OrientationE,o,Grid] = Orient(iE,o,VisitE,OrientationE,Grid)
if VisitE(iE)==1 || size(Grid.Edges(iE).F,2)==0
%   if OrientationE(iE) ~= o
%     error('Möbius Band')
%   end
else
  VisitE(iE)=1;
  iF1=Grid.Edges(iE).F(1);
  if size(Grid.Edges(iE).F,2)==2
    iFEnd=Grid.Edges(iE).F(2);
  else
    iFEnd=0;
  end
  for iE1=1:4
    if iE==Grid.Faces(iF1).E(iE1)
      break
    end
  end
  NF=[Grid.Faces(iF1).N,Grid.Faces(iF1).N(1)];
  o=EdgeDirection(NF(iE1),NF(iE1+1),Grid.Edges(iE).N);
  OrientationE(iE)=o;
  iLoop=0;
  while (iF1~=iFEnd && size(Grid.Edges(iE).F,2)==2) || iLoop==0
    iLoop=iLoop+1;
    NF=[Grid.Faces(iF1).N,Grid.Faces(iF1).N(1)];
    o=EdgeDirection(NF(iE1),NF(iE1+1),Grid.Edges(iE).N);
    iEP1=Shift(iE1,2);
    iEP=Grid.Faces(iF1).E(iEP1);
    if EdgeDirection(NF(iEP1),NF(iEP1+1),Grid.Edges(iEP).N)==o
      if o==-1
        Grid.Edges(iEP).N(1)=NF(iEP1);
        Grid.Edges(iEP).N(2)=NF(iEP1+1);
      else
        Grid.Edges(iEP).N(1)=NF(iEP1+1);
        Grid.Edges(iEP).N(2)=NF(iEP1);
      end
    end
    iE=iEP;
    VisitE(iE)=1;
    OrientationE(iE)=o;
    if size(Grid.Edges(iE).F,2)>1
      if iF1==Grid.Edges(iE).F(1)
        iF1=Grid.Edges(iE).F(2);
      else
        iF1=Grid.Edges(iE).F(1);
      end
    else
      iF1=Grid.Edges(iE).F(1);
    end
    for iE1=1:4
      if iE==Grid.Faces(iF1).E(iE1)
        break
      end
    end
  end
  
  if iF1==iFEnd
    NF=[Grid.Faces(iF1).N,Grid.Faces(iF1).N(1)];
    o=EdgeDirection(NF(iE1),NF(iE1+1),Grid.Edges(iE).N);
    iEP1=Shift(iE1,2);
    iEP=Grid.Faces(iF1).E(iEP1);
    oP=EdgeDirection(NF(iEP1),NF(iEP1+1),Grid.Edges(iEP).N);
    if iEP==Grid.Edges(iEP).E
      if o==oP
        error('Möbius Band')
      end
    else
      if o==oP
        error('Möbius Band')
      end
    end
  end
end
end


function Dir=EdgeDirection(N1,N2,E)
if N1==E(1)
  Dir=1;
else
  Dir=-1;
end
end

function iS=Shift(i,S)
iS=i+S;
if iS>4
  iS=iS-4;
end
end

