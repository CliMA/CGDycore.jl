function Reordering!(Grid)
#  Order faces with HilbertOrdering
#  Old ordering is lexicographic per panel
#
  FaceOrderNew=zeros(Int,Grid.NumFaces)
  NodeOrderNew=zeros(Int,Grid.NumNodes)

end
function Renumbering(Grid)
for iF=1:Grid.NumFaces
  Grid.Faces[iF]=RenumberingFace4(Grid.Faces[iF],Grid);
end
for iE=1:Grid.NumEdges
  Grid.Edges[iE]=PosEdgeInFace(Grid.Edges[iE],Grid);
end
return Grid
end

function PosEdgeInFace(Edge,Grid)
Edge.FE=zeros(2)
for i=1:size(Edge.F,1)
  iF=Edge.F[i];
  if iF > 0
    for iE=1:4
      if Edge.EI==Grid.Edges[Grid.Faces[iF].E[iE]].EI
        Edge.FE[i]=iE;
        break
      end
    end
  end
end
if size(Edge.F,1)>1
  if Edge.FE[1] > Edge.FE[2]
    iTemp=Edge.FE[1];
    Edge.FE[1]=Edge.FE[2];
    Edge.FE[2]=iTemp;
    iTemp=Edge.F[1];
    Edge.F[1]=Edge.F[2];
    Edge.F[2]=iTemp;
    iTemp=Edge.FG[1];
    Edge.FG[1]=Edge.FG[2];
    Edge.FG[2]=iTemp;
  end
end
return Edge
end
function RenumberingFace4(Face,Grid)
local iN
for iN_in=1:length(Face.N)
  iN = iN_in
  N=Face.N[iN];
  num=0;
  for iE=1:length(Face.E)
    if N==Grid.Edges[Face.E[iE]].N[1]
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
  for i=1:length(Face.N)
    Face.N[i]=NTemp[i+iN-1];
    Face.E[i]=ETemp[i+iN-1];
    Face.P[i]=PTemp[i+iN-1];
  end
end
if length(Face.N) == 4
  OrientL=zeros(4,1);
  OrientL[1]=-1;
  OrientL[2]=1;
  OrientL[3]=1;
  OrientL[4]=-1;
  Face.OrientE = Vector{Int}(undef, 0)
  for i=1:4
    if Grid.Edges[Face.E[i]].N[1]==Face.N[i]
      push!(Face.OrientE, 1*OrientL[i]);
    else
      push!(Face.OrientE, -1*OrientL[i]);
    end
  end
end
return Face
end

