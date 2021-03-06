Base.@kwdef mutable struct Face
    N = nothing
    NG = nothing
    E = nothing
    F = nothing
    C = nothing
    n = nothing
    a = nothing
    Mid = nothing
    OrientE = nothing
    Dir = nothing
    Ori = nothing
    SF = nothing
    Kite = nothing
    DF = nothing
    DFInv = nothing
    J = nothing
    lon = nothing
    lat = nothing
    P = nothing
end

function Face(Edges,Grid,Pos,Dir,OrientFace=nothing,P=nothing)
  F = Face(;)
  if Edges[1]==0
    return (F,Grid)
  end
  @assert !isnothing(OrientFace)

  F.DF=zeros(0,0);
  nE=size(Edges,2);
  F.F=Pos;
  F.Dir=Dir;
  F.C=zeros(1,2);
  # TODO: check translation
  for iE=1:nE
    Grid.Edges[Edges[iE]].F=[Grid.Edges[Edges[iE]].F Pos];
  end
  #Sort edges
  F.E=zeros(Int,1,nE);
  F.E[1]=Edges[1];
  N2=Grid.Edges[F.E[1]].N[2];
  for iE=2:nE
    for iE1=iE:nE
      if N2==Grid.Edges[Edges[iE1]].N[1]
        F.E[iE]=Edges[iE1];
        N2=Grid.Edges[Edges[iE1]].N[2];
        Edges[iE1]=Edges[iE];
        Edges[iE]=F.E[iE];
        break
      elseif N2==Grid.Edges[Edges[iE1]].N[2]
        F.E[iE]=Edges[iE1];
        N2=Grid.Edges[Edges[iE1]].N[1];
        Edges[iE1]=Edges[iE];
        Edges[iE]=F.E[iE];
        break
      end
    end
  end
  F.a=0;
  F.N=zeros(Int, 1,nE);
  F.N[1,1:2]=Grid.Edges[F.E[1]].N;
  for iE=2:nE-1
    if F.N[1,iE]==Grid.Edges[F.E[iE]].N[1]
      F.N[1,iE+1]=Grid.Edges[F.E[iE]].N[2];
    else
      F.N[1,iE+1]=Grid.Edges[F.E[iE]].N[1];
    end
  end
  if P == nothing
    F.P=zeros(3,size(F.N,2));
    for i=1:size(F.N,2)
      F.P[:,i]=Grid.Nodes[F.N[i]].P';
    end
  else
    F.P=P;
  end
  PT=[0, 0, 0];
  for i=1:nE-1
    #P=P+cross(Grid.Nodes(F.N(i)).P,Grid.Nodes(F.N(i+1)).P);
    PT=PT+cross(F.P[:,i],F.P[:,i+1]);
  end
  #P=P+cross(Grid.Nodes(F.N(nE)).P,Grid.Nodes(F.N[1]).P);
  PT=PT+cross(F.P[:,nE],F.P[:,1]);
  F.a=0.5*norm(PT);
  F.Mid=zeros(3);
  for i=1:nE
    #F.Mid=F.Mid+Grid.Nodes(F.N(i)).P;
    F.Mid=F.Mid+F.P[:,i];
  end
  F.Mid=F.Mid/nE;

  NumE=size(Edges,2);
  #F.n=cross(Grid.Nodes(F.N(NumE)).P,Grid.Nodes(F.N[1]).P);
  F.n=cross(F.P[:,NumE],F.P[:,1]);
  for i=1:NumE-1
    #F.n=F.n+cross(Grid.Nodes(F.N(i)).P,Grid.Nodes(F.N(i+1)).P);
    F.n=F.n+cross(F.P[:,i],F.P[:,i+1]);
  end
  F.n=F.n/norm(F.n);
  if OrientFace(F.n',F.Mid)<0 # TODO: check translation with Oswald
    #Change Orientation
    NTemp=copy(F.N);
    ETemp=copy(F.E);
    PTemp=copy(F.P);
    for i=1:nE
      F.N[i]=NTemp[nE-i+1];
      F.P[:,i]=PTemp[:,nE-i+1];
    end
    for i=1:nE-1
      F.E[i]=ETemp[nE-i];
    end
    F.n=-F.n;
  end
  return (F,Grid)
end
