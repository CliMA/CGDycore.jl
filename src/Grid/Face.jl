mutable struct Face
  N::Array{Int, 1}
  E::Array{Int, 1}
  F::Int
  FG::Int
  n::Point
  a::Float64
  Mid::Point
  OrientE::Array{Int, 1}
  Type::String
  P::Array{Point, 1}
end

function Face()
  N=zeros(Int,0)
  E=zeros(Int,0)
  F=0
  FG=0
  n=Point()
  a=0.0
  Mid=Point()
  OrientE=zeros(Int,0)
  Type=""
  P=Array{Point}(undef, 0)
  return Face(
    N,
    E,
    F,
    FG,
    n,
    a,
    Mid,
    OrientE,
    Type,
    P,
  )
end  

function Face(Edges::Array{Int, 1},Grid,Pos,Type,OrientFace;P::Array{Float64,2}=[])
  F = Face()
  if Edges[1]==0
    return (F,Grid)
  end

  nE=size(Edges,1);
  F.F=Pos;
  F.Type=Type
  # TODO: check translation
  for iE=1:nE
    Grid.Edges[Edges[iE]].NumF +=1  
    Grid.Edges[Edges[iE]].F[Grid.Edges[Edges[iE]].NumF]=Pos;
  end
  #Sort edges
  F.E=zeros(Int,nE);
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
  F.N=zeros(Int,nE);
  F.N[1:2]=Grid.Edges[F.E[1]].N;
  for iE=2:nE-1
    if F.N[iE]==Grid.Edges[F.E[iE]].N[1]
      F.N[iE+1]=Grid.Edges[F.E[iE]].N[2];
    else
      F.N[iE+1]=Grid.Edges[F.E[iE]].N[1];
    end
  end
  if P == zeros(Float64,0,0)
    F.P=Array{Point}(undef, size(F.N,1))  
    for i=1:size(F.N,1)
      F.P[i]=Grid.Nodes[F.N[i]].P;
    end
  else
    F.P=Array{Point}(undef, size(F.N,1))  
    for i=1:size(F.N,1)
      F.P[i]=Point(P[:,i])
    end
  end
  PT=Point([0.0, 0.0, 0.0]);
  for i=1:nE-1
    PT=PT+cross(F.P[i],F.P[i+1]);
  end
  PT=PT+cross(F.P[nE],F.P[1]);
  F.a=0.5*norm(PT);
  for i=1:nE
    F.Mid=F.Mid+F.P[i];
  end
  F.Mid=F.Mid/Float64(nE);

  NumE=size(Edges,1);
  F.n=cross(F.P[NumE],F.P[1]);
  for i=1:NumE-1
    F.n=F.n+cross(F.P[i],F.P[i+1]);
  end
  F.n=F.n/norm(F.n);
  if OrientFace(F.n,F.Mid) < 0 
    #Change Orientation
    NTemp=copy(F.N);
    ETemp=copy(F.E);
    PTemp=copy(F.P);
    for i=1:nE
      F.N[i]=NTemp[nE-i+1];
      F.P[i]=PTemp[nE-i+1];
    end
    for i=1:nE-1
      F.E[i]=ETemp[nE-i];
    end
    F.n=-F.n;
  end
  return (F,Grid)
end
