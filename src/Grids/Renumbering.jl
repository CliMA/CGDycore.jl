function Reordering!(Grid)
#  Order faces with HilbertOrdering
#  Old ordering is lexicographic per panel
#
  FaceOrderNew=zeros(Int,Grid.NumFaces)
  NodeOrderNew=zeros(Int,Grid.NumNodes)

end

function Renumbering!(Edges,Faces)
  for iF = 1 : size(Faces,1)
    RenumberingFace4!(Faces[iF],Edges)
  end
  for iE = 1 : size(Edges,1)
    PosEdgeInFace!(Edges[iE],Edges,Faces)
  end
end

function PosEdgeInFace!(Edge,Edges,Faces)
  Edge.FE = zeros(2)
  for i = 1:size(Edge.F,1)
    iF = Edge.F[i]
    if iF > 0
      for iE = 1 : 4
        if Edge.EI == Edges[Faces[iF].E[iE]].EI
          Edge.FE[i] = iE
          break
        end
      end
    end
  end
  if size(Edge.F,1)>1
    iF = Edge.F[1]  
    if iF > 0
      EdgeType = Edge.FE[1]  
      if Faces[iF].Orientation * Faces[iF].OrientE[EdgeType] == -1
        iTemp = Edge.FE[1]
        Edge.FE[1] = Edge.FE[2]
        Edge.FE[2] = iTemp
        iTemp = Edge.F[1]
        Edge.F[1] = Edge.F[2]
        Edge.F[2] = iTemp
        iTemp = Edge.FG[1]
        Edge.FG[1] = Edge.FG[2]
        Edge.FG[2] = iTemp
      end
    end
  end
end

function RenumberingFace4!(Face,Edges)
local iN
for iN_in=1:length(Face.N)
  iN = iN_in
  N=Face.N[iN]
  num=0
  for iE=1:length(Face.E)
    if N==Edges[Face.E[iE]].N[1]
      num=num+1
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
  NTemp=[Face.N Face.N]
  ETemp=[Face.E Face.E]
  PTemp=[Face.P Face.P]
  for i=1:length(Face.N)
    Face.N[i]=NTemp[i+iN-1]
    Face.E[i]=ETemp[i+iN-1]
    Face.P[i]=PTemp[i+iN-1]
  end
end
if length(Face.N) == 4
  OrientL=zeros(4,1)
  OrientL[1]=-1
  OrientL[2]=1
  OrientL[3]=1
  OrientL[4]=-1
  Face.OrientE = Vector{Int}(undef, 0)
  for i=1:4
    if Edges[Face.E[i]].N[1]==Face.N[i]
      push!(Face.OrientE, 1)
    else
      push!(Face.OrientE, -1)
    end
  end
end
end

