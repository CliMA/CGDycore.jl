function Orientation!(Edges,Faces)
  NumEdges=size(Edges,1)
  VisitE=zeros(Int, NumEdges)
  OrientationE=zeros(Int, NumEdges)
  o=1
  for iE=1:NumEdges
    (VisitE,OrientationE,o,Edges,Faces)=Orient(iE,o,VisitE,OrientationE,Edges,Faces)
  end
end

function Orient(iE,o,VisitE,OrientationE,Edges,Faces)
  if VisitE[iE]==1 || size(Edges[iE].F,2)==0
  #   if OrientationE[iE] ≠ o
  #     error("Möbius Band")
  #   end
  else
    VisitE[iE]=1
    iF1=Edges[iE].F[1]
    iFEnd=Edges[iE].F[2]
    local iE1
    for iE1_local=1:4
      iE1 = iE1_local
      if iE==Faces[iF1].E[iE1]
        break
      end
    end
    NF=[Faces[iF1].N...,Faces[iF1].N[1]]
    o=EdgeDirection(NF[iE1],NF[iE1+1],Edges[iE].N)
    OrientationE[iE]=o
    iLoop=0
    while (iF1≠iFEnd && Edges[iE].F[2] > 0) || iLoop==0
      iLoop=iLoop+1
      NF=[Faces[iF1].N...,Faces[iF1].N[1]]
      o=EdgeDirection(NF[iE1],NF[iE1+1],Edges[iE].N)
      iEP1=Shift(iE1,2)
      iEP=Faces[iF1].E[iEP1]
      if EdgeDirection(NF[iEP1],NF[iEP1+1],Edges[iEP].N)==o
        if o==-1
          Edges[iEP].N[1]=NF[iEP1]
          Edges[iEP].N[2]=NF[iEP1+1]
        else
          Edges[iEP].N[1]=NF[iEP1+1]
          Edges[iEP].N[2]=NF[iEP1]
        end
      end
      iE=iEP
      VisitE[iE]=1
      OrientationE[iE]=o
      if Edges[iE].F[2] > 0
        if iF1==Edges[iE].F[1]
          iF1=Edges[iE].F[2]
        else
          iF1=Edges[iE].F[1]
        end
      else
        iF1=Edges[iE].F[1]
      end
      for iE1_local=1:4
        iE1 = iE1_local
        if iE==Faces[iF1].E[iE1]
          break
        end
      end
    end
    if iF1==iFEnd
      NF=[Faces[iF1].N...,Faces[iF1].N[1]]
      o=EdgeDirection(NF[iE1],NF[iE1+1],Edges[iE].N)
      iEP1=Shift(iE1,2)
      iEP=Faces[iF1].E[iEP1]
      oP=EdgeDirection(NF[iEP1],NF[iEP1+1],Edges[iEP].N)
      if iEP==Edges[iEP].E
        if o==oP
          error("Möbius Band")
        end
      else
        if o==oP
          error("Möbius Band")
        end
      end
    end
  end
  return VisitE,OrientationE,o,Edges,Faces
end


function EdgeDirection(N1,N2,E)
if N1==E[1]
  Dir=1
else
  Dir=-1
end
return Dir
end

function Shift(i,S)
iS=i+S
if iS>4
  iS=iS-4
end
return iS
end

