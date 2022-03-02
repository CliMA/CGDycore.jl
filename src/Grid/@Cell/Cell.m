classdef Cell
  
  properties
    N
    E
    F
    C
    FOrient
    Mid
  end
  
  methods
    function [C,Grid]=Cell(Faces,Position,Grid,Pos)
      C.C=Pos;
      C.F=Faces;
      nF=size(Faces,2);
      nN=0;
      for i=1:nF
        if strcmp(Position(i),'R')
          Grid.Faces(Faces(i)).C(1)=Pos;
        else
          Grid.Faces(Faces(i)).C(2)=Pos;
        end
        nN=nN+size(Grid.Faces(C.F(i)).N,2);
      end
      N=zeros(1,nN);
      iN=1;
      for i=1:nF
        N(iN:iN+size(Grid.Faces(C.F(i)).N,2)-1)=Grid.Faces(C.F(i)).N;
        iN=iN+size(Grid.Faces(C.F(i)).N,2);
      end
      C.N=unique(N);
      C.Mid=[0 0 0];
      for i=1:size(C.N,2)
        C.Mid=C.Mid+Grid.Nodes(C.N(i)).P;
      end
      C.Mid=C.Mid/size(C.N,2);
      for i=1:nF
        if strcmp(Position(i),'R')
          Grid.Faces(Faces(i)).C(1)=Pos;
          if Grid.Faces(Faces(i)).n'*(Grid.Faces(Faces(i)).Mid-C.Mid)<0
            Grid.Faces(Faces(i)).n=-Grid.Faces(Faces(i)).n;
          end
        else
          Grid.Faces(Faces(i)).C(2)=Pos;
        end
        nN=nN+size(Grid.Faces(C.F(i)).N,2);
      end
      C.FOrient=ones(1,nF);
      for i=1:nF
        if Grid.Faces(C.F(i)).n'*(Grid.Faces(C.F(i)).Mid-C.Mid)'>0
          C.FOrient(i)=1;
        else
          C.FOrient(i)=-1;
        end
      end
    end
  end
end

