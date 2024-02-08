function StiffMatrixFD(backend,FTB,FeF::HDivElement,FeT::ScalarElement,Grid,QuadOrd,Jacobi)
    QQ = Seifert.QuadRule{FTB}(Grid.Type,backend,QuadOrd)
    Weights = QQ.Weights
    Points = QQ.Points
    fFRef  = zeros(FeF.Comp,FeF.DoF,length(Weights))
    fTRef  = zeros(FeT.Comp,FeT.DoF,length(Weights))



    if FeF.HDiv
      for i=1:size(wVec,1)
        fTRef(:,:,i)=FeT.f(xBarVec(i,1),xBarVec(i,2))
        fFRef(:,:,i)=FeF.Divf(xBarVec(i,1),xBarVec(i,2))
      end
      DofT=FeT.Dof
      DofF=FeF.Dof
      Val=zeros(DofT*DofF*Grid.NumFaces,1)
      RowInd=zeros(DofT*DofF*Grid.NumFaces,1)
      ColInd=zeros(DofT*DofF*Grid.NumFaces,1)
      iS=1
      for iF=1:Grid.NumFaces
        MLoc=zeros(DofT,DofF)
        for i=1:size(wVec,1)
          fTLoc=fTRef(:,:,i)
          fFLoc=fFRef(:,:,i)
          MLoc=MLoc+wVec(i)*(fTLoc'*fFLoc)
        end
        for i=1:size(FeT.Faces(iF).Glob,1)
          for j=1:size(FeF.Faces(iF).Glob,1)
            Val(iS)=MLoc(i,j) [S]=StiffMatrixFD(FeT,FeF,Grid,Trans)
            [wVec,xBarVec]=GetQuadMeth(5,'GaussLobattoQuad','Quad')
            fTRef=zeros(FeT.Comp,FeT.Dof,size(wVec,1))
            fFRef=zeros(1,FeF.Dof,size(wVec,1))
            ColInd(iS)=FeF.Faces(iF).Glob(j)
            RowInd(iS)=FeT.Faces(iF).Glob(i)
            iS=iS+1
          end
        end
      end
    else
      for i=1:size(wVec,1)
        fTRef(:,:,i)=FeT.f(xBarVec(i,1),xBarVec(i,2));
        fFRef(:,:,:,i)=FeF.Gradf(xBarVec(i,1),xBarVec(i,2));
      end
      DofT=FeT.Dof;
      DofF=FeF.Dof;
      Val=zeros(DofT*DofF*Grid.NumFaces,1);
      RowInd=zeros(DofT*DofF*Grid.NumFaces,1);
      ColInd=zeros(DofT*DofF*Grid.NumFaces,1);
      iS=1;
      for iF=1:Grid.NumFaces
        MLoc=zeros(DofT,DofF);
        for i=1:size(wVec,1)
          [det,DF]=Trans([xBarVec(i,1) xBarVec(i,2)],Grid.Faces(iF),Grid);
          fTLoc=fTRef(:,:,i);
          fFLoc=fFRef(:,:,i);
          MLoc=MLoc+wVec(i)*(fTLoc'*fFLoc);
        end
        for i=1:size(FeT.Faces(iF).Glob,1)
          for j=1:size(FeF.Faces(iF).Glob,1)
            Val(iS)=MLoc(i,j);
            ColInd(iS)=FeF.Faces(iF).Glob(j);
            RowInd(iS)=FeT.Faces(iF).Glob(i);
            iS=iS+1;
          end
        end
      end
    end
    S=sparse(RowInd(1:iS-1),ColInd(1:iS-1),Val(1:iS-1));
    S=S(1:FeT.NumI,1:FeF.NumI);
    end
    