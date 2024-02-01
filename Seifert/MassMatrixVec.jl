function MassMatrixVec(backend,FTB,Fe,Grid,QuadOrd)
    QuadraturQuad = Seifert.QuadRule{FTB}(Grid.Type,backend,QuadOrd)
    Weights = QuadraturQuad.Weights
    Points = QuadraturQuad.Points
    fRef=zeros(Fe.Comp,Fe.DoF,length(Weights))
    for i = 1 : length(Weights)
        for iComp = 1 : Fe.Comp
            for iD = 1 : Fe.DoF
                fRef[iComp,iD,i] = Fe.phi[iD,iComp,1](Points[i,1])*
                Fe.phi[iD,iComp,2](Points[i,2])
            end
        end
    end
    return fRef

    #=
    [wVec,xBarVec]=GetQuadMeth(QuadOrd,'GaussLegendreQuad',Grid.Type);
    fRef=zeros(2,Fe.Dof,size(wVec,1));
    for i=1:size(wVec,1)
      fRef(:,:,i)=Fe.f(xBarVec(i,1),xBarVec(i,2));
    end
    Dof=Fe.Dof;
    Val=zeros(Dof*Dof*Grid.NumFaces,1);
    RowInd=zeros(Dof*Dof*Grid.NumFaces,1);
    ColInd=zeros(Dof*Dof*Grid.NumFaces,1);
    iS=1;
    
    
    for iF=1:Grid.NumFaces
      MLoc=zeros(Dof,Dof);
      for i=1:size(wVec,1)
        [det,DF]=Trans(xBarVec(i,:),Grid.Faces(iF),Grid);
        fLoc=DF*fRef(:,:,i);
        MLoc=MLoc+1/abs(det)*wVec(i)*(fLoc'*fLoc);
      end
      nFeGlob=size(Fe.Faces(iF).Glob,1);
      RowInd(iS:iS+nFeGlob^2-1)=repmat(Fe.Faces(iF).Glob,nFeGlob,1);
      ColInd(iS:iS+nFeGlob^2-1)=reshape(repmat(Fe.Faces(iF).Glob,1,nFeGlob)',nFeGlob^2,1);
      Val(iS:iS+nFeGlob^2-1)=reshape(MLoc,nFeGlob^2,1);
      iS=iS+nFeGlob^2;
    end
    M=sparse(RowInd(1:iS-1),ColInd(1:iS-1),Val(1:iS-1));
    =#
    end
    
