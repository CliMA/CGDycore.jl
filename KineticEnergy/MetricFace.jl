  J=zeros(OrdPolyX+1,OrdPolyY+1,Param.Grid.NumFaces);
  X=zeros(OrdPolyX+1,OrdPolyY+1,3,Param.Grid.NumFaces);
  dXdx=zeros(OrdPolyX+1,OrdPolyY+1,2,2,Param.Grid.NumFaces);
  if Param.Coriolis
    lat=zeros(OrdPolyX+1,OrdPolyY+1,Param.Grid.NumFaces);
    for iF=1:Param.Grid.NumFaces
      for j=1:OrdPolyY+1
        for i=1:OrdPolyX+1
          [X(i,j,:,iF),J(i,j,iF),dXdx(i,j,:,:,iF),lat(i,j,iF)]...
            =JacobiTrans(xwX(i),xwY(j),Param.Grid.Faces(iF),Param.Grid);
        end
      end
    end
  else
    if Param.JacobiDG
      for iF=1:Param.Grid.NumFaces
        [X(:,:,:,iF),J(:,:,iF),dXdx(:,:,:,:,iF)]...
          =JacobiDG1(FElem,Param.Grid.Faces(iF),Topo,Param);
      end
    else
      for iF=1:Param.Grid.NumFaces
        for j=1:OrdPolyY+1
          for i=1:OrdPolyX+1
            [X(i,j,:,iF),J(i,j,iF),dXdx(i,j,:,:,iF)]...
              =JacobiTrans(xwX(i),xwY(j),Param.Grid.Faces(iF),Param.Grid);
          end
        end
      end
    end
  end
  
  %
  N(Param.Grid.NumEdges).N(1,1)=1;
  T1(Param.Grid.NumEdges).T(1,1)=1;
  VolSurf(Param.Grid.NumEdges).VolSurf=zeros(OrdPolyX+1,1);
  NN(Param.Grid.NumEdges).N(1,1)=1;
  TT1(Param.Grid.NumEdges).T(1,1)=1;
  VVolSurf(Param.Grid.NumEdges).VolSurf=zeros(OrdPolyX+1,1);
  for iE=1:Param.Grid.NumEdges
    iF=Param.Grid.Edges(iE).F(1);
    FE=Param.Grid.Edges(iE).FE(1);
    if FE==1
      N(iE).N=zeros(OrdPolyX+1,2);
      T1(iE).T=zeros(OrdPolyX+1,2);
      VolSurf(iE).VolSurf=zeros(OrdPolyX+1,1);
      T1(iE).T(:,:)=dXdx(:,1,:,1,iF)...
        ./sqrt(dXdx(:,1,1,1,iF).^2+dXdx(:,1,2,1,iF).^2);
      VolSurf(iE).VolSurf(:,1)=sqrt(dXdx(:,1,1,1,iF).^2 ...
                                       +dXdx(:,1,2,1,iF).^2);
    elseif FE==2
      N(iE).N=zeros(OrdPolyY+1,2);
      T1(iE).T=zeros(OrdPolyY+1,2);
      VolSurf(iE).VolSurf=zeros(OrdPolyY+1,1);
      T1(iE).T(:,:)=dXdx(OrdPolyX+1,:,:,2,iF)...
        ./sqrt(dXdx(OrdPolyX+1,:,1,2,iF).^2+dXdx(OrdPolyX+1,:,2,2,iF).^2);
      VolSurf(iE).VolSurf(:,1)=sqrt(dXdx(OrdPolyX+1,:,1,2,iF).^2 ...
                                       +dXdx(OrdPolyX+1,:,2,2,iF).^2);
  elseif FE==3
      N(iE).N=zeros(OrdPolyX+1,2);
      T1(iE).T=zeros(OrdPolyX+1,2);
      VolSurf(iE).VolSurf=zeros(OrdPolyX+1,1);
      T1(iE).T(:,:)=dXdx(:,OrdPolyY+1,:,1,iF)...
      ./sqrt(dXdx(:,OrdPolyY+1,1,1,iF).^2+dXdx(:,OrdPolyY+1,2,1,iF).^2);
      VolSurf(iE).VolSurf(:,1)=sqrt(dXdx(:,OrdPolyY+1,1,1,iF).^2 ...
                                       +dXdx(:,OrdPolyY+1,2,1,iF).^2);     
    elseif FE==4
      N(iE).N=zeros(OrdPolyY+1,2);
      T1(iE).T=zeros(OrdPolyY+1,2);
      VolSurf(iE).VolSurf=zeros(OrdPolyY+1,1);
      T1(iE).T(:,:)=dXdx(1,:,:,2,iF)...
        ./sqrt(dXdx(1,:,1,2,iF).^2+dXdx(1,:,2,2,iF).^2);
      VolSurf(iE).VolSurf(:,1)=sqrt(dXdx(1,:,1,2,iF).^2 ...
                                       +dXdx(1,:,2,2,iF).^2);
    end
    N(iE).N(:,1)=T1(iE).T(:,2);
    N(iE).N(:,2)=-T1(iE).T(:,1);
   
    if size(Param.Grid.Edges(iE).F,2)>1
      iF=Param.Grid.Edges(iE).F(2);
      FE=Param.Grid.Edges(iE).FE(2);
      if FE==1
        NN(iE).N=zeros(OrdPolyX+1,2);
        TT1(iE).T=zeros(OrdPolyX+1,2);
        VVolSurf(iE).VolSurf=zeros(OrdPolyX+1,1);
        TT1(iE).T(:,:)=dXdx(:,1,:,1,iF)...
          ./sqrt(dXdx(:,1,1,1,iF).^2+dXdx(:,1,2,1,iF).^2);
        VVolSurf(iE).VolSurf(:,1)=sqrt(dXdx(:,1,1,1,iF).^2 ...
          +dXdx(:,1,2,1,iF).^2);
      elseif FE==2
        NN(iE).N=zeros(OrdPolyY+1,2);
        TT1(iE).T=zeros(OrdPolyY+1,2);
        VVolSurf(iE).VolSurf=zeros(OrdPolyY+1,1);
        TT1(iE).T(:,:)=dXdx(OrdPolyX+1,:,:,2,iF)...
          ./sqrt(dXdx(OrdPolyX+1,:,1,2,iF).^2+dXdx(OrdPolyX+1,:,2,2,iF).^2);
        VVolSurf(iE).VolSurf(:,1)=sqrt(dXdx(OrdPolyX+1,:,1,2,iF).^2 ...
          +dXdx(OrdPolyX+1,:,2,2,iF).^2);
      elseif FE==3
        NN(iE).N=zeros(OrdPolyX+1,2);
        TT1(iE).T=zeros(OrdPolyX+1,2);
        VVolSurf(iE).VolSurf=zeros(OrdPolyX+1,1);
        TT1(iE).T(:,:)=dXdx(:,OrdPolyY+1,:,1,iF)...
          ./sqrt(dXdx(:,OrdPolyY+1,1,1,iF).^2+dXdx(:,OrdPolyY+1,2,1,iF).^2);
        VVolSurf(iE).VolSurf(:,1)=sqrt(dXdx(:,OrdPolyY+1,1,1,iF).^2 ...
          +dXdx(:,OrdPolyY+1,2,1,iF).^2);
      elseif FE==4
        NN(iE).N=zeros(OrdPolyY+1,2);
        TT1(iE).T=zeros(OrdPolyY+1,2);
        VVolSurf(iE).VolSurf=zeros(OrdPolyY+1,1);
        TT1(iE).T(:,:)=dXdx(1,:,:,2,iF)...
          ./sqrt(dXdx(1,:,1,2,iF).^2+dXdx(1,:,2,2,iF).^2);
        VVolSurf(iE).VolSurf(:,1)=sqrt(dXdx(1,:,1,2,iF).^2 ...
          +dXdx(1,:,2,2,iF).^2);
      end
      NN(iE).N(:,1)=TT1(iE).T(:,2);
      NN(iE).N(:,2)=-TT1(iE).T(:,1);
    end
  end
  i1(1).Ind=1:OrdPolyX+1;
  i1(2).Ind=OrdPolyX+1;
  i1(3).Ind=1:OrdPolyX+1;
  i1(4).Ind=1;
  
  i2(1).Ind=1;
  i2(2).Ind=1:OrdPolyY+1;
  i2(3).Ind=OrdPolyY+1;
  i2(4).Ind=1:OrdPolyY+1;
else
  J=zeros(OrdPolyX+1,OrdPolyY+1,Param.Grid.NumFaces);
  for iF=1:Param.Grid.NumFaces
    for j=1:OrdPolyY+1
      for i=1:OrdPolyX+1
        [X,J(i,j,iF),dXdx]...
          =JacobiTrans(xwX(i),xwY(j),Param.Grid.Faces(iF),Param.Grid);
      end
    end
  end
end

[M]=MassCG(FElem,Param);
end
