function [DG,Param]=InitMethodDG(OrdPoly,Param,Trans)
Param.OrdPoly=OrdPoly;
DG.OrdPoly=OrdPoly;

[DG.w,DG.xw]=GaussLobattoQuad(DG.OrdPoly);
[DG.DW,DG.DS,DG.DV,DG.DVT]=DerivativeMatrixSingle(DG);
[Param.DW,Param.DS,Param.DV,Param.DVT]=DerivativeMatrixSingle(DG);

Param.i1(1).Ind=1:Param.OrdPoly+1;
Param.i1(2).Ind=Param.OrdPoly+1;
Param.i1(3).Ind=1:Param.OrdPoly+1;
Param.i1(4).Ind=1;

Param.i2(1).Ind=1;
Param.i2(2).Ind=1:Param.OrdPoly+1;
Param.i2(3).Ind=Param.OrdPoly+1;
Param.i2(4).Ind=1:Param.OrdPoly+1;

Param.V(1)=1;
Param.V(2)=1;
Param.V(3)=-1;
Param.V(4)=-1;

[Param.w,Param.xw]=GetQuadMeth(Param.OrdPoly+1,'GaussLobattoQuad','Edge');
%[Param.DW,Param.DS,Param.DV,Param.DVT]=DerivativeMatrix(Param);
[Param.DW,Param.DS]=DerivativeMatrixSingle(Param);

Param.MX=Param.w;
Param.MZ=Param.w;

Param.X=zeros(Param.OrdPoly+1,Param.OrdPoly+1,2,Param.Grid.NumFaces);
Param.dXdx=zeros(Param.OrdPoly+1,Param.OrdPoly+1,2,2,Param.Grid.NumFaces);
Param.J=zeros(Param.OrdPoly+1,Param.OrdPoly+1,Param.Grid.NumFaces);

for iF=1:Param.Grid.NumFaces
  [Param.X(:,:,:,iF),Param.J(:,:,iF),Param.dXdx(:,:,:,:,iF)]=...
    JacobiDG1(DG,Param.Grid.Faces(iF),Param.Grid,Trans,Param);
end

Param.N=zeros(Param.OrdPoly+1,2,Param.Grid.NumEdges);
Param.T1=zeros(Param.OrdPoly+1,2,Param.Grid.NumEdges);
Param.VolSurf=zeros(Param.OrdPoly+1,Param.Grid.NumEdges);
for iE=1:Param.Grid.NumEdges
  iF=Param.Grid.Edges(iE).F(1);
  FE=Param.Grid.Edges(iE).FE(1);
  if FE==1
    Param.T1(:,:,iE)=Param.dXdx(:,1,:,1,iF)...
      ./sqrt(Param.dXdx(:,1,1,1,iF).^2+Param.dXdx(:,1,2,1,iF).^2);
    Param.VolSurf(:,iE)=sqrt(Param.dXdx(:,1,1,1,iF).^2 ...
                            +Param.dXdx(:,1,2,1,iF).^2);
  elseif FE==2
    Param.T1(:,:,iE)=Param.dXdx(Param.OrdPoly+1,:,:,2,iF)...
      ./sqrt(Param.dXdx(Param.OrdPoly+1,:,1,2,iF).^2+Param.dXdx(Param.OrdPoly+1,:,2,2,iF).^2);
    Param.VolSurf(:,iE)=sqrt(Param.dXdx(Param.OrdPoly+1,:,1,2,iF).^2 ...
                            +Param.dXdx(Param.OrdPoly+1,:,2,2,iF).^2);
    
  elseif FE==3
    Param.T1(:,:,iE)=Param.dXdx(:,Param.OrdPoly+1,:,1,iF)...
      ./sqrt(Param.dXdx(:,Param.OrdPoly+1,1,1,iF).^2+Param.dXdx(:,Param.OrdPoly+1,2,1,iF).^2);
    Param.VolSurf(:,iE)=sqrt(Param.dXdx(:,Param.OrdPoly+1,1,1,iF).^2 ...
                            +Param.dXdx(:,Param.OrdPoly+1,2,1,iF).^2);
  elseif FE==4
    Param.T1(:,:,iE)=Param.dXdx(1,:,:,2,iF)...
      ./sqrt(Param.dXdx(1,:,1,2,iF).^2+Param.dXdx(1,:,2,2,iF).^2);
    Param.VolSurf(:,iE)=sqrt(Param.dXdx(1,:,1,2,iF).^2 ...
                            +Param.dXdx(1,:,2,2,iF).^2);
  end
  Param.N(:,1,iE)=Param.T1(:,2,iE);
  Param.N(:,2,iE)=-Param.T1(:,1,iE);
end
end
