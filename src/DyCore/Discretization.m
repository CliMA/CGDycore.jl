function [CG,Param]=Discretization(OrdPoly,OrdPolyZ,Jacobi,Param) 
% Discretization

CG.OrdPoly=OrdPoly;
CG.OrdPolyZ=OrdPolyZ;
nz=Param.Grid.nz;

[CG.Faces,CG.NumG,CG.NumI,CG.Glob,CG.FaceGlob,CG.Stencil]...
  =NumberingFemCG(Param.Grid,OrdPoly);
Param.J=zeros(OrdPoly+1,OrdPoly+1,OrdPolyZ+1,Param.Grid.NumFaces,nz);
Param.lat=zeros(OrdPoly+1,OrdPoly+1,Param.Grid.NumFaces);
Param.X=zeros(OrdPoly+1,OrdPoly+1,OrdPolyZ+1,3,Param.Grid.NumFaces,nz);
Param.dXdx=zeros(OrdPoly+1,OrdPoly+1,OrdPolyZ+1,Param.Grid.NumFaces,nz,3,3);
Param.dXdxI=zeros(OrdPoly+1,OrdPoly+1,OrdPolyZ+1,Param.Grid.NumFaces,nz,3,3);
[CG.w,CG.xw]=GaussLobattoQuad(CG.OrdPoly);
[CG.wX,CG.xwX]=GaussLobattoQuad(CG.OrdPoly);
[CG.wY,CG.xwY]=GaussLobattoQuad(CG.OrdPoly);
[CG.wZ,CG.xwZ]=GaussLobattoQuad(CG.OrdPolyZ);
[CG.DW,CG.DS]=DerivativeMatrixSingle(CG.OrdPoly);
[CG.DWZ,CG.DSZ]=DerivativeMatrixSingle(CG.OrdPolyZ);


for iF=1:Param.Grid.NumFaces
  for iz=1:nz
    zI=[Param.Grid.z(iz),Param.Grid.z(iz+1)];
    [X,J,dXdx,dXdxI,lat]=Jacobi(CG,Param.Grid.Faces(iF),zI,@Topo,Param);
    Param.X(:,:,:,:,iF,iz)=X;
    Param.J(:,:,:,iF,iz)=J;
    Param.dXdx(:,:,:,iF,iz,:,:)=reshape(dXdx,OrdPoly+1,OrdPoly+1,OrdPolyZ+1,1,3,3);
    Param.dXdxI(:,:,:,iF,iz,:,:)=reshape(dXdxI,OrdPoly+1,OrdPoly+1,OrdPolyZ+1,1,3,3);
    Param.lat(:,:,iF)=lat;
  end
end
Param.JC=Average(Param.J);
Param.JF=AverageFB(Param.J);
Param.dXdxIC=zeros(OrdPoly+1,OrdPoly+1,Param.Grid.NumFaces,nz,3,3);
Param.dXdxIF=zeros(OrdPoly+1,OrdPoly+1,Param.Grid.NumFaces,nz+1,3);
for i=1:3
  for j=1:3
    Param.dXdxIC(:,:,:,:,i,j)=Average(Param.dXdxI(:,:,:,:,:,i,j));
    Param.dXdxIF(:,:,:,:,i,j)=AverageFB(Param.dXdxI(:,:,:,:,:,i,j));
  end
end

[CG.M,CG.MW]=MassCG(CG,Param);
Param.latN=zeros(CG.NumG,1);
OP=CG.OrdPoly+1;
for iM=1:size(CG.FaceGlob,2)
  Param.latN(reshape(CG.Glob(:,CG.FaceGlob(iM).Ind,:,:)...
    ,OP*OP*size(CG.FaceGlob(iM).Ind,1),1),:,:)...
    =Param.latN(reshape(CG.Glob(:,CG.FaceGlob(iM).Ind,:,:)...
    ,OP*OP*size(CG.FaceGlob(iM).Ind,1),1),:,:)...
    +reshape(Param.lat(:,:,CG.FaceGlob(iM).Ind)...
    .*Param.JC(:,:,CG.FaceGlob(iM).Ind,1)...
    ,OP*OP*size(CG.FaceGlob(iM).Ind,1),1);
end
Param.latN=Param.latN./CG.M(:,1);
end
