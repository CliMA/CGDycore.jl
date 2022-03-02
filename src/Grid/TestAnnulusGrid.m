function TestAnnulusGrid()
close all
nx=10;
ny=10;
RadI=1;
RadO=2;
Boundary.X='Period';
Param.Grid=AnnulusGrid(nx,ny,RadI,RadO,@OrientFaceCart,Boundary);
fig=1;
PlotFaceGrid(Param.Grid,4,@JacobiAnnulus,fig)
DG.OrdPoly=2;
[DG.w,DG.xw]=GaussLobattoQuad(DG.OrdPoly);
[DG.DS,DG.DW]=DerivativeMatrix(DG);

Param.X=zeros(DG.OrdPoly+1,DG.OrdPoly+1,2,Param.Grid.NumFaces);
Param.dXdx=zeros(DG.OrdPoly+1,DG.OrdPoly+1,2,2,Param.Grid.NumFaces);
Param.J=zeros(DG.OrdPoly+1,DG.OrdPoly+1,Param.Grid.NumFaces);
for iF=1:Param.Grid.NumFaces
  [Param.J(:,:,iF),Param.dXdx(:,:,:,:,iF),Param.X(:,:,:,iF)]=...
    JacobiDG(DG,Param.Grid.Faces(iF),Param.Grid,@Annulus);
end
end

