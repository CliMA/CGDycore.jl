function TestCart()
nx=4;
nz=4;
xL=0;
xU=10;
zL=0;
zU=10;
Boundary.BT='FreeSlip';
Boundary.WE='FreeSlip';
Param.hS='AgnesiCart';
Param.hC=1;
Param.x0C=5;
Param.aC=1;

Param.Grid=CartGrid(nx,nz,xU-xL,zU-zL,xL,zL,@OrientFaceCart,Boundary,Param);

DG.OrdPoly=2;
DG.OrdPolyX=2;
DG.OrdPolyY=2;
[DG.w,DG.xw]=GaussLobattoQuad(DG.OrdPoly);
[DG.wX,DG.xwX]=GaussLobattoQuad(DG.OrdPolyX);
[DG.wY,DG.xwY]=GaussLobattoQuad(DG.OrdPolyX);
[DG.DW,DG.DS,DG.DV,DG.DVT,DG.B]=DerivativeMatrixSingle(DG);
[DG.DWX,DG.DSX]=DerivativeMatrixSingle(DG);
[DG.DWY,DG.DSY]=DerivativeMatrixSingle(DG);

Param.X=zeros(DG.OrdPoly+1,DG.OrdPoly+1,3,Param.Grid.NumFaces);
Param.dXdx=zeros(DG.OrdPoly+1,DG.OrdPoly+1,2,2,Param.Grid.NumFaces);
Param.J=zeros(DG.OrdPoly+1,DG.OrdPoly+1,Param.Grid.NumFaces);
for iF=1:Param.Grid.NumFaces
  [Param.X(:,:,:,iF),Param.J(:,:,iF),Param.dXdx(:,:,:,:,iF)]=...
    JacobiDG1(DG,Param.Grid.Faces(iF),Param.Grid,@Cart,Param);
end
Subs=1;
fig=1;
PlotFaceGrid(Param.Grid,Subs,@JacobiCart,fig)
end

