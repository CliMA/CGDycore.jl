function testNHHillY()
close all
nx=2;
ny=60;
lx=2000;
ly=60000;
x0=0;
y0=-30000;
Param.Grav=9.81; 
Param.Coriolis=false;
Param.HyperD=0; %1.e15*(nx/30)^3.2;
Param.StrideDamp=6000;
Param.Relax=5.e-1;
Param.Damping=true;
Param.Flat=false;
Param.hC=400;
Param.y0C=0;
Param.aC=1000;
Param.H=15600;
Param.cS=360;
Param.Grav=9.81d0;
Param.Cpd=1004.0d0;
Param.Cvd=717.0d0;
Param.Rd=Param.Cpd-Param.Cvd;
Param.p0=1.0d5;
Param.Cpv=1885.0d0;
Param.Gamma=Param.Cpd/Param.Cvd;
Param.kappa=Param.Rd/Param.Cpd;
Param.Th0=300;
Param.uMax=0;
Param.vMax=10;
Param.NBr=1.e-2;
Param.TopoS='AgnesiCartY';
fig=1;
level=1;
OrdPoly=3;
CG.OrdPoly=OrdPoly;
Boundary.WE='Period';
Boundary.BT='Period';
Param.hS='';
Param.Grid=CartGrid(nx,ny,lx,ly,x0,y0,@OrientFaceCart,Boundary,Param);

nz=40;
zP=zeros(nz,1);
z=zeros(nz+1,1);
dz=Param.H/nz;
zP(1)=dz/2;
for i=2:nz
  zP(i)=zP(i-1)+dz;
end
for i=2:nz+1
  z(i)=z(i-1)+dz;
end

Param.Grid.nz=nz;
Param.Grid.z=z;
Param.Grid.zP=zP;
Param.Grid.dz=dz;
[CG.Faces,CG.NumG,CG.NumI,CG.Glob,CG.FaceGlob,CG.Stencil]...
  =NumberingFemCG(Param.Grid,OrdPoly);
Param.JF=zeros(OrdPoly+1,OrdPoly+1,Param.Grid.NumFaces,nz+1);
Param.XF=zeros(OrdPoly+1,OrdPoly+1,3,Param.Grid.NumFaces,nz+1);
Param.dXdxF=zeros(OrdPoly+1,OrdPoly+1,Param.Grid.NumFaces,nz+1,3,3);
Param.dXdxIF=zeros(OrdPoly+1,OrdPoly+1,Param.Grid.NumFaces,nz+1,3,3);
[CG.w,CG.xw]=GaussLobattoQuad(CG.OrdPoly);
[CG.wX,CG.xwX]=GaussLobattoQuad(CG.OrdPoly);
[CG.wY,CG.xwY]=GaussLobattoQuad(CG.OrdPoly);
[CG.DW,CG.DS]=DerivativeMatrixSingle(CG);

for iF=1:Param.Grid.NumFaces
  for iz=1:nz+1
    [X,J,dXdx,dXdxI]=JacobiDG3(CG,Param.Grid.Faces(iF),z(iz),@Topo,Param);
    Param.XF(:,:,:,iF,iz)=X;
    Param.JF(:,:,iF,iz)=J;
    Param.dXdxF(:,:,iF,iz,:,:)=reshape(dXdx,OrdPoly+1,OrdPoly+1,1,3,3);
    Param.dXdxIF(:,:,iF,iz,:,:)=reshape(dXdxI,OrdPoly+1,OrdPoly+1,1,3,3);
  end
end
Param.X=0.5*(Param.XF(:,:,:,:,1:nz)+Param.XF(:,:,:,:,2:nz+1));
Param.J=0.5*(Param.JF(:,:,:,1:nz)+Param.JF(:,:,:,2:nz+1));
Param.dXdx(:,:,:,1:nz,:,:)=0.5*(Param.dXdxF(:,:,:,1:nz,:,:)...
                               +Param.dXdxF(:,:,:,2:nz+1,:,:));
Param.dXdxI(:,:,:,1:nz,:,:)=0.5*(Param.dXdxIF(:,:,:,1:nz,:,:)...
                                +Param.dXdxIF(:,:,:,2:nz+1,:,:));                             
[CG.M,CG.MW]=MassCG(CG,Param);

U=zeros(CG.NumG,nz,2);
W=zeros(CG.NumG,nz+1,1);
ProfRho='GravityHill';
ProfTheta='GravityHill';
ProfVel='Const';
Param.ProfVel=ProfVel;
Param.ProfRho=ProfRho;
Param.ProfTheta=ProfTheta;
Rho=Project(@fRho,CG,Param);
Slice.Type='XZ';
Slice.y=y0+ly/ny/2;
Slice.iy=1;
Param.SliceXZ=Slice;
Param.SliceYZ.Type='YZ';
Param.SliceYZ.ix=1;
Param.SliceYZ.x=x0+lx/nx/2;
Param.SliceXY.Type='XY';
Param.SliceXY.iz=1;
%PlotCG(Rho,CG,@TransCart,@Topo,Param,fig,Slice);
[U(:,:,1),U(:,:,2)]=ProjectVec(@fVel,CG,Param);
Th=Project(@fTheta,CG,Param).*Rho;
CFL=0.125;
dtau=.4;
time=0;

IntMethod='RungeKutta';
nIter=20000;
%v = VideoWriter ('Galewsky.avi');
%open (v);
switch IntMethod
  case 'RungeKutta'
    for i=1:nIter
      [FRho,FU,FW,FTh]=FcnNHCurlVec(Rho,U,W,Th,CG,Param);
      RhoNeu=Rho+1/3*dtau*FRho;
      UNeu=U+1/3*dtau*FU;
      WNeu=W+1/3*dtau*FW;
      ThNeu=Th+1/3*dtau*FTh;
      [FRho,FU,FW,FTh]=FcnNHCurlVec(RhoNeu,UNeu,WNeu,ThNeu,CG,Param);
      RhoNeu=Rho+1/2*dtau*FRho;
      UNeu=U+1/2*dtau*FU;
      WNeu=W+1/2*dtau*FW;
      ThNeu=Th+1/2*dtau*FTh;
      [FRho,FU,FW,FTh]=FcnNHCurlVec(RhoNeu,UNeu,WNeu,ThNeu,CG,Param);
      Rho=Rho+dtau*FRho;
      U=U+dtau*FU;
      W=W+dtau*FW;
      Th=Th+dtau*FTh;
      time=time+dtau;
      if mod(i,1000)==0
        
        fig=PlotCG(U(:,:,1),CG,@TransCart,@Topo,Param,fig,Param.SliceYZ);
        fig=PlotCG(U(:,:,2),CG,@TransCart,@Topo,Param,fig,Param.SliceYZ);
        fig=PlotCG(0.5*(W(:,1:nz)+W(:,2:nz+1)),CG,@TransCart...
          ,@Topo,Param,fig,Param.SliceYZ);
%         frame=getframe(gcf);
%         writeVideo (v, frame);
      end
    end
end
%fig=PlotCG(U(:,1),CG,@JacobiSphere2,Param,fig);
%close(v)
end


