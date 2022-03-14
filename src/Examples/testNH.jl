close all
n=10;
Param.nPanel=n;
Param.Grav=9.81; 
Param.Omega=2*pi/(24*3600);
Param.RadEarth=6.37122d+6;
Param.Coriolis=true;
Param.HyperDM=0.0; %1.e15*(n/30)^3.2;
Param.HyperD=1.e15*(n/30)^3.2;
Param.Flat=true;


fig=1;
level=1;
OrdPoly=4;
CG.OrdPoly=OrdPoly;
Param.Grid=CubedGrid(n,@OrientFaceSphere,Param);
nz=2;
z=zeros(nz,1);
dz=500;
z(1)=dz/2;
for i=2:nz
  z(i)=z(i-1)+dz;
end
Param.Grid.nz=nz;
Param.Grid.z=z;
Param.Grid.dz=dz;
[CG.Faces,CG.NumG,CG.NumI,CG.Glob,CG.FaceGlob,CG.Stencil]...
  =NumberingFemCG(Param.Grid,OrdPoly);
Param.J=zeros(OrdPoly+1,OrdPoly+1,Param.Grid.NumFaces,nz);
Param.X=zeros(OrdPoly+1,OrdPoly+1,3,Param.Grid.NumFaces,nz);
Param.dXdx=zeros(OrdPoly+1,OrdPoly+1,Param.Grid.NumFaces,nz,3,3);
Param.lat=zeros(OrdPoly+1,OrdPoly+1,Param.Grid.NumFaces);
[CG.w,CG.xw]=GaussLobattoQuad(CG.OrdPoly);
[CG.wX,CG.xwX]=GaussLobattoQuad(CG.OrdPoly);
[CG.wY,CG.xwY]=GaussLobattoQuad(CG.OrdPoly);
[CG.DW,CG.DS]=DerivativeMatrixSingle(CG);

for iF=1:Param.Grid.NumFaces
  for iz=1:nz
  for j=1:OrdPoly+1
    for i=1:OrdPoly+1
      [Param.X(i,j,:,iF),Param.J(i,j,iF,iz),Param.dXdx(i,j,iF,iz,:,:),Param.lat(i,j,iF)]...
        =JacobiSphere3(CG.xw(i),CG.xw(j),Param.Grid.Faces(iF),z(iz),Param.Grid);
    end
  end
  end
end
[CG.M]=MassCG(CG,Param);
U=zeros(CG.NumG,nz,2);
W=zeros(CG.NumG,nz+1,1);
Profh='Galewsky';
ProfVel='Galewsky';
Param.ProfVel=ProfVel;
Param.Profh=Profh;
Rho=Project(@fh,CG,@JacobiSphere3,Param);
[U(:,:,2),U(:,:,3)]=ProjectVec(@fVel,CG,@JacobiSphere3,Param);
Param.hS=zeros(size(U(:,:,1)));
Th=Rho;
CFL=0.125;
dtau=100;
time=0;

IntMethod='RungeKutta';
nIter=200;
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
      F=FcnShallowCurlVec(UNeu,CG,Param);
      UNeu=U+1/2*dtau*F;
      F=FcnShallowCurlVec(UNeu,CG,Param);
      U=U+dtau*F;
      time=time+dtau;
      if mod(i,200)==0
        Vort=FCurl2(U(:,:,2),U(:,:,3),CG,Param);
        PlotCG(Vort,CG,@JacobiSphere3,Param,fig,level);
%         frame=getframe(gcf);
%         writeVideo (v, frame);
      end
    end
end
%fig=PlotCG(U(:,1),CG,@JacobiSphere2,Param,fig);
%close(v)

for i=1:200
GradVec=FGrad2Vec(U(:,1),CG,Param);
Grad=FGrad2(U(:,1),CG,Param);
DivVec=FDiv2Vec(U(:,2),U(:,3),CG,Param);
Div=FDiv2(U(:,2),U(:,3),CG,Param);
FVec=FCurlNon2Vec(U(:,2),U(:,3),CG,Param);
F=FCurlNon2(U(:,2),U(:,3),CG,Param);
end
bb=3;
