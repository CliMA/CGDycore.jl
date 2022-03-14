function TestHill()

close all
nx=20;
ny=1;
lx=2000;
ly=100;
x0=-1000;
y0=0;
Param.Grav=9.81; 
Param.Coriolis=false;
Param.HyperDM=1.e15*(nx/30)^3.2;
Param.HyperD=1.e15*(nx/30)^3.2;
Param.Flat=false;

Boundary.WE='Period';
Boundary.BT='Period';
Param.Grid=CartGridCartGrid(nx,ny,lx,ly,x0,y0,OrientFaceCart,Boundary,Param);

OrdPoly=3;
CG.OrdPoly=OrdPoly;
[CG.w,CG.xw]=GaussLobattoQuad(CG.OrdPoly);
[CG.wX,CG.xwX]=GaussLobattoQuad(CG.OrdPoly);
[CG.wY,CG.xwY]=GaussLobattoQuad(CG.OrdPoly);
[CG.DW,CG.DS]=DerivativeMatrixSingle(CG);

Param.J=zeros(OrdPoly+1,OrdPoly+1,Param.Grid.NumFaces);
Param.X=zeros(OrdPoly+1,OrdPoly+1,3,Param.Grid.NumFaces);
Param.dXdx=zeros(OrdPoly+1,OrdPoly+1,2,2,Param.Grid.NumFaces);
Param.lat=zeros(OrdPoly+1,OrdPoly+1,Param.Grid.NumFaces);
% Param.J1=zeros(OrdPoly+1,OrdPoly+1,Param.Grid.NumFaces);
% Param.X1=zeros(OrdPoly+1,OrdPoly+1,3,Param.Grid.NumFaces);
% Param.dXdx1=zeros(OrdPoly+1,OrdPoly+1,2,2,Param.Grid.NumFaces);
% Param.lat1=zeros(OrdPoly+1,OrdPoly+1,Param.Grid.NumFaces);
for iF=1:Param.Grid.NumFaces
%   [Param.X1(:,:,:,iF),Param.J1(:,:,iF),Param.dXdx1(:,:,:,:,iF),Param.lat1(:,:,iF)]=...
%     JacobiSphereDG(CG,Param.Grid.Faces(iF),Param.Grid);
  for j=1:OrdPoly+1
    for i=1:OrdPoly+1
      [Param.X(i,j,:,iF),Param.J(i,j,iF),Param.dXdx(i,j,:,:,iF),Param.lat(i,j,iF)]...
        =JacobiSphere2(CG.xw(i),CG.xw(j),Param.Grid.Faces(iF),Param.Grid);
    end
  end
end

%Setting initial Values
Profh='Hill';
ProfVel='Hill';
Param.ProfVel=ProfVel;
Param.Profh=Profh;

[CG.Faces,CG.NumG,CG.NumI]=NumberingFemCG(Param.Grid,OrdPoly);
[CG.M]=MassCG(CG,Param);


U=zeros(CG.NumG,3);
U(:,1)=Project(@fh,CG,@JacobiSphere2,Param);
Param.hS=Project(@fhs,CG,@JacobiSphere2,Param);
[U(:,2),U(:,3)]=ProjectVec(@fVel,CG,@JacobiSphere2,Param);
fig=1;
fig=PlotCG(U(:,2),CG,@JacobiSphere2,Param,fig);

CFL=0.125;
dtau=100;
time=0;

IntMethod='RungeKutta';
nIter=5*3600*24/100;
v = VideoWriter ('Hill.avi');
open (v);
switch IntMethod
  case 'RungeKutta'
    for i=1:nIter
      F=FcnShallowCurl(U,CG,Param);
      UNeu=U+1/3*dtau*F;
      F=FcnShallowCurl(UNeu,CG,Param);
      UNeu=U+1/2*dtau*F;
      F=FcnShallowCurl(UNeu,CG,Param);
      U=U+dtau*F;
      time=time+dtau;
      if mod(i,50)==0
        PlotCG(U(:,2),CG,@JacobiSphere2,Param,fig);
        frame=getframe(gcf);
        writeVideo (v, frame);
      end
    end
end
%fig=PlotCG(U(:,1),CG,@JacobiSphere2,Param,fig);
close(v)
end
