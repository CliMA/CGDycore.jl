function testNHBubbleX()
close all
clear all

% Physical parameters
Param=PhysParameters();

% Grid
nx=20;
ny=2;
lx=20000;
ly=2000;
x0=-10000;
y0=0;
Boundary.WE='Period';
Boundary.BT='Period';
Param.hS='';
Param.Grid=CartGrid(nx,ny,lx,ly,x0,y0,@OrientFaceCart,Boundary,Param);
Param.TopoS='';

Param.H=10000;
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

% Model
Param.ModelType='Curl';
%Param.ModelType='Div';
Param.Coriolis=false;
Param.Thermo='';
Param.Source=false;
Param.HyperD=0; %1.e15*(nx/30)^3.2;
Param.StrideDamp=6000;
Param.Relax=5.e-1;
Param.Damping=false;
Param.Flat=false;
Param.Coriolis=false;
Param.Buoyancy=true;
Param.xC0=0;
Param.zC0=2000;
Param.rC0=2000;
Param.Th0=300;
Param.DeltaTh=2;
Param.uMax=20;
Param.vMax=0;
Param.NBr=1.e-2;
Param.Equation='Compressible';



% Discretization
OrdPoly=4;
OrdPolyZ=1;
[CG,Param]=Discretization(OrdPoly,OrdPolyZ,@JacobiDG3,Param);
LRef=11*1.e5;
Param.HyperVisc=true;
Param.HyperDCurl=1.e4; %1.e14*(dx/LRef)^3.2;
Param.HyperDGrad=1.e4;
Param.HyperDDiv=1.e4;

% Initial conditions 
Param.NumV=5;
U=zeros(CG.NumG,nz,Param.NumV);
Param.ProfRho='WarmBubble2D';
Param.ProfTheta='WarmBubble2D';
Param.ProfVel='Const';
Param.RhoPos=1;
Param.uPos=2;
Param.vPos=3;
Param.wPos=4;
Param.ThPos=5;
U(:,:,Param.RhoPos)=Project(@fRho,CG,Param);
[U(:,:,Param.uPos),U(:,:,Param.vPos)]=ProjectVec(@fVel,CG,Param);
U(:,:,Param.ThPos)=Project(@fTheta,CG,Param).*U(:,:,Param.RhoPos);
% Output
Slice.Type='XZ';
Slice.y=y0+ly/ny/2;
Slice.iy=1;
Param.SliceXZ=Slice;
Param.SliceXY.Type='XY';
Param.SliceXY.iz=1;
Param.fig=1;


%IntMethod='Rosenbrock';
IntMethod='RungeKutta';
if strcmp(IntMethod,'Rosenbrock')
  dtau=500;
else
  dtau=.4;
end
nIter=10000;
%v = VideoWriter ('Galewsky.avi');
%open (v);
Param.RK=RungeKuttaMethod('RK4');
%Param.RK=RungeKuttaMethod('RK1');
Param.ROS=RosenbrockMethod('ROSRK3');
CFL=0.125;
time=0;

switch IntMethod
  case 'RungeKutta'
    for i=1:nIter
      i
      U=RungeKuttaExplicit(U,dtau,@FcnNHCurlVec,CG,Param);
      
      if mod(i,1000)==0
        Param.fig=PlotCG(U(:,:,Param.uPos),CG,@TransCart,@Topo,Param,Param.fig,Param.SliceXZ);
        W=zeros(size(U,1),nz+1);
        W(:,2:nz+1)=U(:,:,Param.wPos);
        W(:,1)=BoundaryW(U(:,:,Param.uPos),U(:,:,Param.vPos),CG,Param);
        Param.fig=PlotCG(0.5*(W(:,1:nz)+W(:,2:nz+1)),CG,@TransCart...
          ,@Topo,Param,Param.fig,Param.SliceXZ);
      end
      time=time+dtau;
    end
end
%fig=PlotCG(U(:,1),CG,@JacobiSphere2,Param,fig);
%close(v)
end


