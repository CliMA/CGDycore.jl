function testNHHeldSuarezSphere
clear all
close all

% Physical parameters
Param=PhysParameters();

% Grid
nz=10;
Param.nPanel=4;
Param.H=30000;
Param.Grid=CubedGrid(Param.nPanel,@OrientFaceSphere,Param);


Param.Grid.nz=nz;
Param.Grid.zP=zeros(nz,1);
Param.Grid.z=zeros(nz+1,1);
Param.Grid.dz=Param.H/nz;
Param.Grid.zP(1)=Param.Grid.dz/2;
for i=2:nz
  Param.Grid.zP(i)=Param.Grid.zP(i-1)+Param.Grid.dz;
end
for i=2:nz+1
  Param.Grid.z(i)=Param.Grid.z(i-1)+Param.Grid.dz;
end


% Model
Param.ModelType='Curl';
Param.Deep=false;
Param.HeightLimit=30000.0;
Param.T0E=310.0;
Param.T0P=240.0;
Param.B=2.0;
Param.K=3.0;
Param.LapseRate=0.005;
Param.U0=-0.5;
Param.PertR=1.0/6.0;
Param.Up=0.0; %1.0;
Param.PertExpR=0.1;
Param.PertLon=pi/9.0;
Param.PertLat=2.0*pi/9.0;
Param.PertZ=15000.0;

Param.StrideDamp=6000;
Param.Relax=1.e-4;
Param.Damping=false;
Param.Coriolis=true;
Param.CoriolisType='Sphere';
Param.Buoyancy=true;
Param.Source=true;
Param.Th0=300;
Param.uMax=10;
Param.vMax=0;
Param.NBr=1.e-2;
Param.DeltaT=1;
Param.ExpDist=5;
Param.T0=300;
Param.Equation='Compressible';
Param.TopoS='';
Param.lat0=0;
Param.lon0=pi/2;
Param.ProfVel='BaroWaveSphere';
Param.ProfRho='BaroWaveSphere';
Param.ProfTheta='BaroWaveSphere';
Param.NumV=5;
Param.RhoPos=1;
Param.uPos=2;
Param.vPos=3;
Param.wPos=4;
Param.ThPos=5;
Param.Thermo='';%'Energy'
%Held Suarez
Param.day=3600*24;
Param.k_a=1/(40 * Param.day);
Param.k_f=1/Param.day;
Param.k_s=1/(4*Param.day);
Param.DeltaT_y=60;
Param.DeltaTh_z=10;
Param.T_equator=315;
Param.T_min=200;
Param.sigma_b=7/10;
Param.z_D=20.0e3;

% Discretization
OrdPoly=4;
OrdPolyZ=1;
[CG,Param]=Discretization(OrdPoly,OrdPolyZ,@JacobiSphere3,Param);
LRef=11*1.e5;
dx=2*pi*Param.RadEarth/4/Param.nPanel/OrdPoly;
Param.Upwind=false;
Param.HyperVisc=true;
Param.HyperDCurl=2.e17; %1.e14*(dx/LRef)^3.2;
Param.HyperDGrad=2.e17;
Param.HyperDDiv=2.e17; % Scalars


% Output
Param.Flat=true;
Param.level=1;
Param.fig=1;
Param.vtk=1;
Param.SliceXY.Type='XY';
Param.SliceXY.iz=3;
Param.vtkFileName='HeldSuarez';
Param.RadPrint=Param.H;
Param.vtk=0;
vtkGrid=vtkCGGrid(CG,@TransSphere,@Topo,Param);
Param.cNames(1).s='Rho';
Param.cNames(2).s='u';
Param.cNames(3).s='v';
Param.cNames(4).s='w';
Param.cNames(5).s='Th';
cOut=zeros(CG.NumG,nz,5);

% Initial conditions
U=zeros(CG.NumG,nz,Param.NumV);
U(:,:,Param.RhoPos)=Project(@fRho,CG,Param);
[U(:,:,Param.uPos),U(:,:,Param.vPos)]=ProjectVec(@fVel,CG,Param);
U(:,:,Param.ThPos)=Project(@fTheta,CG,Param).*U(:,:,Param.RhoPos);
PresStart=Pressure(U(:,:,Param.ThPos),U(:,:,Param.ThPos),U(:,:,Param.ThPos),Param);
ThB=Project(@fThetaBGrd,CG,Param);
if strcmp(Param.Thermo,'Energy')
  U(:,:,Param.ThPos)=PotToEnergy(U,CG,Param);
end


% Integration
CFL=0.125;
dtau=500;
time=0;

IntMethod='Rosenbrock';
%IntMethod='RungeKutta';
if strcmp(IntMethod,'Rosenbrock')
  dtau=600;
else
  dtau=8;
end
nIter=20000;
%v = VideoWriter ('Galewsky.avi');
%open (v);
Param.RK=RungeKuttaMethod('RK4');
Param.ROS=RosenbrockMethod('ROSRK3');
SimDays=1000;
PrintDay=10;
nIter=24*3600*SimDays/dtau;
PrintInt=24*3600*PrintDay/dtau;
% Print initial conditions
cOut(:,:,1)=U(:,:,Param.RhoPos);
cOut(:,:,2)=U(:,:,Param.uPos);
cOut(:,:,3)=U(:,:,Param.vPos);
cOut(:,1,4)=0.5*U(:,1,Param.wPos);
cOut(:,2:nz-1,4)=0.5*(U(:,1:nz-2,Param.wPos)+U(:,2:nz-1,Param.wPos));
cOut(:,nz,4)=0.5*U(:,nz-1,Param.wPos);
cOut(:,:,5)=U(:,:,Param.ThPos)./U(:,:,Param.RhoPos);
Param.vtk=vtkCG(cOut,CG,Param,vtkGrid,Param.vtk);
%
switch IntMethod
  case 'Rosenbrock'
    tic
    for i=1:nIter
      i
      U=RosenbrockSchur(U,dtau,@FcnNHCurlVec,@JacSchur,CG,Param);
      time=time+dtau;
      if mod(i,PrintInt)==0
        cOut(:,:,1)=U(:,:,Param.RhoPos);
        cOut(:,:,2)=U(:,:,Param.uPos);
        cOut(:,:,3)=U(:,:,Param.vPos);
        cOut(:,1,4)=0.5*U(:,1,Param.wPos);
        cOut(:,2:nz-1,4)=0.5*(U(:,1:nz-2,Param.wPos)+U(:,2:nz-1,Param.wPos));
        cOut(:,nz,4)=0.5*U(:,nz-1,Param.wPos);
        cOut(:,:,5)=U(:,:,Param.ThPos)./U(:,:,Param.RhoPos);
        Param.vtk=vtkCG(cOut,CG,Param,vtkGrid,Param.vtk);
      end
    end
    toc
  case 'RungeKutta'
    for i=1:nIter
      i
      U=RungeKuttaExplicit(U,dtau,@FcnNHCurlVec,CG,Param);
      
      time=time+dtau;
      if mod(i,1000)==0
        Param.fig=PlotCG(U(:,:,Param.ThPos)./U(:,:,Param.RhoPos)...
          ,CG,@TransSphere,@Topo,Param,Param.fig,Param.SliceXY);
        Param.fig=PlotCG(U(:,:,Param.uPos),CG,@TransSphere,@Topo,Param,Param.fig,Param.SliceXY);
        Param.fig=PlotCG(U(:,:,Param.vPos),CG,@TransSphere,@Topo,Param,Param.fig,Param.SliceXY);
        W=zeros(size(U,1),nz+1);
        W(:,2:nz+1)=U(:,:,Param.wPos);
        Param.fig=PlotCG(0.5*(W(:,1:nz)+W(:,2:nz+1)),CG,@TransSphere...
          ,@Topo,Param,Param.fig,Param.SliceXY);
      end
    end
end
%fig=PlotCG(U(:,1),CG,@JacobiSphere2,Param,fig);
%close(v)
end


