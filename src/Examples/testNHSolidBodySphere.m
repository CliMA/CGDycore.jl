function testNHSolidBodySphere
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
Param.Relax=1.e-3;
Param.Damping=false;
Param.Coriolis=true;
Param.CoriolisType='Sphere';
Param.Buoyancy=true;
Param.Source=false;
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
Param.ProfVel='SolidBody';
Param.ProfRho='SolidBody';
Param.ProfTheta='SolidBody';
Param.NumV=5;
Param.RhoPos=1;
Param.uPos=2;
Param.vPos=3;
Param.wPos=4;
Param.ThPos=5;
Param.Thermo='';%'Energy';



% Output
Param.Flat=true;
Param.level=3;
Param.fig=1;
Param.SliceXY.Type='XY';
Param.SliceXY.iz=1;


% Discretization

OrdPoly=4;
[CG,Param]=Discretization(OrdPoly,@JacobiSphere3,Param);
LRef=11*1.e5;
dx=2*pi*Param.RadEarth/4/Param.nPanel/OrdPoly;
Param.HyperD=2.e17; %1.e14*(dx/LRef)^3.2;


% Initial conditions

U=zeros(CG.NumG,nz,Param.NumV);
U(:,:,Param.RhoPos)=Project(@fRho,CG,Param);
[U(:,:,Param.uPos),U(:,:,Param.vPos)]=ProjectVec(@fVel,CG,Param);
U(:,:,Param.ThPos)=Project(@fTheta,CG,Param).*U(:,:,Param.RhoPos);
PresStart=Pressure(U(:,:,Param.ThPos),U(:,:,Param.ThPos),U(:,:,Param.ThPos),Param);
Param.fig=PlotCG(U(:,:,Param.ThPos)./U(:,:,Param.RhoPos)...
  ,CG,@TransSphere,@Topo,Param,Param.fig,Param.SliceXY);
if strcmp(Param.Thermo,'Energy')
  U(:,:,Param.ThPos)=PotToEnergy(U,CG,Param);
end


% Integration
CFL=0.125;
dtau=500;
time=0;

IntMethod='Rosenbrock'; %'RungeKutta';
%IntMethod='RungeKutta';
nIter=20000;
%v = VideoWriter ('Galewsky.avi');
%open (v);
Param.RK=RungeKuttaMethod('RK4');
Param.ROS=RosenbrockMethod('ROSRK3');
SimDays=100;
PrintDay=5;
nIter=round(24*3600*SimDays/dtau);
PrintInt=round(24*3600*PrintDay/dtau);
switch IntMethod
  case 'Rosenbrock'
    for i=1:nIter
      i
      sum(sum(U(:,:,Param.wPos)))
      U=Rosenbrock(U,dtau,@FcnNHCurlVec,@Jac,CG,Param);
      time=time+dtau;
      if mod(i,4000)==0
        Pres=Pressure(U(:,:,Param.ThPos),U(:,:,Param.ThPos),U(:,:,Param.ThPos),Param);
        if strcmp(Param.Thermo,'Energy')
          Th=EnergyToPot(U,CG,Param);
        else
          Th=U(:,:,Param.ThPos);
        end
        Param.fig=PlotCG(Th./U(:,:,Param.RhoPos)...
          ,CG,@TransSphere,@Topo,Param,Param.fig,Param.SliceXY);
        Param.fig=PlotCG(U(:,:,Param.uPos),CG,@TransSphere,@Topo,Param,Param.fig,Param.SliceXY);
        Param.fig=PlotCG(U(:,:,Param.vPos),CG,@TransSphere,@Topo,Param,Param.fig,Param.SliceXY);
        W=zeros(size(U,1),nz+1);
        W(:,2:nz+1)=U(:,:,Param.wPos);
        Param.fig=PlotCG(0.5*(W(:,1:nz)+W(:,2:nz+1)),CG,@TransSphere...
          ,@Topo,Param,Param.fig,Param.SliceXY);
      end
    end
  case 'RungeKutta'
    for i=1:nIter
      U=RungeKuttaExplicit(U,dtau,@FcnNHCurlVec,CG,Param);
      
      time=time+dtau;
      if mod(i,1)==0
        Param.fig=PlotCG(U(:,:,Param.ThPos)./U(:,:,Param.RhoPos)...
          ,CG,@TransSphere,@Topo,Param,Param.fig,Param.SliceXY);
        Param.fig=PlotCG(U(:,:,Param.uPos),CG,@TransSphere,@Topo,Param,Param.fig,Param.SliceXZ);
        Param.fig=PlotCG(U(:,:,Param.vPos),CG,@TransSphere,@Topo,Param,Param.fig,Param.SliceXZ);
        W=zeros(size(U,1),nz+1);
        W(:,2:nz+1)=U(:,:,Param.wPos);
        Param.fig=PlotCG(0.5*(W(:,1:nz)+W(:,2:nz+1)),CG,@TransSpheres...
          ,@Topo,Param,Param.fig,Param.SliceXZ);
        frame=getframe(gcf);
        writeVideo (v, frame);
      end
    end
end
%fig=PlotCG(U(:,1),CG,@JacobiSphere2,Param,fig);
%close(v)
end


