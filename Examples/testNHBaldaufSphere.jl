function testNHBaldaufSphere()
close all
clear all

% Physical parameters
Param=PhysParameters();

% Grid
nz=10;
Param.nPanel=8;
Param.H=10000;
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
LRef=11*1.e5;
Param.HyperD=0; %1.e14*(dx/LRef)^3.2;
Param.StrideDamp=6000;
Param.Relax=5.e-1;
Param.Damping=true;
Param.Coriolis=false;
Param.CoriolisType='Sphere';
Param.Buoyancy=true;
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
Param.ProfVel='Baldauf';
Param.ProfRho='Baldauf';
Param.ProfTheta='Baldauf';
Param.NumV=5;
Param.RhoPos=1;
Param.uPos=2;
Param.vPos=3;
Param.wPos=4;
Param.ThPos=5;
Param.Thermo='Energy';



% Output
Param.Flat=true;
Param.level=1;
Param.fig=1;
Param.SliceXY.Type='XY';
Param.SliceXY.iz=4;


% Discretization

OrdPoly=4;
[CG,Param]=Discretization(OrdPoly,@JacobiSphere3,Param);

% Initial conditions

U=zeros(CG.NumG,nz,Param.NumV);
U(:,:,Param.RhoPos)=Project(@fRho,CG,Param);
[U(:,:,Param.uPos),U(:,:,Param.vPos)]=ProjectVec(@fVel,CG,Param);
U(:,:,Param.ThPos)=Project(@fTheta,CG,Param).*U(:,:,Param.RhoPos);
ThB=Project(@fThetaBGrd,CG,Param);
Param.fig=PlotCG(U(:,:,Param.ThPos)./U(:,:,Param.RhoPos)-ThB...
  ,CG,@TransSphere,@Topo,Param,Param.fig,Param.SliceXY);
if strcmp(Param.Thermo,'Energy')
  U(:,:,Param.ThPos)=PotToEnergy(U,CG,Param);
end


% Integration
CFL=0.125;
dtau=240;
time=0;

IntMethod='Rosenbrock'; %'RungeKutta';
%IntMethod='RungeKutta';
nIter=20000;
%v = VideoWriter ('Galewsky.avi');
%open (v);
Param.RK=RungeKuttaMethod('RK4');
Param.ROS=RosenbrockMethod('ROSRK3');
switch IntMethod
  case 'Rosenbrock'
    for i=1:nIter
      U=Rosenbrock(U,dtau,@FcnNHCurlVec,@Jac,CG,Param);
      time=time+dtau;
      if mod(i,10)==0
        if strcmp(Param.Thermo,'Energy')
          Th=EnergyToPot(U,CG,Param);
        else
          Th=U(:,:,Param.ThPos);
        end
        Param.fig=PlotCG(Th./U(:,:,Param.RhoPos)-ThB...
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


