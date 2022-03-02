function testNHBaroWaveCart()
close all
clear all

% Grid
nz=20;
Param.Lx=30000*1.e3;
Param.Ly=6000*1.e3;
Param.H=30*1.e3;
dx=400*1.e3;
dy=400*1.e3;
nx=Param.Lx/dx;
ny=2*Param.Ly/dy;
x0=0;
y0=-Param.Ly;

Param.Boundary.WE='Period';
Param.Boundary.BT='Period';
Param.hS='';
Param.Grid=CartGrid(nx,ny,Param.Lx,2*Param.Ly,x0,y0,@OrientFaceCart,Param.Boundary,Param);

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
Param.Flat=false;
Param.Coriolis=true;
Param.CoriolisType='Beta-Plane';
Param.Buoyancy=true;
Param.Th0=300;
Param.uMax=10;
Param.vMax=0;
Param.NBr=1.e-2;
Param.lat0=pi/2;
Param.Equation='Compressible';
Param.TopoS='';
Param.Omega=2*pi/(24*3600);
Param.RadEarth=6.37122d+6;
Param.f0=2*sin(Param.lat0)*Param.Omega;
Param.beta0=2*cos(Param.lat0)*Param.Omega/Param.RadEarth;
Param.y0=Param.Ly/2;
Param.xC=2000*1.e3;
Param.yC=2500*1.e3;
Param.Lp=600*1.e3;
Param.uP=1;
Param.b=2;
Param.Lapse=0.005;
Param.T0=288;
Param.u0=35;
Param.ProfVel='BaroWaveCart';
Param.ProfRho='BaroWaveCart';
Param.ProfTheta='BaroWaveCart';
Param.NumV=5;
Param.RhoPos=1;
Param.uPos=2;
Param.vPos=3;
Param.wPos=4;
Param.ThPos=5;
Param.Thermo='Energy';


% Physical parameters
Param.Grav=9.81; 
Param.cS=360;
Param.Grav=9.81d0;
Param.Cpd=1004.0d0;
Param.Cvd=717.0d0;
Param.Rd=Param.Cpd-Param.Cvd;
Param.p0=1.0d5;
Param.Cpv=1885.0d0;
Param.Gamma=Param.Cpd/Param.Cvd;
Param.kappa=Param.Rd/Param.Cpd;

% Output
Param.Flat=false;
Param.level=1;
Param.fig=1;
Param.SliceYZ.Type='YZ';
Param.SliceYZ.x=x0+Param.Lx/nx/2;
Param.SliceYZ.ix=10;
Param.SliceXY.Type='XY';
Param.SliceXY.z=1500;
Param.SliceXY.iz=1;


% Discretization

OrdPoly=3;
[CG,Param]=Discretization(OrdPoly,@JacobiDG3,Param);

% Initial conditions

U=zeros(CG.NumG,nz,Param.NumV);
U(:,:,Param.RhoPos)=Project(@fRho,CG,Param);
[U(:,:,Param.uPos),U(:,:,Param.vPos)]=ProjectVec(@fVel,CG,Param);
U(:,:,Param.ThPos)=Project(@fTheta,CG,Param).*U(:,:,Param.RhoPos);
if strcmp(Param.Thermo,'Energy')
  U(:,:,Param.ThPos)=PotToEnergy(U,CG,Param);
end
ThB=Project(@fThetaBGrd,CG,Param);
Param.fig=PlotCG(U(:,:,Param.uPos)...
  ,CG,@TransCart,@Topo,Param,Param.fig,Param.SliceXY);


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
      if mod(i,360)==0
        if strcmp(Param.Thermo,'Energy')
          Th=EnergyToPot(U,CG,Param);
        else
          Th=U(:,:,Param.ThPos);
        end
        Param.fig=PlotCG(Th./U(:,:,Param.RhoPos)-ThB...
          ,CG,@TransCart,@Topo,Param,Param.fig,Param.SliceXY);
%         Param.fig=PlotCG(U(:,:,Param.uPos),CG,@TransCart,@Topo,Param,Param.fig,Param.SliceXY);
%         Param.fig=PlotCG(U(:,:,Param.vPos),CG,@TransCart,@Topo,Param,Param.fig,Param.SliceXY);
        W=zeros(size(U,1),nz+1);
        W(:,2:nz+1)=U(:,:,Param.wPos);
        Param.fig=PlotCG(0.5*(W(:,1:nz)+W(:,2:nz+1)),CG,@TransCart...
          ,@Topo,Param,Param.fig,Param.SliceXY);
      end
    end
  case 'RungeKutta'
    for i=1:nIter
      U=RungeKuttaExplicit(U,dtau,@FcnNHCurlVec,CG,Param);
      
      time=time+dtau;
      if mod(i,1)==0
        Param.fig=PlotCG(U(:,:,Param.ThPos)./U(:,:,Param.RhoPos)...
          ,CG,@TransCart,@Topo,Param,Param.fig,Param.SliceXY);
        Param.fig=PlotCG(U(:,:,Param.uPos),CG,@TransCart,@Topo,Param,Param.fig,Param.SliceXZ);
        Param.fig=PlotCG(U(:,:,Param.vPos),CG,@TransCart,@Topo,Param,Param.fig,Param.SliceXZ);
        W=zeros(size(U,1),nz+1);
        W(:,2:nz+1)=U(:,:,Param.wPos);
        Param.fig=PlotCG(0.5*(W(:,1:nz)+W(:,2:nz+1)),CG,@TransCart...
          ,@Topo,Param,Param.fig,Param.SliceXZ);
        frame=getframe(gcf);
        writeVideo (v, frame);
      end
    end
end
%fig=PlotCG(U(:,1),CG,@JacobiSphere2,Param,fig);
%close(v)
end


