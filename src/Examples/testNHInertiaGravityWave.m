function testNHInertiaGravityWave()
close all

% Grid
nx=60;
ny=2;
nz=10;
Param.Lx=300*1.e3;
Param.Ly=2*1.e3;
x0=0;
y0=0;
Param.H=10*1.e3;

Boundary.WE='Period';
Boundary.BT='Period';
Param.hS='';
Param.Grid=CartGrid(nx,ny,Param.Lx,Param.Ly,x0,y0,@OrientFaceCart,Boundary,Param);

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
Param.ProfRho='InertiaGravityWave';
Param.ProfTheta='InertiaGravityWave';
Param.ProfVel='Const';
Param.Damping=false;
Param.Coriolis=false;
Param.Buoyancy=true;
Param.HyperD=0; %1.e15*(nx/30)^3.2;
Param.xC=Param.Lx/3;
Param.a=5*1.e3;
Param.Th0=300;
Param.DeltaTh=1.e-2;
Param.uMax=0;
Param.vMax=0;
Param.NBr=1.e-2;
Param.Equation='Compressible';
Param.TopoS='';
Param.NumV=5;
Param.RhoPos=1;
Param.uPos=2;
Param.vPos=3;
Param.wPos=4;
Param.ThPos=5;
Param.Thermo='Energy';


% Physical parameters
Param.Grav=9.81; 
Param.Coriolis=false;
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
Param.SliceXZ.Type='XZ';
Param.SliceXZ.y=y0+Param.Ly/ny/2;
Param.SliceXZ.iy=1;

% Discretization

OrdPoly=3;
[CG,Param]=Discretization(OrdPoly,Param);

% Initial conditions

U=zeros(CG.NumG,nz,Param.NumV);
U(:,:,Param.RhoPos)=Project(@fRho,CG,Param);
[U(:,:,Param.uPos),U(:,:,Param.vPos)]=ProjectVec(@fVel,CG,Param);
U(:,:,Param.ThPos)=Project(@fTheta,CG,Param).*U(:,:,Param.RhoPos);
if strcmp(Param.Thermo,'Energy')
  U(:,:,Param.ThPos)=PotToEnergy(U,CG,Param);
end
ThB=Project(@fThetaBGrd,CG,Param);
Param.fig=PlotCG(U(:,:,Param.ThPos)./U(:,:,Param.RhoPos)-ThB...
  ,CG,@TransCart,@Topo,Param,Param.fig,Param.SliceXZ);


% Integration
CFL=0.125;
dtau=.40;
time=0;

%IntMethod='Rosenbrock'; %'RungeKutta';
IntMethod='RungeKutta';
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
      if mod(i,2000)==0
        Param.fig=PlotCG(U(:,:,Param.ThPos)./U(:,:,Param.RhoPos)-ThB...
          ,CG,@TransCart,@Topo,Param,Param.fig,Param.SliceXZ);
%         Param.fig=PlotCG(U(:,:,Param.uPos),CG,@TransCart,@Topo,Param,Param.fig,Param.SliceXZ);
%         Param.fig=PlotCG(U(:,:,Param.vPos),CG,@TransCart,@Topo,Param,Param.fig,Param.SliceXZ);
%         W=zeros(size(U,1),nz+1);
%         W(:,2:nz+1)=U(:,:,Param.wPos);
%         Param.fig=PlotCG(0.5*(W(:,1:nz)+W(:,2:nz+1)),CG,@TransCart...
%           ,@Topo,Param,Param.fig,Param.SliceXZ);
      end
    end
  case 'RungeKutta'
    for i=1:nIter
      U=RungeKuttaExplicit(U,dtau,@FcnNHCurlVec,CG,Param);
      
      time=time+dtau;
      if mod(i,1000)==0
        if strcmp(Param.Thermo,'Energy')
          Th=EnergyToPot(U,CG,Param);
        else
          Th=U(:,:,Param.ThPos);
        end
        Param.fig=PlotCG(Th./U(:,:,Param.RhoPos)-ThB...
          ,CG,@TransCart,@Topo,Param,Param.fig,Param.SliceXZ);
        Param.fig=PlotCG(U(:,:,Param.uPos),CG,@TransCart,@Topo,Param,Param.fig,Param.SliceXZ);
        W=zeros(size(U,1),nz+1);
        W(:,2:nz+1)=U(:,:,Param.wPos);
        Param.fig=PlotCG(0.5*(W(:,1:nz)+W(:,2:nz+1)),CG,@TransCart...
          ,@Topo,Param,Param.fig,Param.SliceXZ);
%         frame=getframe(gcf);
%         writeVideo (v, frame);
      end
    end
end
%fig=PlotCG(U(:,1),CG,@JacobiSphere2,Param,fig);
%close(v)
end


