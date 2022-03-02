function testNHHill_L()
close all
clear all

% Physical parameters
Param=PhysParameters();

nx=60;
ny=2;
lx=6000000;
ly=200000;
x0=-3000000;
y0=0;
Boundary.WE='Period';
Boundary.BT='Period';
Param.hS='';
Param.Grid=CartGrid(nx,ny,lx,ly,x0,y0,@OrientFaceCart,Boundary,Param);


Param.H=15600;
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
Param.Relax=1.e-3; %1.e-1;
Param.Damping=true;
Param.Flat=false;
Param.Coriolis=false;
Param.Buoyancy=true;
Param.hC=400;
Param.x0C=0;
Param.aC=100000;
Param.Th0=300;
Param.uMax=10;
Param.vMax=0;
Param.NBr=1.e-2;
Param.Equation='Compressible';
Param.TopoS='AgnesiCartX';



% Discretization
OrdPoly=4;
OrdPolyZ=1;
[CG,Param]=Discretization(OrdPoly,OrdPolyZ,@JacobiDG3,Param);
LRef=11*1.e5;
Param.HyperVisc=true;
Param.HyperDCurl=2.e4; %1.e14*(dx/LRef)^3.2;
Param.HyperDGrad=2.e4;
Param.HyperDDiv=2.e4;
Param.Upwind=false;

% Initial conditions 
Param.NumV=5;
U=zeros(CG.NumG,nz,Param.NumV);
Param.ProfRho='GravityHill';
Param.ProfTheta='GravityHill';
Param.ProfVel='Const';
Param.RefProfile=false;
Param.RhoPos=1;
Param.uPos=2;
Param.vPos=3;
Param.wPos=4;
Param.ThPos=5;
U(:,:,Param.RhoPos)=Project(@fRho,CG,Param);
[U(:,:,Param.uPos),U(:,:,Param.vPos)]=ProjectVec(@fVel,CG,Param);
U(:,:,Param.ThPos)=Project(@fTheta,CG,Param).*U(:,:,Param.RhoPos);


% Output
Param.Flat=false;
Param.vtkFileName='Hill_L';
Param.vtk=0;
vtkGrid=vtkCGGrid(CG,@TransCart,@Topo,Param);
Param.cNames(1).s='u';
Param.cNames(2).s='w';
Param.cNames(3).s='Th';
Param.cNames(4).s='Rho';



IntMethod='Rosenbrock';
%IntMethod='RungeKutta';
if strcmp(IntMethod,'Rosenbrock')
  dtau=20;
else
  dtau=.4;
end
nIter=6000;
%v = VideoWriter ('Galewsky.avi');
%open (v);
Param.RK=RungeKuttaMethod('RK4');
Param.ROSRK3=RosenbrockMethod('ROSRK3');
Param.ROS=RosenbrockMethod('TROSWLASSP3P4S2C');
%Param.ROS=RosenbrockMethod('RODAS_N');
Param.ROS=RosenbrockMethod('RK3_H');
CFL=0.125;
time=0;
Param.EndTime=216000;
nIter=Param.EndTime/dtau;
PrintTime=10000;
PrintInt=floor(PrintTime/dtau);
PrintInt=100;
Param.vtk=vtkOutput(U,vtkGrid,CG,Param);
nIter=1000;
switch IntMethod
  case 'Rosenbrock'
    tic
    for i=1:nIter
      i
      U=RosenbrockSchur(U,dtau,@FcnNHCurlVec,@JacSchur,CG,Param);
      time=time+dtau;
      if mod(i,PrintInt)==0
        Param.vtk=vtkOutput(U,vtkGrid,CG,Param);
      end
    end
    toc
  case 'RungeKutta'
    for i=1:nIter
      i
      U=RungeKuttaExplicit(U,dtau,@FcnNHCurlVec,CG,Param);
      time=time+dtau;
      if mod(i,PrintInt)==0
        Param.vtk=vtkOutput(U,vtkGrid,CG,Param);
      end
    end
end
%fig=PlotCG(U(:,1),CG,@JacobiSphere2,Param,fig);
%close(v)
end



