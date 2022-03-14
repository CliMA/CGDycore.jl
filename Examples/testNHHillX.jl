function testNHHillX()
close all
clear all

% Physical parameters
Param=PhysParameters();

nx=60;
ny=2;
lx=60000;
ly=2000;
x0=-30000;
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
Param.Relax=5.e-3;
Param.Damping=true;
Param.Flat=false;
Param.Coriolis=false;
Param.Buoyancy=true;
Param.hC=400;
Param.x0C=0;
Param.aC=1000;
Param.Th0=300;
Param.uMax=10;
Param.vMax=0;
Param.NBr=1.e-2;
Param.Equation='Compressible';
Param.TopoS='AgnesiCartX';
Param.EndTime=2160;


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
Param.vtk=0;
Param.vtkFileName='Hill';
Param.vtk=0;
vtkGrid=vtkCGGrid(CG,@TransCart,@Topo,Param);
Param.cNames(1).s='u';
Param.cNames(2).s='w';
cOut=zeros(CG.NumG,nz,2);


%IntMethod='Rosenbrock';
IntMethod='RungeKutta';
if strcmp(IntMethod,'Rosenbrock')
  dtau=500;
else
  dtau=.4;
end
nIter=6000;
%v = VideoWriter ('Galewsky.avi');
%open (v);
Param.RK=RungeKuttaMethod('RK4');
Param.ROS=RosenbrockMethod('ROSRK3');
CFL=0.125;
time=0;
nIter=Param.EndTime/dtau;
PrintInt=100/dtau;
cOut(:,:,1)=U(:,:,Param.uPos);
W=BoundaryWOutput(U,CG,Param);
cOut(:,1,2)=0.5*(U(:,1,Param.wPos)+W);
cOut(:,2:nz-1,2)=0.5*(U(:,1:nz-2,Param.wPos)+U(:,2:nz-1,Param.wPos));
cOut(:,nz,2)=0.5*U(:,nz-1,Param.wPos);
Param.vtk=vtkCG(cOut,CG,Param,vtkGrid,Param.vtk);
switch IntMethod
  case 'RungeKutta'
    for i=1:nIter
      i
      U=RungeKuttaExplicit(U,dtau,@FcnNHCurlVec,CG,Param);
      time=time+dtau;
      if mod(i,PrintInt)==0
        cOut(:,:,1)=U(:,:,Param.uPos);
        W=BoundaryWOutput(U,CG,Param);
        cOut(:,1,2)=0.5*(U(:,1,Param.wPos)+W);
        cOut(:,2:nz-1,2)=0.5*(U(:,1:nz-2,Param.wPos)+U(:,2:nz-1,Param.wPos));
        cOut(:,nz,2)=0.5*U(:,nz-1,Param.wPos);
        Param.vtk=vtkCG(cOut,CG,Param,vtkGrid,Param.vtk);
      end
    end
end
%fig=PlotCG(U(:,1),CG,@JacobiSphere2,Param,fig);
%close(v)
end


