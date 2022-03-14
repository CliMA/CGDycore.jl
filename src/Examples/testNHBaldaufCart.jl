function testNHBaldaufCart()
close all

% Physical parameters
Param=PhysParameters();

% Grid
nx=60;
ny=2;
Param.Lx=300*1.e3;
Param.Ly=2*Param.Lx/nx;
x0=0;
y0=0;
Param.H=10*1.e3;
Param.OrdPoly=4;
%Horizontal grid size
dx=Param.Lx/nx/(Param.OrdPoly+1);
dz=min(dx/2,Param.H/10);
nz=ceil(Param.H/dz);
nz=40;

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
Param.ProfRho='BaldaufCart';
Param.ProfRhoBGrd='BaldaufCart';
Param.ProfpBGrd='BaldaufCart';
Param.ProfTheta='BaldaufCart';
Param.ProfVel='Const';
Param.Damping=false;
Param.Coriolis=false;
Param.Buoyancy=true;
Param.Source=false;
Param.xc=Param.Lx/3;
Param.d=5000;
Param.a=5*1.e3;
Param.T0=250;
Param.DeltaT=1.e-2;
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
Param.Thermo='';


% Discretization
OrdPolyZ=1;
[CG,Param]=Discretization(Param.OrdPoly,OrdPolyZ,@JacobiDG3,Param);
Param.HyperVisc=false;
Param.HyperDCurl=2.e5; %1.e14*(dx/LRef)^3.2;
Param.HyperDGrad=2.e5;
Param.HyperDDiv=2.2e5;
Param.Upwind=false;

% Output
Param.Flat=false;
Param.vtkFileName='BaldaufCart';
Param.vtk=0;
vtkGrid=vtkCGGrid(CG,@TransCart,@Topo,Param);
Param.cNames(1).s='Rho';
Param.cNames(2).s='u';
Param.cNames(3).s='w';
Param.cNames(4).s='Th';
Param.cNames(5).s='TPrime';
Param.cNames(6).s='Pres';




% Initial conditions

U=zeros(CG.NumG,nz,Param.NumV);
U(:,:,Param.RhoPos)=Project(@fRho,CG,Param);
[U(:,:,Param.uPos),U(:,:,Param.vPos)]=ProjectVec(@fVel,CG,Param);
U(:,:,Param.ThPos)=Project(@fTheta,CG,Param).*U(:,:,Param.RhoPos);
Param.TBGrd=Project(@fTBGrd,CG,Param);
Param.T=Project(@fT,CG,Param);
pBGrd=Project(@fpBGrd,CG,Param);
RhoBGrd=Project(@fRhoBGrd,CG,Param);
OP=Param.OrdPoly+1;
NF=Param.Grid.NumFaces;
RhoBGrd=reshape(RhoBGrd(reshape(CG.Glob,OP*OP*NF,1),:)...
  ,OP,OP,NF,nz);
Param.RhoBGrdF=0.5*(RhoBGrd(:,:,:,1:nz-1)+RhoBGrd(:,:,:,2:nz));
Param.pBGrd=reshape(pBGrd(reshape(CG.Glob,OP*OP*NF,1),:)...
  ,OP,OP,NF,nz);

% Integration
dtau=.40;
time=0;
%IntMethod='Rosenbrock'; %'RungeKutta';
IntMethod='RungeKutta';
Param.EndTime=1800;
nIter=Param.EndTime/dtau;
PrintTime=100;
PrintInt=PrintTime/dtau;
Param.RK=RungeKuttaMethod('RK4');
Param.ROS=RosenbrockMethod('ROSRK3');
Param.vtk=vtkOutput(U,vtkGrid,CG,Param);
switch IntMethod
  case 'Rosenbrock'
    for i=1:nIter
      U=Rosenbrock(U,dtau,@FcnNHCurlVec,@Jac,CG,Param);
      time=time+dtau;
      if mod(i,PrintInt)==0
        Param.vtk=vtkOutput(U,vtkGrid,CG,Param);
      end
    end
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

end


