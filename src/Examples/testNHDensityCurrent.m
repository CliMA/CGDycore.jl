function testNHDensityCurrent()
close all
clear all

% Physical parameters
Param=PhysParameters();

nx=60;
ny=2;
lx=2*25600;
ly=2000;
x0=-25600;
y0=0;
Boundary.WE='Period';
Boundary.BT='Period';
Param.hS='';
Param.Grid=CartGrid(nx,ny,lx,ly,x0,y0,@OrientFaceCart,Boundary,Param);


Param.H=6400;
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
Param.Thermo='';
Param.Source=false;
Param.Damping=false;
Param.Flat=false;
Param.Coriolis=false;
Param.Buoyancy=true;
Param.xC0=0;
Param.zC0=3000;
Param.xrC0=4000;
Param.zrC0=2000;
Param.T0=300;
Param.DeltaT=-15.0;
Param.uMax=0;
Param.vMax=0;
Param.NBr=1.e-2;
Param.Equation='Compressible';
Param.TopoS='';
Param.EndTime=900;


% Discretization
OrdPoly=4;
OrdPolyZ=1;
[CG,Param]=Discretization(OrdPoly,OrdPolyZ,@JacobiDG3,Param);
LRef=11*1.e5;
Param.HyperVisc=true;
Param.HyperDCurl=2.e5; %1.e14*(dx/LRef)^3.2;
Param.HyperDGrad=2.e5;
Param.HyperDDiv=2.2e5; %2.e2;

% Initial conditions 
Param.NumV=5;
U=zeros(CG.NumG,nz,Param.NumV);
Param.ProfRho='DensityCurrent';
Param.ProfTheta='DensityCurrent';
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
Param.vtkFileName='DensityCurrent';
vtkGrid=vtkCGGrid(CG,@TransCart,@Topo,Param);
Param.cNames(1).s='Rho';
Param.cNames(2).s='u';
Param.cNames(3).s='w';
Param.cNames(4).s='Th';
cOut=zeros(CG.NumG,nz,4);

%IntMethod='Rosenbrock';
IntMethod='RungeKutta';
dtau=.2;
Param.EndTime=900;
PrintTime=30;
%v = VideoWriter ('Galewsky.avi');
%open (v);
Param.RK=RungeKuttaMethod('RK4');
CFL=0.125;
time=0;
nIter=Param.EndTime/dtau;
PrintInt=PrintTime/dtau;
cOut(:,:,1)=U(:,:,Param.RhoPos);
cOut(:,:,2)=U(:,:,Param.uPos);
cOut(:,1,3)=0.5*U(:,1,Param.wPos);
cOut(:,2:nz-1,3)=0.5*(U(:,1:nz-2,Param.wPos)+U(:,2:nz-1,Param.wPos));
cOut(:,nz,3)=0.5*U(:,nz-1,Param.wPos);
cOut(:,:,4)=U(:,:,Param.ThPos)./U(:,:,Param.RhoPos);
Param.vtk=vtkCG(cOut,CG,Param,vtkGrid,Param.vtk);
switch IntMethod
  case 'RungeKutta'
    for i=1:nIter
      i
      U=RungeKuttaExplicit(U,dtau,@FcnNHCurlVec,CG,Param);
      time=time+dtau;
      if mod(i,PrintInt)==0
        cOut(:,:,1)=U(:,:,Param.RhoPos);
        cOut(:,:,2)=U(:,:,Param.uPos);
        cOut(:,1,3)=0.5*U(:,1,Param.wPos);
        cOut(:,2:nz-1,3)=0.5*(U(:,1:nz-2,Param.wPos)+U(:,2:nz-1,Param.wPos));
        cOut(:,nz,3)=0.5*U(:,nz-1,Param.wPos);
        cOut(:,:,4)=U(:,:,Param.ThPos)./U(:,:,Param.RhoPos);
        Param.vtk=vtkCG(cOut,CG,Param,vtkGrid,Param.vtk);
      end
    end
end
%fig=PlotCG(U(:,1),CG,@JacobiSphere2,Param,fig);
%close(v)
end


