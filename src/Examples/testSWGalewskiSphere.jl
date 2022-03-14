function testSWGalewskiSphere
clear all
close all

% Physical parameters
Param=PhysParameters();

% Grid
nz=1;
Param.nPanel=32;
Param.H=2;
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
Param.Damping=false;
Param.Coriolis=true;
Param.CoriolisType='Sphere';
Param.Buoyancy=false;
Param.Source=false;
Param.Equation='Shallow';
Param.TopoS='';
Param.ProfVel='Galewsky';
Param.ProfRho='Galewsky';
Param.ProfTheta='Galewsky';
Param.NumV=5;
Param.RhoPos=1;
Param.uPos=2;
Param.vPos=3;
Param.wPos=4;
Param.ThPos=5;
Param.Thermo='';%'Energy'

% Discretization
OrdPoly=4;
OrdPolyZ=1;
[CG,Param]=Discretization(OrdPoly,OrdPolyZ,@JacobiSphere3,Param);
LRef=11*1.e5;
dx=2*pi*Param.RadEarth/4/Param.nPanel/OrdPoly;
Param.HyperVisc=true;
Param.HyperDCurl=2.e14;
Param.HyperDGrad=2.e14;
Param.HyperDDiv=0; % Scalars
Param.Upwind=false;

% Output
Param.Flat=true;
Param.level=1;
Param.fig=1;
Param.vtk=0;
Param.SliceXY.Type='XY';
Param.SliceXY.iz=1;
Param.vtkFileName='Galewsky';
Param.RadPrint=Param.H;
vtkGrid=vtkCGGrid(CG,@TransSphere,@Topo,Param);
NumOut=3;
Param.cNames(1).s='Vort';
Param.cNames(2).s='u';
Param.cNames(3).s='Height';
cOut=zeros(CG.NumG,nz,NumOut);


% Initial conditions

U=zeros(CG.NumG,nz,Param.NumV);
U(:,:,Param.RhoPos)=Project(@fRho,CG,Param);
[U(:,:,Param.uPos),U(:,:,Param.vPos)]=ProjectVec(@fVel,CG,Param);
U(:,:,Param.ThPos)=Project(@fTheta,CG,Param).*U(:,:,Param.RhoPos);

% Integration
time=0;
IntMethod='RungeKutta';
dtau=100;
%v = VideoWriter ('Galewsky.avi');
%open (v);
Param.RK=RungeKuttaMethod('RK3');
SimDays=6;
PrintDay=.5;
nIter=24*3600*SimDays/dtau;
PrintInt=24*3600*PrintDay/dtau;
% Print initial conditions

cOut(:,:,1)=FVort2VecDSS(U(:,:,Param.uPos),U(:,:,Param.vPos),CG,Param);
cOut(:,:,2)=U(:,:,Param.uPos);
cOut(:,:,3)=U(:,:,Param.RhoPos);

Param.vtk=vtkCG(cOut,CG,Param,vtkGrid,Param.vtk);

%
switch IntMethod
  case 'RungeKutta'
    tic
    for i=1:nIter
      i
      U=RungeKuttaExplicit(U,dtau,@FcnNHCurlVec,CG,Param);
      time=time+dtau;
      if mod(i,PrintInt)==0
        cOut(:,:,1)=FVort2VecDSS(U(:,:,Param.uPos),U(:,:,Param.vPos),CG,Param);
        cOut(:,:,2)=U(:,:,Param.uPos);
        cOut(:,:,3)=U(:,:,Param.RhoPos);
        Param.vtk=vtkCG(cOut,CG,Param,vtkGrid,Param.vtk);
      end
    end
    toc
end
%close(v)
end


