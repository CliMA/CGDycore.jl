# function testNHBubbleX()
using CGDycore

# Physical parameters
Param=CGDycore.PhysParameters();


# Grid
nx=20;
ny=2;
lx=20000;
ly=2000;
x0=-10000;
y0=0;
Boundary = (;WE="Period", BT="Period")
Param.hS="";
Param.Grid=CGDycore.CartGrid(nx,ny,lx,ly,x0,y0,CGDycore.OrientFaceCart,Boundary,Param);
Param.TopoS="";

Param.H=10000;
nz=40;
zP=zeros(nz,1);
z=zeros(nz+1,1);
dz=Param.H/nz;
zP[1]=dz/2;
for i=2:nz
  zP[i]=zP[i-1]+dz;
end
for i=2:nz+1
  z[i]=z[i-1]+dz;
end

Param.Grid.nz=nz;
Param.Grid.z=z;
Param.Grid.zP=zP;
Param.Grid.dz=dz;

# Model
Param.ModelType="Curl";
#Param.ModelType="Div";
Param.Coriolis=false;
Param.Thermo="";
Param.Source=false;
Param.Damping=false;
Param.Flat=false;
Param.Coriolis=false;
Param.Buoyancy=true;
Param.xC0=0;
Param.zC0=2000;
Param.rC0=2000;
Param.Th0=300;
Param.DeltaTh=2;
Param.uMax=20;
Param.vMax=0;
Param.NBr=1.e-2;
Param.Equation="Compressible";
Param.RefProfile=false


# Discretization
OrdPoly=4;
OrdPolyZ=1;
(CG,Param)=CGDycore.Discretization(OrdPoly,OrdPolyZ,CGDycore.JacobiDG3,Param);
Param.HyperVisc=true;
Param.HyperDCurl=1.e4;
Param.HyperDGrad=1.e4;
Param.HyperDDiv=1.e4;
Param.Upwind=false

# Initial conditions 
Param.NumV=5;
U=zeros(CG.NumG,nz,Param.NumV);
Param.ProfRho="WarmBubble2D";
Param.ProfTheta="WarmBubble2D";
Param.ProfVel="Const";
Param.RhoPos=1;
Param.uPos=2;
Param.vPos=3;
Param.wPos=4;
Param.ThPos=5;
U[:,:,Param.RhoPos]=CGDycore.Project(CGDycore.fRho,CG,Param);
(U[:,:,Param.uPos],U[:,:,Param.vPos])=CGDycore.ProjectVec(CGDycore.fVel,CG,Param);
U[:,:,Param.ThPos]=CGDycore.Project(CGDycore.fTheta,CG,Param).*U[:,:,Param.RhoPos];

# Output
Param.vtkFileName="Bubble";
Param.vtk=0;
vtkGrid=CGDycore.vtkCGGrid(CG,CGDycore.TransCart,CGDycore.Topo,Param);
Param.cNames = [
  "Rho",
  "u",
  "v",
  "w",
  "Th"
]


time=[0];
#IntMethod="Rosenbrock";
IntMethod="RungeKutta";
if strcmp(IntMethod,"Rosenbrock")
  dtau=500;
else
  dtau=.4;
end
nIter=10000;
Param.RK=CGDycore.RungeKuttaMethod("RK4");
Param.ROS=CGDycore.RosenbrockMethod("ROSRK3");
CFL=0.125;
time=0;
SimTime=1000;
PrintTime=100;
nIter=SimTime/dtau;
PrintInt=PrintTime/dtau;
# Print initial conditions
Param.vtk=CGDycore.vtkOutput(U,vtkGrid,CG,Param);
#

str = IntMethod
if str == "RungeKutta"
    for i=1:nIter
      @info "Iteration: $i"
      U .= CGDycore.RungeKuttaExplicit(U,dtau,CGDycore.FcnNHCurlVec,CG,Param);
      if mod(i,PrintInt)==0
        Param.vtk=CGDycore.vtkOutput(U,vtkGrid,CG,Param);
      end
    end
end


