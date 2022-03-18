#function testNHHill_L()

using CGDycore

# Physical parameters
Param=CGDycore.PhysParameters();

nx=60;
ny=2;
lx=6000000;
ly=200000;
x0=-3000000;
y0=0;

Boundary = (;WE="Period", BT="Period")
Param.Grid=CGDycore.CartGrid(nx,ny,lx,ly,x0,y0,CGDycore.OrientFaceCart,Boundary,Param);
Param.TopoS="AgnesiCartX";

Param.H=15600;
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
Param.Coriolis=false;
Param.Thermo="";
Param.Source=false;
Param.StrideDamp=6000;
Param.Relax=1.e-1;
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
Param.ProfRho="GravityHill";
Param.ProfTheta="GravityHill";
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
Param.vtkFileName="Hill_L";
Param.vtk=0;
vtkGrid=CGDycore.vtkCGGrid(CG,CGDycore.TransCart,CGDycore.Topo,Param);
Param.cNames = [
  "Rho",
  "u",
  "v",
  "w",
  "Th"
]


IntMethod="Rosenbrock";
#IntMethod="RungeKutta";
if strcmp(IntMethod,"Rosenbrock")
  dtau=30;
else
  dtau=.4;
end
nIter=6000;
Param.RK=CGDycore.RungeKuttaMethod("RK4");
Param.ROS=CGDycore.RosenbrockMethod("SSP-Knoth");
CFL=0.125;
time=0;
Param.EndTime=216000;
nIter=Param.EndTime/dtau;
PrintTime=10000;
PrintInt=floor(PrintTime/dtau);
Param.vtk=CGDycore.vtkOutput(U,vtkGrid,CG,Param);
nIter=20

# Print initial conditions
Param.vtk=CGDycore.vtkOutput(U,vtkGrid,CG,Param);
#

str = IntMethod
if str == "Rosenbrock"
    @time begin
      for i=1:nIter
        @info "Iteration: $i"
        @show dtau
        U .= CGDycore.RosenbrockSchur(U,dtau,CGDycore.FcnNHCurlVec,CGDycore.JacSchur,CG,Param);
        if mod(i,PrintInt)==0
          Param.vtk=CGDycore.vtkOutput(U,vtkGrid,CG,Param);
        end
      end
    end

elseif str == "RungeKutta"
    for i=1:nIter
      @info "Iteration: $i"

      U .= CGDycore.RungeKuttaExplicit(U,dtau,CGDycore.FcnNHCurlVec,CG,Param);

      time[1] += dtau;
      if mod(i,PrintInt)==0
        Param.vtk=CGDycore.vtkOutput(U,vtkGrid,CG,Param);
      end
    end
else
  error("Bad str")
end




