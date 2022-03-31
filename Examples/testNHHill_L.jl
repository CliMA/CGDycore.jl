#function testNHHill_L()

using CGDycore

OrdPoly = 4
OrdPolyZ=1
nx=60;
ny=2;
nz = 40
NF = nx *ny

# Cache
cache=CGDycore.Cache(OrdPoly, OrdPolyZ, nz, NF)

# Physical parameters
Param=CGDycore.PhysParameters(cache);

# Grid
Lx=6000000;
Ly=200000;
x0=-3000000;
y0=0;
Param.Lx=Lx
Param.Ly=Ly

Boundary = (;WE="Period", BT="Period")
Param.Grid=CGDycore.CartGrid(nx,ny,Lx,Ly,x0,y0,CGDycore.OrientFaceCart,Boundary,Param);
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
Param.Relax=1.e0;
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
time=[0];
Param.EndTime=216000;
nIter=Param.EndTime/dtau;
PrintTime=10000;
PrintInt=floor(PrintTime/dtau);

#

OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
nz=Param.Grid.nz;
Param.CacheF1=zeros(OP,OP,NF,nz+1);
Param.CacheF2=zeros(OP,OP,NF,nz+1);
Param.CacheF3=zeros(OP,OP,NF,nz+1);
Param.CacheF4=zeros(OP,OP,NF,nz+1);
Param.CacheF5=zeros(OP,OP,NF,nz+1);
Param.CacheF6=zeros(OP,OP,NF,nz+1);
Param.CacheC1 = view(Param.CacheF1,:,:,:,1:nz)
Param.CacheC2 = view(Param.CacheF2,:,:,:,1:nz)
Param.CacheC3 = view(Param.CacheF3,:,:,:,1:nz)
Param.CacheC4 = view(Param.CacheF4,:,:,:,1:nz)
Param.CacheC5 = view(Param.CacheF5,:,:,:,1:nz)
Param.CacheC6 = view(Param.CacheF6,:,:,:,1:nz)
Param.Cache1=zeros(CG.NumG,nz)
Param.Cache2=zeros(CG.NumG,nz)
Param.Cache3=zeros(CG.NumG,nz)
Param.Cache4=zeros(CG.NumG,nz)
Param.Pres=zeros(OP,OP,NF,nz)
Param.KE=zeros(OP,OP,NF,nz)
Param.FCG=zeros(OP,OP,NF,nz,size(U,3))
Param.Vn=zeros(size(U))
Param.RhoCG=zeros(OP,OP,NF,nz)
Param.v1CG=zeros(OP,OP,NF,nz)
Param.v2CG=zeros(OP,OP,NF,nz)
Param.wCG=zeros(OP,OP,NF,nz+1)
Param.wCCG=zeros(OP,OP,NF,nz+1)
Param.ThCG=zeros(OP,OP,NF,nz)
str = IntMethod
if str == "Rosenbrock"
  Param.J = CGDycore.JacStruct(CG.NumG,nz)
  Param.k=zeros(size(U)..., Param.ROS.nStage);
  Param.fV=zeros(size(U))
elseif str == "RungeKutta"
  Param.f=zeros(size(U)..., Param.RK.nStage);
end

# Print initial conditions
Param.vtk=CGDycore.vtkOutput(U,vtkGrid,CG,Param);

if str == "Rosenbrock"
    @time begin
      for i=1:nIter
        Δt = @elapsed begin
          CGDycore.RosenbrockSchur!(U,dtau,CGDycore.FcnNHCurlVec!,CGDycore.JacSchur,CG,Param);
          time[1] += dtau;
          if mod(i,PrintInt)==0
            Param.vtk=CGDycore.vtkOutput(U,vtkGrid,CG,Param);
          end
        end
        percent = i/nIter*100
        @info "Iteration: $i took $Δt, $percent% complete"
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




