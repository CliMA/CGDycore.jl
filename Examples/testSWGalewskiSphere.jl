#function testSWGalewskiSphere

using CGDycore

# Physical parameters
Param=CGDycore.PhysParameters();

# Grid
nz=1;
Param.nPanel=32;
Param.H=2;
Param.Grid=CGDycore.CubedGrid(Param.nPanel,CGDycore.OrientFaceSphere,Param);


Param.Grid.nz=nz;
Param.Grid.zP=zeros(nz,1);
Param.Grid.z=zeros(nz+1,1);
Param.Grid.dz=Param.H/nz;
Param.Grid.zP[1]=Param.Grid.dz/2;
for i=2:nz
  Param.Grid.zP[i]=Param.Grid.zP[i-1]+Param.Grid.dz;
end
for i=2:nz+1
  Param.Grid.z[i]=Param.Grid.z[i-1]+Param.Grid.dz;
end


# Model
Param.ModelType="Curl";
Param.Deep=false;
Param.HeightLimit=30000.0;
Param.Damping=false;
Param.Coriolis=true;
Param.CoriolisType="Sphere";
Param.Buoyancy=false;
Param.Source=false;
Param.Equation="Shallow";
Param.TopoS="";
Param.ProfVel="Galewsky";
Param.ProfRho="Galewsky";
Param.ProfTheta="Galewsky";
Param.NumV=5;
Param.RhoPos=1;
Param.uPos=2;
Param.vPos=3;
Param.wPos=4;
Param.ThPos=5;
Param.Thermo="";#"Energy"
Param.alphaG=1.0/3.0;
Param.betaG=1.0/15.0;
Param.hH=120.0;
Param.H0G=10000.0;
Param.uM=80.0;
Param.lat0G=pi/7.0;
Param.lat1G=pi/2.0-Param.lat0G;
Param.eN=exp(-4.0/(Param.lat1G-Param.lat0G)^2.0);
Param.RefProfile=false

# Discretization
OrdPoly=4;
OrdPolyZ=1;
(CG,Param)=CGDycore.Discretization(OrdPoly,OrdPolyZ,CGDycore.JacobiSphere3,Param);
LRef=11*1.e5;
dx=2*pi*Param.RadEarth/4/Param.nPanel/OrdPoly;
Param.HyperVisc=true;
Param.HyperDCurl=2.e14;
Param.HyperDGrad=2.e14;
Param.HyperDDiv=0; # Scalars
Param.Upwind=false;

# Output
Param.Flat=true;
Param.vtk=0;
Param.vtkFileName="Galewsky";
Param.RadPrint=Param.H;
vtkGrid=CGDycore.vtkCGGrid(CG,CGDycore.TransSphere,CGDycore.Topo,Param);
Param.cNames = [
  "Rho",
  "u",
  "v",
  "Th",
  "Vort"
]

# Initial conditions
Param.NumV=5;
U=zeros(CG.NumG,nz,Param.NumV);
U[:,:,Param.RhoPos]=CGDycore.Project(CGDycore.fRho,CG,Param);
(U[:,:,Param.uPos],U[:,:,Param.vPos])=CGDycore.ProjectVec(CGDycore.fVel,CG,Param);
U[:,:,Param.ThPos]=CGDycore.Project(CGDycore.fTheta,CG,Param).*U[:,:,Param.RhoPos];

# Integration
time=0;
IntMethod="RungeKutta";
dtau=100;
Param.RK=CGDycore.RungeKuttaMethod("RK3");
SimDays=6;
PrintDay=.5;
nIter=24*3600*SimDays/dtau;
PrintInt=24*3600*PrintDay/dtau;
# Print initial conditions
Param.vtk=CGDycore.vtkOutput(U,vtkGrid,CG,Param);
#
str = IntMethod
if str == "RungeKutta"
    @time begin
      for i=1:nIter
        @info "Iteration: $i"
        U .= CGDycore.RungeKuttaExplicit(U,dtau,CGDycore.FcnNHCurlVec,CG,Param);
        if mod(i,PrintInt)==0
          Param.vtk=CGDycore.vtkOutput(U,vtkGrid,CG,Param);
        end
      end
end
end


