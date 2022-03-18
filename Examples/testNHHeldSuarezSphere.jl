#function testNHHeldSuarezSphere
using CGDycore

# Physical parameters
Param=CGDycore.PhysParameters();

Param.Upwind = false
Param.RefProfile = false
# Grid
nz=10;
Param.nPanel=4;
Param.H=30000;
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
Param.T0E=310.0;
Param.T0P=240.0;
Param.B=2.0;
Param.K=3.0;
Param.LapseRate=0.005;
Param.U0=-0.5;
Param.PertR=1.0/6.0;
Param.Up=0.0; #1.0;
Param.PertExpR=0.1;
Param.PertLon=pi/9.0;
Param.PertLat=2.0*pi/9.0;
Param.PertZ=15000.0;

Param.StrideDamp=6000;
Param.Relax=1.e-4;
Param.Damping=false;
Param.Coriolis=true;
Param.CoriolisType="Sphere";
Param.Buoyancy=true;
Param.Source=true;
Param.Th0=300;
Param.uMax=10;
Param.vMax=0;
Param.NBr=1.e-2;
Param.DeltaT=1;
Param.ExpDist=5;
Param.T0=300;
Param.Equation="Compressible";
Param.TopoS="";
Param.lat0=0;
Param.lon0=pi/2;
Param.ProfVel="BaroWaveSphere";
Param.ProfRho="BaroWaveSphere";
Param.ProfTheta="BaroWaveSphere";
Param.NumV=5;
Param.RhoPos=1;
Param.uPos=2;
Param.vPos=3;
Param.wPos=4;
Param.ThPos=5;
Param.Thermo="";#"Energy"
#Held Suarez
Param.day=3600*24;
Param.k_a=1/(40 * Param.day);
Param.k_f=1/Param.day;
Param.k_s=1/(4*Param.day);
Param.DeltaT_y=60;
Param.DeltaTh_z=10;
Param.T_equator=315;
Param.T_min=200;
Param.sigma_b=7/10;
Param.z_D=20.0e3;

# Discretization
OrdPoly=4;
OrdPolyZ=1;
(CG,Param)=CGDycore.Discretization(OrdPoly,OrdPolyZ,CGDycore.JacobiSphere3,Param);
LRef=11*1.e5;
dx=2*pi*Param.RadEarth/4/Param.nPanel/OrdPoly;
Param.Upwind=false;
Param.HyperVisc=true;
Param.HyperDCurl=2.e17; #1.e14*(dx/LRef)^3.2;
Param.HyperDGrad=2.e17;
Param.HyperDDiv=2.e17; # Scalars


# Output
Param.RadPrint=Param.H;
Param.Flat=true;
Param.vtkFileName="HeldSuarez";
Param.vtk=0;
vtkGrid=CGDycore.vtkCGGrid(CG,CGDycore.TransSphere,CGDycore.Topo,Param);
Param.cNames = [
  "Rho",
  "u",
  "v",
  "w",
  "Th"
]

# Initial conditions
U=zeros(CG.NumG,nz,Param.NumV);
U[:,:,Param.RhoPos]=CGDycore.Project(CGDycore.fRho,CG,Param);
(U[:,:,Param.uPos],U[:,:,Param.vPos])=CGDycore.ProjectVec(CGDycore.fVel,CG,Param);
U[:,:,Param.ThPos]=CGDycore.Project(CGDycore.fTheta,CG,Param).*U[:,:,Param.RhoPos];


# Integration
CFL=0.125;
dtau=500;
time=[0];

IntMethod="Rosenbrock";
# IntMethod="RungeKutta";
if strcmp(IntMethod,"Rosenbrock")
  dtau=600;
else
  dtau=8;
end
Param.RK=CGDycore.RungeKuttaMethod("RK4");
Param.ROS=CGDycore.RosenbrockMethod("SSP-Knoth");
SimDays=1000;
# SimDays=1;
PrintDay=10;
nIter=24*3600*SimDays/dtau;
PrintInt=24*3600*PrintDay/dtau;
# Print initial conditions
Param.vtk=CGDycore.vtkOutput(U,vtkGrid,CG,Param);
#
str = IntMethod
if str == "Rosenbrock"
    @time begin
      for i=1:nIter
        Δt = @elapsed begin
          U .= CGDycore.RosenbrockSchur(U,dtau,CGDycore.FcnNHCurlVec,CGDycore.JacSchur,CG,Param);
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
Param.vtk=CGDycore.vtkOutput(U,vtkGrid,CG,Param)

@info "Success!"
