#function testNHBaroWaveSphere

using CGDycore

OrdPoly = 4
OrdPolyZ=1
nz = 20
nPanel = 8
NF = 6 * nPanel * nPanel


# Physical parameters
Phys=CGDycore.PhysParameters();

#ModelParameters
Param=(T0E=310.0,
       T0P=240.0,
       B=2.0,
       K=3.0,
       LapseRate=0.005,
       U0=-0.5,
       PertR=1.0/6.0,
       Up=1.0,
       PertExpR=0.1,
       PertLon=pi/9.0,
       PertLat=2.0*pi/9.0,
       PertZ=15000.0,
       NBr=1.e-2,
       DeltaT=1,
       ExpDist=5,
       T0=300,
       Deep=false,
       )

Model = CGDycore.Model(Param)
Model.Coriolis=true
Model.CoriolisType="Sphere";

# Grid
H = 30000.0
Topography=(TopoS="",H=H,Rad=Phys.RadEarth)
Grid=CGDycore.Grid(nz,Topography)
Grid=CGDycore.CubedGrid(nPanel,CGDycore.OrientFaceSphere,Phys.RadEarth,Grid);
CGDycore.AddVerticalGrid!(Grid,nz,H)

Output=CGDycore.Output(Topography)

Global = CGDycore.Global(Grid,Model,Phys,Output)
Global.Metric=CGDycore.Metric(OrdPoly+1,OrdPolyZ+1,Grid.NumFaces,nz)


# Discretization
(CG,Global)=CGDycore.Discretization(OrdPoly,OrdPolyZ,CGDycore.JacobiSphere3,Global);
Model.HyperVisc=true;
Model.HyperDCurl=2.e17/4; #1.e14*(dx/LRef)^3.2;
Model.HyperDGrad=2.e17/4;
Model.HyperDDiv=2.e17/4; # Scalars


# Initial conditions 
Model.NumV=5;
U=zeros(nz,CG.NumG,Model.NumV);
Model.ProfRho="BaroWaveSphere"
Model.ProfTheta="BaroWaveSphere"
Model.ProfVel="BaroWaveSphere"
Model.RhoPos=1;
Model.uPos=2;
Model.vPos=3;
Model.wPos=4;
Model.ThPos=5;
U[:,:,Model.RhoPos]=CGDycore.Project(CGDycore.fRho,CG,Global);
(U[:,:,Model.uPos],U[:,:,Model.vPos])=CGDycore.ProjectVec(CGDycore.fVel,CG,Global);
U[:,:,Model.ThPos]=CGDycore.Project(CGDycore.fTheta,CG,Global).*U[:,:,Model.RhoPos];

# Output
Output.vtkFileName="BaroWave";
Output.vtk=0;
Output.Flat=true
Output.nPanel=nPanel
Output.RadPrint=H
Output.H=H
Output.cNames = [
  "Rho",
  "u",
  "v",
  "w",
  "Th"
]
vtkGrid=CGDycore.vtkCGGrid(CG,CGDycore.TransSphere,CGDycore.Topo,Global);


IntMethod="RungeKutta";
IntMethod="Rosenbrock";
if IntMethod == "Rosenbrock"
  dtau=200;
else
  dtau=1;
end
Global.ROS=CGDycore.RosenbrockMethod("SSP-Knoth");
Global.RK=CGDycore.RungeKuttaMethod("RK4");
time=[0.0];
SimDays=10;
PrintDay=.25;
nIter=24*3600*SimDays/dtau;
PrintInt=24*3600*PrintDay/dtau;


Global.Cache=CGDycore.CacheCreate(CG.OrdPoly+1,Global.Grid.NumFaces,CG.NumG,Global.Grid.nz,Model.NumV)
str = IntMethod
if str == "Rosenbrock"
  Global.J = CGDycore.JStruct(CG.NumG,nz)
  Global.Cache.k=zeros(size(U)..., Global.ROS.nStage);
  Global.Cache.fV=zeros(size(U))
elseif str == "RungeKutta"
  Global.Cache.f=zeros(size(U)..., Global.RK.nStage);
end

# Print initial conditions
Global.Output.vtk=CGDycore.vtkOutput(U,vtkGrid,CG,Global);

if str == "Rosenbrock"
    @time begin
      for i=1:nIter
        Δt = @elapsed begin
          CGDycore.RosenbrockSchur!(U,dtau,CGDycore.FcnNHCurlVec!,CGDycore.JacSchur!,CG,Global);
          time[1] += dtau;
          if mod(i,PrintInt)==0
            Global.Output.vtk=CGDycore.vtkOutput(U,vtkGrid,CG,Global);
          end
        end
        percent = i/nIter*100
        @info "Iteration: $i took $Δt, $percent% complete"
      end
    end

elseif str == "RungeKutta"
    @time begin
      for i=1:nIter
        Δt = @elapsed begin
          @time CGDycore.RungeKuttaExplicit!(U,dtau,CGDycore.FcnNHCurlVec!,CG,Global);

          time[1] += dtau;
          if mod(i,PrintInt)==0
            Global.Output.vtk=CGDycore.vtkOutput(U,vtkGrid,CG,Global);
          end
        end
        percent = i/nIter*100
        @info "Iteration: $i took $Δt, $percent% complete"
      end
    end
else
  error("Bad str")
end
