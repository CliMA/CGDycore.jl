#function testNHBaroWaveSphere

using CGDycore
using SimpleGraphAlgorithms

OrdPoly = 4
OrdPolyZ=1
nz = 60
nPanel = 20
NF = 6 * nPanel * nPanel
NumV = 5


# Physical parameters
Phys=CGDycore.PhysParameters();

#ModelParameters
X = 166.7
Param=(
       PertZ=15000.0,
       Deep=false,
       RadEarth=Phys.RadEarth/X,
       TEq=300.0,
       uEq=20,
       )

Model = CGDycore.Model(Param)
Model.Coriolis=false
Model.CoriolisType="Sphere";

# Grid
H = 20000.0
Topography=(TopoS="SchaerSphereRidge",
            H=H,
            Rad=Param.RadEarth,
            d0=5000.0,
            ksi0=4000.0,
            PertLon=pi,
            PertLat=0.0,
            h0=250.0,
            RadEarth=Param.RadEarth,
            )
Grid=CGDycore.Grid(nz,Topography)
Grid=CGDycore.CubedGrid(nPanel,CGDycore.OrientFaceSphere,Param.RadEarth,Grid);
Grid.colors = CGDycore.Coloring(Grid)

CGDycore.AddVerticalGrid!(Grid,nz,H)

Output=CGDycore.Output(Topography)

NumTr=0
Global = CGDycore.Global(Grid,Model,Phys,Output,OrdPoly+1,nz,NumV,NumTr,())
#Global.ThreadCache=CGDycore.CreateCache(OrdPoly+1,nz,NumV)
Global.Metric=CGDycore.Metric(OrdPoly+1,OrdPolyZ+1,Grid.NumFaces,nz)


# Discretization
(CG,Global)=CGDycore.Discretization(OrdPoly,OrdPolyZ,CGDycore.JacobiSphere3,Global);
Model.HyperVisc=true;
Model.HyperDCurl=2.e6
Model.HyperDGrad=2.e6
Model.HyperDDiv=2.e6
Model.StrideDamp=10000;
Model.Relax=1.0e-3;
Model.Damping=true;


# Initial conditions 
Model.NumV=NumV;
U=zeros(nz,CG.NumG,Model.NumV);
Model.ProfRho="SchaerSphere"
Model.ProfTheta="SchaerSphere"
Model.ProfVel="SchaerSphere"
Model.RhoPos=1;
Model.uPos=2;
Model.vPos=3;
Model.wPos=4;
Model.ThPos=5;
U[:,:,Model.RhoPos]=CGDycore.Project(CGDycore.fRho,CG,Global);
(U[:,:,Model.uPos],U[:,:,Model.vPos])=CGDycore.ProjectVec(CGDycore.fVel,CG,Global);
U[:,:,Model.ThPos]=CGDycore.Project(CGDycore.fTheta,CG,Global).*U[:,:,Model.RhoPos];

# Output
Output.vtkFileName="SchaerSphere";
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
Output.OrdPrint=CG.OrdPoly
vtkGrid=CGDycore.vtkCGGrid(CG,CGDycore.TransSphereX,CGDycore.Topo,Global);


IntMethod="RungeKutta";
IntMethod="Rosenbrock";
if IntMethod == "Rosenbrock"
  dtau=.5;
else
  dtau=1;
end
Global.ROS=CGDycore.RosenbrockMethod("SSP-Knoth");
Global.RK=CGDycore.RungeKuttaMethod("RK4");
time=[0.0];
SimHours=1;
PrintHours=.05;
nIter=SimHours*3600/dtau;
PrintInt=PrintHours*3600/dtau;


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
