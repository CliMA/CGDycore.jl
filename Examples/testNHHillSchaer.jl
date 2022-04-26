#function testNHBaroWaveSphere

using CGDycore
using SimpleGraphAlgorithms

OrdPoly = 4
OrdPolyZ=1
nx = 60
ny = 2
nz = 80
NF = nx * ny
NumV = 5


# Physical parameters
Phys=CGDycore.PhysParameters();

#ModelParameters
Param=(
       Deep=false,
       NBr=1.e-2,
       Th0=300.0,
       uMax=20,
       vMax=0,
       TEq=300.0,
       )

Model = CGDycore.Model(Param)
Model.Coriolis=false
Model.CoriolisType="Sphere";

# Grid
Lx=2*30000.0;
Ly=6000.0;
H=25000.0;
x0=-30000.0;
y0=0.0;

Boundary = (;WE="Period", SN="Period")
Topography=(TopoS="SchaerCart",
            H=H,
            d0=5000.0,
            ksi0=4000.0,
            h0=250.0,
            )
Grid=CGDycore.Grid(nz,Topography)
Boundary = (;WE="Period", SN="Period")
Grid=CGDycore.CartGrid(nx,ny,Lx,Ly,x0,y0,CGDycore.OrientFaceCart,Boundary,Grid);
Grid.colors = CGDycore.Coloring(Grid)

CGDycore.AddVerticalGrid!(Grid,nz,H)

Output=CGDycore.Output(Topography)

NumTr=0
Global = CGDycore.Global(Grid,Model,Phys,Output,OrdPoly+1,nz,NumV,NumTr,())
#Global.ThreadCache=CGDycore.CreateCache(OrdPoly+1,nz,NumV)
Global.Metric=CGDycore.Metric(OrdPoly+1,OrdPolyZ+1,Grid.NumFaces,nz)


# Discretization
(CG,Global)=CGDycore.Discretization(OrdPoly,OrdPolyZ,CGDycore.JacobiDG3,Global);
Model.HyperVisc=true;
Model.HyperDCurl=2.e6
Model.HyperDGrad=2.e6
Model.HyperDDiv=2.e6
Model.StrideDamp=10000;
Model.Relax=1.0e-4;
Model.Damping=true;


# Initial conditions 
Model.NumV=NumV;
U=zeros(nz,CG.NumG,Model.NumV);
Model.ProfRho="IsoThermal"
Model.ProfTheta="IsoThermal"
Model.ProfVel="Const"
Model.RhoPos=1;
Model.uPos=2;
Model.vPos=3;
Model.wPos=4;
Model.ThPos=5;
U[:,:,Model.RhoPos]=CGDycore.Project(CGDycore.fRho,CG,Global);
(U[:,:,Model.uPos],U[:,:,Model.vPos])=CGDycore.ProjectVec(CGDycore.fVel,CG,Global);
U[:,:,Model.ThPos]=CGDycore.Project(CGDycore.fTheta,CG,Global).*U[:,:,Model.RhoPos];

# Output
Output.vtkFileName="SchaerCartIsoT";
Output.vtk=0;
Output.Flat=true
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
vtkGrid=CGDycore.vtkCGGrid(CG,CGDycore.TransCartX,CGDycore.Topo,Global);


IntMethod="RungeKutta";
IntMethod="Rosenbrock";
if IntMethod == "Rosenbrock"
  dtau=.4;
else
  dtau=.4;
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
