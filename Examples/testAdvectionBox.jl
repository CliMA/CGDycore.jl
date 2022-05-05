#function testNHDensityCurrent()

using CGDycore

OrdPoly = 4
OrdPolyZ=1
nx=60;
ny=2;
nz = 40
NF = nx *ny
NumV = 5
NumTr = 1
Lx=6000000.0;
H=15600.0
x0=-3000000.0;


# Physical parameters
Phys=CGDycore.PhysParameters();

#ModelParameters
Param=(NBr=1.e-2,
       Th0=300,
       T0=300.0,
       uMax=10,
       vMax=0,
       xC = x0 + 0.25 * Lx,
       zC = 0.2 * H,
       xH = 0.1 * Lx,
       zH = 0.1 * H,
       )
Model = CGDycore.Model(Param)

# Grid
Ly=200000.0;
y0=0.0;

Boundary = (;WE="Period", SN="Period")
Topography=(TopoS="AgnesiCartX",hC=400.0,x0C=0.0,aC=100000.0,H=H)
Grid=CGDycore.Grid(nz,Topography)
Grid=CGDycore.CartGrid(nx,ny,Lx,Ly,x0,y0,CGDycore.OrientFaceCart,Boundary,Grid);
CGDycore.AddVerticalGrid!(Grid,nz,H)

Output=CGDycore.Output(Topography)

Global = CGDycore.Global(Grid,Model,Phys,Output,OrdPoly+1,nz,NumV,NumTr,())

Global.Metric=CGDycore.Metric(OrdPoly+1,OrdPolyZ+1,Grid.NumFaces,nz)


# Discretization
(CG,Global)=CGDycore.Discretization(OrdPoly,OrdPolyZ,CGDycore.JacobiDG3,Global);
Model.HyperVisc=true;
Model.HyperDCurl=1.e4;
Model.HyperDGrad=1.e4;
Model.HyperDDiv=1.e4;
Model.Upwind=true


# Initial conditions 
Model.Damping=true
Model.StrideDamp=6000;
Model.Relax=1.e0;
Model.NumV = NumV
Model.NumTr = NumTr
U=zeros(nz,CG.NumG,Model.NumV + Model.NumTr);
Model.ProfRho="GravityHill";
Model.ProfTheta="GravityHill";
Model.ProfVel="Const";
Model.ProfTr="Cylinder";
Model.RhoPos=1;
Model.uPos=2;
Model.vPos=3;
Model.wPos=4;
Model.ThPos=5;
U[:,:,Model.RhoPos]=CGDycore.Project(CGDycore.fRho,CG,Global);
(U[:,:,Model.uPos],U[:,:,Model.vPos])=CGDycore.ProjectVec(CGDycore.fVel,CG,Global);
U[:,:,Model.ThPos]=CGDycore.Project(CGDycore.fTheta,CG,Global).*U[:,:,Model.RhoPos];
U[:,:,Model.NumV + 1]=CGDycore.Project(CGDycore.fTr,CG,Global).*U[:,:,Model.RhoPos];

# Output
Output.vtkFileName="Hill_L";
Output.vtk=0;
Output.cNames = [
  "Rho",
  "u",
  "v",
  "w",
  "Th",
  "Tr1",
]
Output.OrdPrint=CG.OrdPoly
vtkGrid=CGDycore.vtkCGGrid(CG,CGDycore.TransCartX,CGDycore.Topo,Global);


IntMethod="RungeKutta";
IntMethod="Rosenbrock";
if IntMethod == "Rosenbrock"
  dtau=30;
else
  dtau=0.4;
end
Global.ROS=CGDycore.RosenbrockMethod("SSP-Knoth");
Global.RK=CGDycore.RungeKuttaMethod("RK4");
time=[0.0];
EndTime=216000;
nIter=EndTime/dtau;
PrintTime=10000;
PrintInt=floor(PrintTime/dtau)

Global.Cache=CGDycore.CacheCreate(CG.OrdPoly+1,Global.Grid.NumFaces,CG.NumG,Global.Grid.nz,Model.NumV,Model.NumTr)
str = IntMethod
if str == "Rosenbrock"
  Global.J = CGDycore.JStruct(CG.NumG,nz)
  Global.Cache.k=zeros(size(U[:,:,1:NumV])..., Global.ROS.nStage);
  Global.Cache.fV=zeros(size(U))
  Global.Cache.VS=zeros(size(U[:,:,NumV+1:end])..., Global.ROS.nStage+1);
  Global.Cache.fS=zeros(size(U[:,:,NumV+1:end])..., Global.ROS.nStage);
  Global.Cache.fRhoS=zeros(size(U[:,:,1])..., Global.ROS.nStage);
  Global.Cache.RhoS=zeros(size(U[:,:,1])..., Global.ROS.nStage+1);
elseif str == "RungeKutta"
  Global.Cache.f=zeros(size(U)..., Global.RK.nStage);
end

# Print initial conditions
Global.Output.vtk=CGDycore.vtkOutput(U,vtkGrid,CG,Global);


if str == "Rosenbrock"
    @time begin
      for i=1:nIter
        Δt = @elapsed begin
          CGDycore.RosenbrockSchurSSP!(U,dtau,CGDycore.FcnNHCurlVec!,CGDycore.JacSchur!,CG,Global);
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
          CGDycore.RungeKuttaExplicit!(U,dtau,CGDycore.FcnNHCurlVec!,CG,Global);

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
