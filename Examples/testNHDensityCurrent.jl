#function testNHDensityCurrent()

using CGDycore

OrdPoly = 4
OrdPolyZ=1
nx=80;
ny=2;
nz = 60
NF = nx *ny


# Physical parameters
Phys=CGDycore.PhysParameters();

#ModelParameters
Param=(T0=300.0,
       DeltaT=-15.0,
       xC0=0.0,
       zC0=3000,
       xrC0=4000.0,
       zrC0=2000.0,
       uMax=0,
       vMax=0,
       )
Model = CGDycore.Model(Param)

# Grid
Lx=2*25600.0;
Ly=2000.0;
H=6400.0;
x0=-25600.0;
y0=0.0;

Boundary = (;WE="Period", SN="Period")
Topography=(TopoS="",hC=400.0,x0C=0.0,aC=1000.0,H=H)
Grid=CGDycore.Grid(nz,Topography)
Grid=CGDycore.CartGrid(nx,ny,Lx,Ly,x0,y0,CGDycore.OrientFaceCart,Boundary,Grid);
CGDycore.AddVerticalGrid!(Grid,nz,H)

Output=CGDycore.Output(Topography)

Global = CGDycore.Global(Grid,Model,Phys,Output)
Global.Metric=CGDycore.Metric(OrdPoly+1,OrdPolyZ+1,Grid.NumFaces,nz)


# Discretization
(CG,Global)=CGDycore.Discretization(OrdPoly,OrdPolyZ,CGDycore.JacobiDG3,Global);
Model.HyperVisc=true;
Model.HyperDCurl=1.e6;
Model.HyperDGrad=1.e6;
Model.HyperDDiv=1.e6;
Model.Upwind=true


# Initial conditions 
Model.NumV=5;
U=zeros(nz,CG.NumG,Model.NumV);
Model.ProfRho="DensityCurrent";
Model.ProfTheta="DensityCurrent";
Model.ProfVel="Const";
Model.RhoPos=1;
Model.uPos=2;
Model.vPos=3;
Model.wPos=4;
Model.ThPos=5;
U[:,:,Model.RhoPos]=CGDycore.Project(CGDycore.fRho,CG,Global);
(U[:,:,Model.uPos],U[:,:,Model.vPos])=CGDycore.ProjectVec(CGDycore.fVel,CG,Global);
U[:,:,Model.ThPos]=CGDycore.Project(CGDycore.fTheta,CG,Global).*U[:,:,Model.RhoPos];

# Output
Output.vtkFileName="DensityCurrentU";
Output.vtk=0;
Output.cNames = [
  "Rho",
  "u",
  "v",
  "w",
  "Th"
]
vtkGrid=CGDycore.vtkCGGrid(CG,CGDycore.TransCart,CGDycore.Topo,Global);


#IntMethod="Rosenbrock";
IntMethod="RungeKutta";
if IntMethod == "Rosenbrock"
  dtau=30;
else
  dtau=0.3;
end
Global.ROS=CGDycore.RosenbrockMethod("SSP-Knoth");
Global.RK=CGDycore.RungeKuttaMethod("RK4");
time=[0.0];
EndTime=900;
nIter=EndTime/dtau;
PrintTime=90;
PrintInt=floor(PrintTime/dtau);


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
          CGDycore.RosenbrockSchur!(U,dtau,CGDycore.FcnNHCurlVec!,CGDycore.JacSchur,CG,Global);
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
