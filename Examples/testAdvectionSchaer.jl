#function testNHDensityCurrent()

using CGDycore

OrdPoly = 1
OrdPolyZ=1
nx=75;  # 300
ny=2;
nz = 50
NF = nx *ny
NumV = 4
NumTr = 1
Lx = 100000.0
H = 25000.0
x0 = -50000.0


# Physical parameters
Phys=CGDycore.PhysParameters();

#ModelParameters
Param=(xC = -50000.0,
       Ax = 25000.0,
       zC = 15000.0, #9000.0
       Az = 3000.0,
       uMax=10,
       vMax=0,
       q0 = 1.0,
       z1 = 4000.0,
       z2 = 9000.0, #5000.0
       )
Model = CGDycore.Model(Param)

# Grid
Ly=200000.0;
y0=0.0;

Boundary = (;WE="Period", SN="Period")
Topography=(TopoS="AdvectionSchaer",h0=3000.0,lambda=8000.0,a=25000.0,H=H)
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
Model.NumV = NumV
Model.NumTr = NumTr
U=zeros(nz,CG.NumG,Model.NumV+Model.NumTr);
Model.ProfRho="AdvectionSchaer"
Model.ProfVel="AdvectionSchaer"
Model.ProfTr="AdvectionSchaer"
Model.RhoPos=1;
Model.uPos=2;
Model.vPos=3;
Model.wPos=4;
U[:,:,Model.RhoPos]=CGDycore.Project(CGDycore.fRho,CG,Global);
(U[:,:,Model.uPos],U[:,:,Model.vPos])=CGDycore.ProjectVec(CGDycore.fVel,CG,Global);
U[:,:,Model.NumV+1]=CGDycore.Project(CGDycore.fTr,CG,Global).*U[:,:,Model.RhoPos];

# Output
Output.vtkFileName="AdvSchaer"
Output.vtk=0;
Output.cNames = [
  "Rho",
  "u",
  "v",
  "w",
  "Tr1",
]
Output.OrdPrint=CG.OrdPoly
vtkGrid=CGDycore.vtkCGGrid(CG,CGDycore.TransCartX,CGDycore.Topo,Global);


IntMethod="RungeKutta";
if IntMethod == "Rosenbrock"
  dtau=30;
else
  dtau=1;
end
Global.RK=CGDycore.RungeKuttaMethod("RK1");
time=[0.0];
EndTime=10;
nIter=EndTime/dtau;
PrintTime=1;
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
          CGDycore.RungeKuttaExplicit!(U,dtau,CGDycore.FcnTracer!,CG,Global);

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
