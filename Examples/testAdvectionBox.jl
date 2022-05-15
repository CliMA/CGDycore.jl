#function testNHDensityCurrent()

using CGDycore

OrdPoly = 4
OrdPolyZ = 1
nx = 16
ny = 16
nz = 64
NF = nx *ny
NumV = 4
NumTr = 1
Lx = 4 * pi
Ly = 4 * pi
H = 4 * pi
x0 = -2 * pi
y0 = -2 * pi


# Physical parameters
Phys=CGDycore.PhysParameters();

#ModelParameters
Param=(xC = 0.0,
       yC = 0.0,
       uMax = pi / 2,
       time = 0,
       EndTime = 2 * pi,
       TimeDependent = true,
       H = H,
       xB1 = x0 + Lx / 4,
       yB1 = y0 + Ly / 2,
       zB1 = H / 2,
       xB2 = x0 + 3.0 * Lx / 4,
       yB2 = y0 + Ly / 2,
       zB2 = H / 2,
       r0 = Lx / 6,)
Model = CGDycore.Model(Param)

# Grid

Boundary = (;WE="FreeSlip", SN="FreeSlip")
Topography=(TopoS="",h0=3000.0,lambda=8000.0,a=25000.0,H=H)
Grid=CGDycore.Grid(nz,Topography)
Grid=CGDycore.CartGrid(nx,ny,Lx,Ly,x0,y0,CGDycore.OrientFaceCart,Boundary,Grid);
CGDycore.AddVerticalGrid!(Grid,nz,H)

Output=CGDycore.Output(Topography)

Global = CGDycore.Global(Grid,Model,Phys,Output,OrdPoly+1,nz,NumV,NumTr,())

Global.Metric=CGDycore.Metric(OrdPoly+1,OrdPolyZ+1,Grid.NumFaces,nz)


# Discretization
(CG,Global)=CGDycore.Discretization(OrdPoly,OrdPolyZ,CGDycore.JacobiDG3,Global);
Model.HyperVisc=true;
Model.HyperDCurl=0.e4;
Model.HyperDGrad=0.e4;
Model.HyperDDiv=1.e-3;
Model.Upwind=true


# Initial conditions 
Model.NumV = NumV
Model.NumTr = NumTr
U=zeros(nz,CG.NumG,Model.NumV+Model.NumTr);
Model.Problem="AdvectionTestDeform"
Model.ProfRho="AdvectionTestDeform"
Model.ProfVel="AdvectionTestDeform"
Model.ProfVelW="AdvectionTestDeform"
Model.ProfTr="AdvectionTestDeform"
Model.RhoPos=1;
Model.uPos=2;
Model.vPos=3;
Model.wPos=4;
U[:,:,Model.RhoPos]=CGDycore.Project(CGDycore.fRho,0.0,CG,Global);
(U[:,:,Model.uPos],U[:,:,Model.vPos])=CGDycore.ProjectVec(CGDycore.fVel,0.0,CG,Global);
U[:,:,Model.wPos]=CGDycore.ProjectW(CGDycore.fVelW,0.0,CG,Global);
U[:,:,Model.NumV+1]=CGDycore.Project(CGDycore.fTr,0.0,CG,Global).*U[:,:,Model.RhoPos];

# Output
Output.vtkFileName="AdvDeform"
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


IntMethod="SSPRungeKutta";
dtau=Param.EndTime / 1000 #800

Global.RK=CGDycore.RungeKuttaMethod("RK4");
Global.SSP=CGDycore.SSPRungeKuttaMethod("SSP-Knoth");

time=[0.0];
EndTime=Param.EndTime
nIter=EndTime/dtau;
PrintTime=0.01 * EndTime
PrintInt=floor(PrintTime/dtau)

Global.Cache=CGDycore.CacheCreate(CG.OrdPoly+1,Global.Grid.NumFaces,CG.NumG,Global.Grid.nz,Model.NumV,Model.NumTr)
str = IntMethod
if str == "SSPRungeKutta"
  Global.Cache.fV=zeros(size(U))
  Global.Cache.VS=zeros(size(U[:,:,NumV+1:end])..., Global.SSP.nStage+1);
  Global.Cache.fS=zeros(size(U[:,:,NumV+1:end])..., Global.SSP.nStage);
  Global.Cache.fRhoS=zeros(size(U[:,:,1])..., Global.SSP.nStage);
  Global.Cache.RhoS=zeros(size(U[:,:,1])..., Global.SSP.nStage+1);
elseif str == "RungeKutta"
  Global.Cache.f=zeros(size(U)..., Global.RK.nStage);
end

# Print initial conditions
Global.Output.vtk=CGDycore.vtkOutput(U,vtkGrid,CG,Global);


if str == "SSPRungeKutta"
    @time begin
      for i=1:nIter
        Δt = @elapsed begin
          CGDycore.SSPRungeKutta!(time[1],U,dtau,CGDycore.FcnTracer!,CG,Global);
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
          CGDycore.RungeKuttaExplicit!(time[1],U,dtau,CGDycore.FcnTracer!,CG,Global);

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
Global.Output.vtk=CGDycore.vtkOutput(U,vtkGrid,CG,Global);
