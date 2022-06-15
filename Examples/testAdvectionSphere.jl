#function testNHDensityCurrent()

using CGDycore

OrdPoly = 4
nz = 12

OrdPolyZ=1
nPanel = 8
NF = 6 * nPanel * nPanel
NumV = 4
NumTr = 2
H = 1.2e4



# Physical parameters
  Phys=CGDycore.PhysParameters();

#ModelParameters
  T_0 = 300.0
  ScaleHeight = Float64(Phys.Rd * T_0 / Phys.Grav)
  tau = 1036800.0
  omega_0 = 23000 * pi / tau
  Param=(xC = 0.0,
         H = H,
         R_t = Phys.RadEarth / 2.0,
         Z_t = 1000.0,
         z_c = 5.0e3,
         p_top = 25494.4,
         T_0 = T_0,
         ScaleHeight = ScaleHeight,
         b = 0.2,
         Lon_c1 = 150.0/360.0 * 2 * pi,
         Lon_c2 = 210.0/360.0 * 2 * pi,
         Lat_c = 0.0,
         tau = tau,
         omega_0 = omega_0,
         TimeDependent = true,)

  Model = CGDycore.Model(Param)

# Grid
  Topography=(TopoS="",H=H,Rad=Phys.RadEarth)
  Grid=CGDycore.Grid(nz,Topography)
  Grid=CGDycore.CubedGrid(nPanel,CGDycore.OrientFaceSphere,Phys.RadEarth,Grid)
  CGDycore.AddVerticalGrid!(Grid,nz,H)
  Output=CGDycore.Output(Topography)
  Global = CGDycore.Global(Grid,Model,Phys,Output,OrdPoly+1,nz,NumV,NumTr,())
  Global.Metric=CGDycore.Metric(OrdPoly+1,OrdPolyZ+1,Grid.NumFaces,nz)

# Discretization
  (CG,Global)=CGDycore.Discretization(OrdPoly,OrdPolyZ,CGDycore.JacobiSphere3,Global);
  Model.HyperVisc=true
  Model.HyperDDiv=1.e16
  Model.Upwind=true
  Model.HorLimit=true

# Initial conditions 
Model.NumV = NumV
Model.NumTr = NumTr
U=zeros(nz,CG.NumG,Model.NumV+Model.NumTr);
Model.Problem="AdvectionSphereDCMIP"
Model.ProfRho="AdvectionSphereDCMIP"
Model.ProfVel="AdvectionSphereDCMIP"
Model.ProfVelW="AdvectionSphereDCMIP"
Model.RhoPos=1;
Model.uPos=2;
Model.vPos=3;
Model.wPos=4;
U[:,:,Model.RhoPos]=CGDycore.Project(CGDycore.fRho,0.0,CG,Global);
(U[:,:,Model.uPos],U[:,:,Model.vPos])=CGDycore.ProjectVec(CGDycore.fVel,0.0,CG,Global);
U[:,:,Model.wPos]=CGDycore.ProjectW(CGDycore.fVelW,0.0,CG,Global)
Model.ProfTr="AdvectionSphereDCMIPQ1"
U[:,:,Model.NumV+1]=CGDycore.Project(CGDycore.fTr,0.0,CG,Global).*U[:,:,Model.RhoPos];
Model.ProfTr="AdvectionSphereDCMIPQ2"
U[:,:,Model.NumV+2]=CGDycore.Project(CGDycore.fTr,0.0,CG,Global).*U[:,:,Model.RhoPos];

# Output
  Output.vtkFileName="AdvSphere"
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
   "Tr1",
   "Tr2",
  ]
  Output.OrdPrint=CG.OrdPoly
  vtkGrid=CGDycore.vtkCGGrid(CG,CGDycore.TransSphereX,CGDycore.Topo,Global);
  Output.Flat=false
  vtkGridC=CGDycore.vtkGridElem(CG,CGDycore.TransSphereX,CGDycore.Topo,Global);
  Output.Flat=true
# Proc=zeros(NF)
# for i=1:NF
#  if i <= NF/2   
#    Proc[i] = 1
#  else
#    Proc[i] = 2
#  end  
# end 
# CGDycore.vtkElem(Proc,vtkGridC);
# OwnGrid = CGDycore.ConstructSubGrid(Grid,Proc,1)
# stop 

IntMethod="RungeKutta";
IntMethod="SSPRungeKutta";
dtau = 10.0 * 60.0

Global.RK=CGDycore.RungeKuttaMethod("SSP3");
Global.SSP=CGDycore.SSPRungeKuttaMethod("SSP-Knoth");

time=[0.0];
SimDays=12
PrintDay=.5
nIter=24*3600*SimDays/dtau
PrintInt=ceil(24*3600*PrintDay/dtau)

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
