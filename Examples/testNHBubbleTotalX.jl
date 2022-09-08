using CGDycore
using MPI
using Base

Base.@kwdef struct ParamStruct
  xC0=0.0
  zC0=2000.0
  rC0=2000.0
  Th0=300.0
  DeltaTh=2.0
  uMax=0.0
  vMax=0.0
  Stretch = false
end  

MPI.Init()
comm = MPI.COMM_WORLD
Proc = MPI.Comm_rank(comm) + 1
ProcNumber = MPI.Comm_size(comm)
print("$Proc: \n")
print("$ProcNumber: \n")

OrdPoly = 4
nx=40
ny=2
nz=80
OrdPolyZ=1
NumV = 5
NumTr = 0
Parallel = true


# Physical parameters
Phys=CGDycore.PhysParameters()


#ModelParameters
Param = ParamStruct()
Model = CGDycore.Model()
  Model.Equation="Compressible"
  Model.NumV=NumV
  Model.NumTr=NumTr
  Model.Problem="WarmBubble2Dx"
  Model.ProfRho="WarmBubble2Dx"
  Model.ProfTheta="WarmBubble2Dx"
  Model.ProfVel="Const"
  Model.RhoPos=1
  Model.uPos=2
  Model.vPos=3
  Model.wPos=4
  Model.ThPos=5
  Model.HorLimit = false
  Model.Upwind = true
  Model.Thermo = "TotalEnergy"

# Grid
Lx=2*10000.0
Ly=2000.0
H=10000.0
x0=-10000.0
y0=0.0;

Boundary = CGDycore.Boundary()
Boundary.WE="Period"
Boundary.SN="Period"
Boundary.BT=""
Topography=(TopoS="",H=H)

Grid=CGDycore.Grid(nz,Topography)
Grid=CGDycore.CartGrid(nx,ny,Lx,Ly,x0,y0,CGDycore.OrientFaceCart,Boundary,Grid)

if Parallel
  CellToProc = CGDycore.Decompose(Grid,ProcNumber)
  SubGrid = CGDycore.ConstructSubGrid(Grid,CellToProc,Proc)
  if Param.Stretch
    sigma = 1.0
    lambda = 3.16
    CGDycore.AddStretchICONVerticalGrid!(SubGrid,nz,H,sigma,lambda)
  else  
    CGDycore.AddVerticalGrid!(SubGrid,nz,H)
  end  
  Exchange = CGDycore.InitExchange(SubGrid,OrdPoly,CellToProc,Proc,ProcNumber,Parallel)
  Output=CGDycore.Output(Topography)
  Global = CGDycore.Global(SubGrid,Model,Phys,Output,Exchange,OrdPoly+1,nz,NumV,NumTr,())
  Global.Metric=CGDycore.Metric(OrdPoly+1,OrdPolyZ+1,SubGrid.NumFaces,nz)
else
  CellToProc=zeros(0)
  Proc = 0
  ProcNumber = 0
  CGDycore.AddVerticalGrid!(SubGrid,nz,H)
  Exchange = CGDycore.InitExchange(Grid,OrdPoly,CellToProc,Proc,ProcNumber,Parallel)
  Output=CGDycore.Output(Topography)
  Global = CGDycore.Global(Grid,Model,Phys,Output,Exchange,OrdPoly+1,nz,NumV,NumTr,())
  Global.Metric=CGDycore.Metric(OrdPoly+1,OrdPolyZ+1,Grid.NumFaces,nz)
end  
  (CG,Global)=CGDycore.Discretization(OrdPoly,OrdPolyZ,CGDycore.JacobiDG3,Global)
  Model.HyperVisc=true
  Model.HyperDCurl=1.e6
  Model.HyperDGrad=1.e6
  Model.HyperDDiv=1.e6

# Output
  Output.OrdPrint=CG.OrdPoly
  vtkGrid=CGDycore.vtkCGGrid(CG,CGDycore.TransCartX,CGDycore.Topo,Global)


  U = zeros(Float64,nz,CG.NumG,Model.NumV+Model.NumTr)
  U[:,:,Model.RhoPos]=CGDycore.Project(CGDycore.fRho,0.0,CG,Global,Param)
  (U[:,:,Model.uPos],U[:,:,Model.vPos])=CGDycore.ProjectVec(CGDycore.fVel,0.0,CG,Global,Param)
  if Global.Model.Thermo == "InternalEnergy"
    U[:,:,Model.ThPos]=CGDycore.Project(CGDycore.fIntEn,0.0,CG,Global,Param).*U[:,:,Model.RhoPos]
  elseif Global.Model.Thermo == "TotalEnergy"
    U[:,:,Model.ThPos]=CGDycore.Project(CGDycore.fTotEn,0.0,CG,Global,Param).*U[:,:,Model.RhoPos]
  else   
    U[:,:,Model.ThPos]=CGDycore.Project(CGDycore.fTheta,0.0,CG,Global,Param).*U[:,:,Model.RhoPos]
  end  
  if NumTr>0
    U[:,:,Model.RhoVPos+Model.NumV]=CGDycore.Project(CGDycore.fQv,0.0,CG,Global,Param).*U[:,:,Model.RhoPos]
  end   

# Output
  Output.vtkFileName=string("Bubble2DTotalx_")
  Output.vtk=0
  Output.Flat=false
  Output.H=H
  Output.cNames = [
    "Rho",
    "u",
    "v",
    "w",
    "Pres",
    "Th",
]
  Output.OrdPrint=CG.OrdPoly
  @show "Compute vtkGrid"
  Global.vtkCache = CGDycore.vtkInit(Output.OrdPrint,CGDycore.TransCartX,CG,Global)

  IntMethod="RungeKutta"
  IntMethod="RosenbrockD"
  IntMethod="LinIMEX"
  IntMethod="Rosenbrock"
  IntMethod="RungeKutta"
  if IntMethod == "Rosenbrock" || IntMethod == "RosenbrockD" || IntMethod == "RosenbrockSSP" || IntMethod == "LinIMEX"
    dtau = 0.2
  else
    dtau=0.3
  end
  Global.ROS=CGDycore.RosenbrockMethod("RODAS")
  Global.ROS=CGDycore.RosenbrockMethod("M1HOMME")
  Global.ROS=CGDycore.RosenbrockMethod("SSP-Knoth")
  Global.ROS=CGDycore.RosenbrockMethod("ROSEul")
  Global.RK=CGDycore.RungeKuttaMethod("RK4")
  Global.LinIMEX=CGDycore.LinIMEXMethod("ARS343")
  Global.LinIMEX=CGDycore.LinIMEXMethod("AR2")
  Global.LinIMEX=CGDycore.LinIMEXMethod("M1HOMME")

# Simulation period
  time=[0.0]
  SimSec=1000
  PrintSec=30
  PrintStartSec = 0
  nIter=ceil(SimSec/dtau)
  PrintInt=ceil(PrintSec/dtau)
  PrintStartInt=ceil(PrintStartSec/dtau)

  Global.Cache=CGDycore.CacheCreate(CG.OrdPoly+1,Global.Grid.NumFaces,CG.NumG,Global.Grid.nz,Model.NumV,Model.NumTr)

  if IntMethod == "Rosenbrock" || IntMethod == "RosenbrockD"
    Global.J = CGDycore.JStruct(CG.NumG,nz,Model.NumTr)
    Global.Cache.k=zeros(size(U[:,:,1:NumV+NumTr])..., Global.ROS.nStage);
    Global.Cache.fV=zeros(size(U))
    Global.Cache.Vn=zeros(size(U))
  elseif IntMethod == "RosenbrockSSP"
    Global.J = CGDycore.JStruct(CG.NumG,nz,Model.NumTr)
    Global.Cache.k=zeros(size(U[:,:,1:NumV])..., Global.ROS.nStage);
    Global.Cache.fV=zeros(size(U))
    Global.Cache.VS=zeros(size(U[:,:,NumV+1:end])..., Global.ROS.nStage+1);
    Global.Cache.fS=zeros(size(U[:,:,NumV+1:end])..., Global.ROS.nStage);
    Global.Cache.fRhoS=zeros(size(U[:,:,1])..., Global.ROS.nStage);
    Global.Cache.RhoS=zeros(size(U[:,:,1])..., Global.ROS.nStage+1);
  elseif IntMethod == "LinIMEX"
    Global.J = CGDycore.JStruct(CG.NumG,nz,Model.NumTr)
    Global.Cache.Ymyn=zeros(size(U[:,:,1:NumV+NumTr])..., Global.LinIMEX.nStage-1);
    Global.Cache.fV=zeros(size(U))
    Global.Cache.f=zeros(size(U)..., Global.LinIMEX.nStage)
    Global.Cache.fV=zeros(size(U))
    Global.Cache.Vn=zeros(size(U))
  elseif IntMethod == "RungeKutta"
    Global.Cache.f=zeros(size(U)..., Global.RK.nStage)
  end

# Print initial conditions
  @show "Print initial conditions"
  CGDycore.unstructured_vtkSphere(U,CGDycore.TransCartX,CG,Global,Proc,ProcNumber)

  @show "Choose integration method"
  if IntMethod == "Rosenbrock"
    @time begin
      for i=1:nIter
        Δt = @elapsed begin
          @show dtau
          CGDycore.RosenbrockSchur!(U,dtau,CGDycore.FcnNHCurlVecI!,CGDycore.JacSchur!,CG,Global,Param);
          MassRho = CGDycore.GlobalIntegral(U[:,:,1],CG,Global)
          @show MassRho
          time[1] += dtau
          if mod(i,PrintInt) == 0 && i >= PrintStartInt
            CGDycore.unstructured_vtkSphere(U,CGDycore.TransCartX,CG,Global,Proc,ProcNumber)
          end
        end
        percent = i/nIter*100
        @info "Iteration: $i took $Δt, $percent% complete"
      end
    end
  elseif IntMethod == "RosenbrockD"
    @time begin
      for i=1:nIter
        Δt = @elapsed begin
          CGDycore.RosenbrockDSchur!(U,dtau,CGDycore.FcnNHCurlVec!,CGDycore.JacSchur!,CG,Global);
          time[1] += dtau
          if mod(i,PrintInt) == 0 && i >= PrintStartInt
            CGDycore.unstructured_vtkSphere(U,CGDycore.TransCartX,CG,Global,Proc,ProcNumber)
          end
        end
        percent = i/nIter*100
        @info "Iteration: $i took $Δt, $percent% complete"
      end
    end  
  elseif IntMethod == "RosenbrockSSP"
    @time begin
      for i=1:nIter
        Δt = @elapsed begin
          CGDycore.RosenbrockSchurSSP!(U,dtau,CGDycore.FcnNHCurlVec!,CGDycore.JacSchur!,CG,Global);
          time[1] += dtau
          if mod(i,PrintInt) == 0 && i >= PrintStartInt
            CGDycore.unstructured_vtkSphere(U,CGDycore.TransCartX,CG,Global,Proc,ProcNumber)
          end
        end
        percent = i/nIter*100
        @info "Iteration: $i took $Δt, $percent% complete"
      end
    end
  elseif IntMethod == "LinIMEX"
    @time begin
      for i=1:nIter
        Δt = @elapsed begin
          CGDycore.LinIMEXSchur!(U,dtau,CGDycore.FcnNHCurlVec!,CGDycore.JacSchur!,CG,Global);
          time[1] += dtau
          if mod(i,PrintInt) == 0 && i >= PrintStartInt
            CGDycore.unstructured_vtkSphere(U,CGDycore.TransCartX,CG,Global,Proc,ProcNumber)
          end
        end
        percent = i/nIter*100
        @info "Iteration: $i took $Δt, $percent% complete"
      end
    end
  elseif IntMethod == "RungeKutta"
    @time begin
      for i=1:nIter
        Δt = @elapsed begin
          CGDycore.RungeKuttaExplicit!(U,dtau,CGDycore.FcnNHCurlVecI!,CG,Global,Param)

          time[1] += dtau
          if mod(i,PrintInt)==0 && i >= PrintStartInt
            CGDycore.unstructured_vtkSphere(U,CGDycore.TransCartX,CG,Global,Proc,ProcNumber)
          end
        end
        percent = i/nIter*100
        @info "Iteration: $i took $Δt, $percent% complete"
      end
    end
  else
    error("Bad IntMethod")
  end
