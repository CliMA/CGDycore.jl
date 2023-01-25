using CGDycore
using MPI
using Base

Base.@kwdef struct ParamStruct
  Deep=false
  NBr=1.e-2
  Th0=300.0
  uMax=20
  vMax=0
  TEq=300.0
  Stretch = false
end  

MPI.Init()
comm = MPI.COMM_WORLD
Proc = MPI.Comm_rank(comm) + 1
ProcNumber = MPI.Comm_size(comm)
print("$Proc: \n")
print("$ProcNumber: \n")

OrdPoly = 4
nx = 60
ny = 2
nz = 80
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
  Model.Problem="SchaerCart"
  Model.ProfRho="IsoThermal"
  Model.ProfTheta="IsoThermal"
  Model.ProfVel="Const"
  Model.RhoPos=1
  Model.uPos=2
  Model.vPos=3
  Model.wPos=4
  Model.ThPos=5
  Model.HorLimit = false
  Model.Upwind = true
  Model.Thermo = "" #"TotalEnergy"
  Model.Damping=true
  Model.StrideDamp=10000
  Model.Relax=1.0e-2

# Grid
Lx=2*30000.0;
Ly=6000.0;
H=25000.0;
x0=-30000.0;
y0=0.0;

Boundary = CGDycore.Boundary()
Boundary.WE="Period"
Boundary.SN="Period"
Boundary.BT=""
Topography=(TopoS="SchaerCart",
            H=H,
            P1=5000.0,
            P2=4000.0,
            P3=250.0,
            )

TimeStepper=CGDycore.TimeStepper()

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
  Exchange = CGDycore.InitExchangeCG(SubGrid,OrdPoly,CellToProc,Proc,ProcNumber,Parallel)
  Output=CGDycore.Output(Topography)
  Global = CGDycore.Global(SubGrid,Model,TimeStepper,Phys,Output,Exchange,OrdPoly+1,nz,NumV,NumTr,())
  Global.Metric=CGDycore.Metric(OrdPoly+1,OrdPolyZ+1,SubGrid.NumFaces,nz)
else
  CellToProc=zeros(0)
  Proc = 0
  ProcNumber = 0
  CGDycore.AddVerticalGrid!(SubGrid,nz,H)
  Exchange = CGDycore.InitExchange(Grid,OrdPoly,CellToProc,Proc,ProcNumber,Parallel)
  Output=CGDycore.Output(Topography)
  Global = CGDycore.Global(Grid,Model,TimeStepper,Phys,Output,Exchange,OrdPoly+1,nz,NumV,NumTr,())
  Global.Metric=CGDycore.Metric(OrdPoly+1,OrdPolyZ+1,Grid.NumFaces,nz)
end  
  (CG,Global)=CGDycore.DiscretizationCG(OrdPoly,OrdPolyZ,CGDycore.JacobiDG3,Global)
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
  U[:,:,Model.ThPos]=CGDycore.Project(CGDycore.fTheta,0.0,CG,Global,Param).*U[:,:,Model.RhoPos]
  if NumTr>0
    U[:,:,Model.RhoVPos+Model.NumV]=CGDycore.Project(CGDycore.fQv,0.0,CG,Global,Param).*U[:,:,Model.RhoPos]
  end   

# Output
  Output.vtkFileName=string("SchaerCart_")
  Output.vtk=0
  Output.Flat=false
  Output.H=H
  Output.cNames = [
    "Rho",
    "u",
    "v",
    "w",
    "Th",
]
  Output.OrdPrint=CG.OrdPoly
  @show "Compute vtkGrid"
  Global.vtkCache = CGDycore.vtkInit(Output.OrdPrint,CGDycore.TransCartX,CG,Global)

  IntMethod="RosenbrockD"
  IntMethod="LinIMEX"
  IntMethod="Rosenbrock"
  IntMethod="RungeKutta"
  IntMethod="IMEX"
  IntMethod="MIS"
  if IntMethod == "Rosenbrock" || IntMethod == "RosenbrockD" || 
    IntMethod == "RosenbrockSSP" || IntMethod == "LinIMEX" || IntMethod == "IMEX"
    dtau = 0.4
  elseif IntMethod == "MIS"  
    dtau = 6.0
    dtauFast = 0.4
  else
    dtau = 0.2
  end
  Global.ROS=CGDycore.RosenbrockMethod("RODAS")
  Global.ROS=CGDycore.RosenbrockMethod("M1HOMME")
  Global.ROS=CGDycore.RosenbrockMethod("SSP-Knoth")
  Global.RK=CGDycore.RungeKuttaMethod("RK3")
  Global.RK=CGDycore.RungeKuttaMethod("RK4")
  Global.IMEX=CGDycore.IMEXMethod("ARS343")
  Global.MIS=CGDycore.MISMethod("MISRK3")
  Global.MIS=CGDycore.MISMethod("MISRKJeb")
  Global.LinIMEX=CGDycore.LinIMEXMethod("ARS343")
  Global.LinIMEX=CGDycore.LinIMEXMethod("M1HOMME")
  Global.LinIMEX=CGDycore.LinIMEXMethod("AR2")

# Simulation period
  time=[0.0];
  SimHours=.5;
  PrintHours=.05;
  nIter=SimHours*3600/dtau;
  PrintInt=PrintHours*3600/dtau;
  PrintStartSec = 0
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
  elseif IntMethod == "IMEX"
    Global.J = CGDycore.JStruct(CG.NumG,nz,Model.NumTr)
    Global.Cache.fV=zeros(size(U))
    Global.Cache.R=zeros(size(U))
    Global.Cache.dZ=zeros(size(U))
    Global.Cache.Y=zeros(size(U[:,:,1:NumV+NumTr])..., Global.IMEX.nStage);
    Global.Cache.Z=zeros(size(U[:,:,1:NumV+NumTr])..., Global.IMEX.nStage);
    Global.Cache.Vn=zeros(size(U))
  elseif IntMethod == "MIS"
    Global.Cache.f=zeros(size(U)..., Global.MIS.nStage)
    Global.Cache.VS=zeros(size(U)..., Global.MIS.nStage - 1)
    Global.J = CGDycore.JStruct(CG.NumG,nz,Model.NumTr)
    Global.Cache.k=zeros(size(U[:,:,1:NumV+NumTr])..., Global.ROS.nStage);
    Global.Cache.fV=zeros(size(U))
    Global.Cache.Vn=zeros(size(U))
    Global.Cache.R=zeros(size(U))
    Global.Cache.dZ=zeros(size(U))
  end

  if Global.Model.Thermo == "TotalEnergy"
    E = zeros(Float64,nz,CG.NumG)
    CGDycore.Energy!(E,U[:,:,Model.ThPos],U[:,:,Model.RhoPos],U[:,:,Model.NumV+1:Model.NumV+Model.NumTr],
      U[:,:,Model.uPos:Model.wPos],CG,Global)
    @views @. U[:,:,Model.ThPos] = E
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
          CGDycore.LinIMEXSchur!(U,dtau,CGDycore.FcnNHCurlVecI!,CGDycore.JacSchur!,CG,Global,Param);
          time[1] += dtau
          if mod(i,PrintInt) == 0 && i >= PrintStartInt
            CGDycore.unstructured_vtkSphere(U,CGDycore.TransCartX,CG,Global,Proc,ProcNumber)
          end
        end
        percent = i/nIter*100
        @info "Iteration: $i took $Δt, $percent% complete"
      end
    end
  elseif IntMethod == "IMEX"
    @time begin
      for i=1:nIter
        Δt = @elapsed begin
          CGDycore.IMEXSchur!(U,dtau,CGDycore.FcnNHCurlExp1DVecI!,CGDycore.FcnNHCurlImp1DGlobalVecI!,CGDycore.JacSchur!,CG,Global,Param);
          time[1] += dtau
          if mod(i,PrintInt) == 0 && i >= PrintStartInt
            CGDycore.unstructured_vtkSphere(U,CGDycore.TransCartX,CG,Global,Proc,ProcNumber)
          end
        end
        percent = i/nIter*100
        @info "Iteration: $i took $Δt, $percent% complete"
      end
    end
  elseif IntMethod == "MIS"
    @time begin
      for i=1:nIter
        Δt = @elapsed begin
          CGDycore.MISSchur!(U,dtau,dtauFast,CGDycore.FcnNHCurlExp3DVecI!,CGDycore.FcnNHCurlImp3DVecI!,CGDycore.JacSchur!,CG,Global,Param);
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
          CGDycore.RungeKuttaExplicit!(U,dtau,CGDycore.FcnNHCurlExp1DVecI!,CGDycore.FcnNHCurlImp1DVecI!,CG,Global,Param)
#         CGDycore.RungeKuttaExplicit!(U,dtau,CGDycore.FcnNHCurlExp3DVecI!,CGDycore.FcnNHCurlImp3DVecI!,CG,Global,Param)
#         CGDycore.RungeKuttaExplicit!(U,dtau,CGDycore.FcnNHCurlVecI!,CG,Global,Param)

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
