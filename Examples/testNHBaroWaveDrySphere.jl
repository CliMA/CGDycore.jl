using CGDycore
using MPI
using Base

Base.@kwdef struct ParamStruct
  T0E=310.0
  T0P=240.0
  B=2.0
  K=3.0
  LapseRate=0.005
  U0=-0.5
  PertR=1.0/6.0
  Up=1.0
  PertExpR=0.1
  PertLon=pi/9.0
  PertLat=2.0*pi/9.0
  PertZ=15000.0
  NBr=1.e-2
  DeltaT=1
  ExpDist=5
  T0=300
  T_init = 315
  lapse_rate = -0.008
  Deep=false
  pert = 0.1
  uMax = 1.0
  vMax = 0.0
  DeltaT_y=0
  DeltaTh_z=-5
  T_equator=315
  T_min=200
  sigma_b=7/10
  z_D=20.0e3
  #      Moist
  q_0 = 0.018                # Maximum specific humidity (default: 0.018)
  q_t = 1.0e-12
  Stretch = true
end  

MPI.Init()
comm = MPI.COMM_WORLD
Proc = MPI.Comm_rank(comm) + 1
ProcNumber = MPI.Comm_size(comm)
print("$Proc: \n")
print("$ProcNumber: \n")

OrdPoly = 3
nz = 45

OrdPolyZ=1
nPanel = 16
NF = 6 * nPanel * nPanel
NumV = 5
NumTr = 0
Parallel = true


# Physical parameters
Phys=CGDycore.PhysParameters()


#ModelParameters
Param = ParamStruct()
Model = CGDycore.Model()
# Initial conditions
  Model.Equation="Compressible"
  Model.NumV=NumV
  Model.NumTr=NumTr
  Model.Problem="BaroWaveSphere"
  Model.ProfRho="BaroWaveSphere"
  Model.ProfTheta="BaroWaveSphere"
  Model.ProfVel="BaroWaveSphere"
  Model.RhoPos=1
  Model.uPos=2
  Model.vPos=3
  Model.wPos=4
  Model.ThPos=5
  Model.HorLimit = false
  Model.Source = false
  Model.Upwind = true
  Model.Damping = false
  Model.StrideDamp=20000.0
  Model.Relax = 1.0/100.0
  Model.Coriolis=true
  Model.CoriolisType="Sphere"
  Model.VerticalDiffusion = false
  Model.Thermo = "" #"InternalEnergy" #"" #"TotalEnergy"

# Grid
H = 30000.0
#H = 45000.0
Topography=(TopoS="",H=H,Rad=Phys.RadEarth)




Grid=CGDycore.Grid(nz,Topography)
Grid=CGDycore.CubedGrid(nPanel,CGDycore.OrientFaceSphere,Phys.RadEarth,Grid)
P0Sph = [   0.0,-0.5*pi,Phys.RadEarth]
P1Sph = [2.0*pi, 0.5*pi,Phys.RadEarth]
CGDycore.HilbertFaceSphere!(Grid,P0Sph,P1Sph)
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
  sigma = 1.0
  lambda = 3.16
  CGDycore.AddStretchICONVerticalGrid!(Grid,nz,H,sigma,lambda)
  Exchange = CGDycore.InitExchange(Grid,OrdPoly,CellToProc,Proc,ProcNumber,Parallel)
  Output=CGDycore.Output(Topography)
  Global = CGDycore.Global(Grid,Model,Phys,Output,Exchange,OrdPoly+1,nz,NumV,NumTr,())
  Global.Metric=CGDycore.Metric(OrdPoly+1,OrdPolyZ+1,Grid.NumFaces,nz)
end  
  (CG,Global)=CGDycore.Discretization(OrdPoly,OrdPolyZ,CGDycore.JacobiSphere3,Global)
  Model.HyperVisc=true
  Model.HyperDCurl=1.e14 #7.e15
  Model.HyperDGrad=1.e14 #7.e15
  Model.HyperDDiv=1.e14 #7.e15

# Output
  Output.OrdPrint=CG.OrdPoly

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
  Output.vtkFileName=string("BaroWaveDry_")
  Output.vtk=0
  Output.Flat=true
  Output.nPanel=nPanel
  Output.RadPrint=H
  Output.H=H
  Output.cNames = [
    "Rho",
    "u",
    "v",
    "w",
    "Th",
    "Pres",
]
  Output.OrdPrint=CG.OrdPoly
  Global.vtkCache = CGDycore.vtkInit(Output.OrdPrint,CGDycore.TransSphereX,CG,Global)

  IntMethod="RungeKutta"
  IntMethod="RosenbrockD"
  IntMethod="Rosenbrock"
  IntMethod="LinIMEX"
  IntMethod="IMEX"
  IntMethod="MIS"
  if IntMethod == "Rosenbrock" || IntMethod == "RosenbrockD" || IntMethod == "RosenbrockSSP" || IntMethod == "LinIMEX" || IntMethod == "IMEX"
    dtau = 400
  elseif IntMethod == "MIS" 
    dtau = 1200.0
    dtauFast =  200.0   
  else
    dtau=3
  end
  Global.ROS=CGDycore.RosenbrockMethod("RODAS")
  Global.ROS=CGDycore.RosenbrockMethod("M1HOMME")
  Global.ROS=CGDycore.RosenbrockMethod("SSP-Knoth")
  Global.RK=CGDycore.RungeKuttaMethod("RK4")
  Global.IMEX=CGDycore.IMEXMethod("ARS343")
  Global.MIS=CGDycore.MISMethod("MISRK3")
  Global.MIS=CGDycore.MISMethod("MISRKJeb")
  Global.LinIMEX=CGDycore.LinIMEXMethod("AR2")
  Global.LinIMEX=CGDycore.LinIMEXMethod("ARS343")
  Global.LinIMEX=CGDycore.LinIMEXMethod("M1HOMME")

# Simulation period
  time=[0.0]
  SimDays=10
  PrintDay=.5
  PrintStartDay = 0
  nIter=ceil(24*3600*SimDays/dtau)
  PrintInt=ceil(24*3600*PrintDay/dtau)
  PrintStartInt=ceil(24*3600*PrintStartDay/dtau)
  
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
  elseif IntMethod == "MIS"
    Global.Cache.f=zeros(size(U)..., Global.MIS.nStage)
    Global.Cache.VS=zeros(size(U)..., Global.MIS.nStage - 1)
    Global.J = CGDycore.JStruct(CG.NumG,nz,Model.NumTr)
    Global.Cache.k=zeros(size(U[:,:,1:NumV+NumTr])..., Global.ROS.nStage);
    Global.Cache.fV=zeros(size(U))
    Global.Cache.Vn=zeros(size(U))
    Global.Cache.R=zeros(size(U))
    Global.Cache.dZ=zeros(size(U))  
  elseif IntMethod == "IMEX"
    Global.J = CGDycore.JStruct(CG.NumG,nz,Model.NumTr)
    Global.Cache.fV=zeros(size(U))
    Global.Cache.R=zeros(size(U))
    Global.Cache.dZ=zeros(size(U))
    Global.Cache.Y=zeros(size(U[:,:,1:NumV+NumTr])..., Global.IMEX.nStage);
    Global.Cache.Z=zeros(size(U[:,:,1:NumV+NumTr])..., Global.IMEX.nStage);
    Global.Cache.Vn=zeros(size(U))  
  end

 # Boundary values
  if Model.SurfaceFlux
    Global.Cache.TSurf=CGDycore.ProjectSurf(CGDycore.fTSurf,0.0,CG,Global,Param)
  end


# Print initial conditions
  @show "Print initial conditions"
  CGDycore.unstructured_vtkSphere(U,CGDycore.TransSphereX,CG,Global,Proc,ProcNumber)

  @show "Choose integration method"
  if IntMethod == "Rosenbrock"
    @time begin
      for i=1:nIter
        Δt = @elapsed begin
          CGDycore.RosenbrockSchur!(U,dtau,CGDycore.FcnNHCurlVecI!,CGDycore.JacSchur!,CG,Global,Param);
          time[1] += dtau
          if mod(i,PrintInt) == 0 && i >= PrintStartInt
            CGDycore.unstructured_vtkSphere(U,CGDycore.TransSphereX,CG,Global,Proc,ProcNumber)
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
            Global.Output.vtk=CGDycore.vtkOutput(U,vtkGrid,CG,Global)
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
            Global.Output.vtk=CGDycore.vtkOutput(U,vtkGrid,CG,Global)
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
            CGDycore.unstructured_vtkSphere(U,CGDycore.TransSphereX,CG,Global,Proc,ProcNumber)
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
            CGDycore.unstructured_vtkSphere(U,CGDycore.TransSphereX,CG,Global,Proc,ProcNumber)
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
            CGDycore.unstructured_vtkSphere(U,CGDycore.TransSphereX,CG,Global,Proc,ProcNumber)
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
          @time CGDycore.RungeKuttaExplicit!(U,dtau,CGDycore.FcnNHCurlVec!,CG,Global)

          time[1] += dtau
          if mod(i,PrintInt)==0 && i >= PrintStartInt
            Global.Output.vtk=CGDycore.vtkOutput(U,vtkGrid,CG,Global)
          end
        end
        percent = i/nIter*100
        @info "Iteration: $i took $Δt, $percent% complete"
      end
    end
  else
    error("Bad IntMethod")
  end
