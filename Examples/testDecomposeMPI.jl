using CGDycore
using MPI

#function testDecompose()
MPI.Init()
comm = MPI.COMM_WORLD
Proc = MPI.Comm_rank(comm) + 1
ProcNumber = MPI.Comm_size(comm)
print("$Proc: \n")
print("$ProcNumber: \n")

OrdPoly = 4
nz = 20

OrdPolyZ=1
nPanel = 24
NF = 6 * nPanel * nPanel
NumV = 5
NumTr = 0
Parallel = true


# Physical parameters
Phys=CGDycore.PhysParameters()

#ModelParameters
day=3600*24
Param=(T0E=310.0,
       T0P=240.0,
       B=2.0,
       K=3.0,
       LapseRate=0.005,
       U0=-0.5,
       PertR=1.0/6.0,
       Up=1.0,
       PertExpR=0.1,
       PertLon=pi/9.0,
       PertLat=2.0*pi/9.0,
       PertZ=15000.0,
       NBr=1.e-2,
       DeltaT=1,
       ExpDist=5,
       T0=300,
       T_init = 315,
       lapse_rate = -0.008,
       Deep=false,
       k_a=1/(40 * day),
       k_f=1/day,
       k_s=1/(4 * day),
       pert = 0.1,
       uMax = 1.0,
       vMax = 0.0,
       DeltaT_y=0,
       DeltaTh_z=-5,
       T_equator=315,
       T_min=200,
       sigma_b=7/10,
       z_D=20.0e3,
#      Boundary layer       
       C_E = 0.0044,
       p_pbl = 85000.0,
       p_strato = 10000.0,
#      Surface flux       
       CTr = 0.004,
#      Moist
       q_0 = 0.018,                # Maximum specific humidity (default: 0.018)
       q_t = 1e-12,
       )
Model = CGDycore.Model(Param)
# Initial conditions
  Model.Equation="Compressible"
  Model.NumV=NumV
  Model.NumTr=NumTr
  Model.Problem="HeldSuarez"
  Model.ProfRho="BaroWaveSphere"
  Model.ProfTheta="BaroWaveSphere"
  Model.ProfVel="BaroWaveSphere"
  Model.RhoPos=1
  Model.uPos=2
  Model.vPos=3
  Model.wPos=4
  Model.ThPos=5
  Model.RhoVPos = 1
  Model.RhoCPos = 2
  Model.HorLimit = false
  Model.Source = false
  Model.Upwind = true
  Model.Microphysics = false
  Model.RelCloud = 0.01
  Model.Damping = false
  Model.StrideDamp=20000.0
  Model.Relax = 1.0/100.0
  Model.Coriolis=true
  Model.CoriolisType="Sphere"

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
  CGDycore.AddVerticalGrid!(SubGrid,nz,H)
  Exchange = CGDycore.InitExchange(SubGrid,OrdPoly,CellToProc,Proc,ProcNumber,Parallel)
  Output=CGDycore.Output(Topography)
  Global = CGDycore.Global(SubGrid,Model,Phys,Output,Exchange,OrdPoly+1,nz,NumV,NumTr,())
  Global.Metric=CGDycore.Metric(OrdPoly+1,OrdPolyZ+1,SubGrid.NumFaces,nz)
else
  CellToProc=zeros(0)
  Proc = 0
  ProcNumber = 0
  CGDycore.AddVerticalGrid!(Grid,nz,H)
  Exchange = CGDycore.InitExchange(Grid,OrdPoly,CellToProc,Proc,ProcNumber,Parallel)
  Output=CGDycore.Output(Topography)
  Global = CGDycore.Global(Grid,Model,Phys,Output,Exchange,OrdPoly+1,nz,NumV,NumTr,())
  Global.Metric=CGDycore.Metric(OrdPoly+1,OrdPolyZ+1,Grid.NumFaces,nz)
end  
  Grid = nothing
  (CG,Global)=CGDycore.Discretization(OrdPoly,OrdPolyZ,CGDycore.JacobiSphere3,Global)
  Model.HyperVisc=true
  Model.HyperDCurl=2.e17/4/2^4 #1.e14*(dx/LRef)^3.2;
  Model.HyperDGrad=2.e17/4/2^4
  Model.HyperDDiv=2.e17/4/2^4 # Scalars

# Output
  Output.OrdPrint=CG.OrdPoly
  vtkGrid=CGDycore.vtkCGGrid(CG,CGDycore.TransSphereX,CGDycore.Topo,Global)


  U = zeros(Float64,nz,CG.NumG,Model.NumV+Model.NumTr)
  U[:,:,Model.RhoPos]=CGDycore.Project(CGDycore.fRho,0.0,CG,Global)
  (U[:,:,Model.uPos],U[:,:,Model.vPos])=CGDycore.ProjectVec(CGDycore.fVel,0.0,CG,Global)
  U[:,:,Model.ThPos]=CGDycore.Project(CGDycore.fTheta,0.0,CG,Global).*U[:,:,Model.RhoPos]
  if NumTr>0
    U[:,:,Model.RhoVPos+Model.NumV]=CGDycore.Project(CGDycore.fQv,0.0,CG,Global).*U[:,:,Model.RhoPos]
  end   

# Output
  Output.vtkFileName=string("BaroWaveMoist",string(Proc),"_")
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
]
  Output.OrdPrint=CG.OrdPoly
  @show "Compute vtkGrid"
  vtkGrid=CGDycore.vtkCGGrid(CG,CGDycore.TransSphereX,CGDycore.Topo,Global)

  IntMethod="RungeKutta"
  IntMethod="RosenbrockD"
  IntMethod="LinIMEX"
  IntMethod="RungeKutta"
  IntMethod="Rosenbrock"
  if IntMethod == "Rosenbrock" || IntMethod == "RosenbrockD" || IntMethod == "RosenbrockSSP" || IntMethod == "LinIMEX"
    dtau = 50
  else
    dtau=3
  end
  Global.ROS=CGDycore.RosenbrockMethod("RODAS")
  Global.ROS=CGDycore.RosenbrockMethod("M1HOMME")
  Global.ROS=CGDycore.RosenbrockMethod("SSP-Knoth")
  Global.RK=CGDycore.RungeKuttaMethod("RK4")
  Global.LinIMEX=CGDycore.LinIMEXMethod("ARS343")
  Global.LinIMEX=CGDycore.LinIMEXMethod("AR2")
  Global.LinIMEX=CGDycore.LinIMEXMethod("M1HOMME")

# Simulation period
  time=[0.0]
  SimDays=10
  PrintDay=10
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
  end

# Boundary values
 if Model.NumTr > 0
   @views @. Global.Cache.cTrS[:,:,:,1] = 0.0
   @views @. Global.Cache.cTrS[:,:,1:20,2] = 0.0
 end

# Print initial conditions
  @show "Print initial conditions"
  Global.Output.vtk=CGDycore.vtkOutput(U,vtkGrid,CG,Global)

  @show "Choose integration method"
  if IntMethod == "Rosenbrock"
    @time begin
      for i=1:nIter
        Δt = @elapsed begin
          CGDycore.RosenbrockSchur!(U,dtau,CGDycore.FcnNHCurlVec!,CGDycore.JacSchur!,CG,Global);
          time[1] += dtau
          if mod(i,PrintInt) == 0 && i >= PrintStartInt
            Global.Output.vtk=CGDycore.vtkOutput(U,vtkGrid,CG,Global)
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
          CGDycore.LinIMEXSchur!(U,dtau,CGDycore.FcnNHCurlVec!,CGDycore.JacSchur!,CG,Global);
          time[1] += dtau
          if mod(i,PrintInt) == 0 && i >= PrintStartInt
            Global.Output.vtk=CGDycore.vtkOutput(U,vtkGrid,CG,Global)
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