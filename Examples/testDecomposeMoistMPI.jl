using CGDycore
using MPI
using Base

Base.@kwdef struct ParamStruct
  day = 3600.0 * 24.0
  k_a=1.0/(40.0 * day)
  k_f=1.0/day
  k_s=1.0/(4.0 * day)
  DeltaT_y=0.0
  DeltaTh_z=-5.0
  T_equator=315.0
  T_min=200.0
  sigma_b=7.0/10.0
  CE = 0.0044
  CH = 0.0044
  CTr = 0.004
  p_pbl = 85000.0
  p_strato = 10000.0
  T_virt_surf = 290.0
  T_min_ref = 220.0
  H_t = 8.e3
  q_0 = 0.018                # Maximum specific humidity (default: 0.018)
  q_t = 1e-12
  T0E = 310.0
  T0P = 240.0
  B=2.0
  K=3.0
  LapseRate=0.005
  DeltaTS = 29.0
  TSMin = 271.0
  DeltaLat = 26.0 * pi / 180.0
end  
# k_a::Float64
# k_f::Float64
# k_s::Float64
# DeltaT_y::Float64
# DeltaTh_z::Float64
# T_equator::Float64
# T_min::Float64
# sigma_b::Float64
# CE::Float64
# CH::Float64
# CTr::Float64
# p_pbl::Float64
# p_strato::Float64
# T_virt_surf::Float64
# T_min_ref::Float64
# H_t::Float64
# q_0::Float64
# q_t::Float64
# T0E::Float64
# T0P::Float64
#nd

#Base.@kwdef struct ParamStruct1
#  B=2.0
#  K=3.0
#  LapseRate=0.005
#end  
#function ParamStruct()
#  day = 3600.0 * 24.0
#  k_a=1.0/(40.0 * day)
#  k_f=1.0/day
#  k_s=1.0/(4.0 * day)
##  DeltaT_y=0.0
#  DeltaTh_z=-5.0
#  T_equator=315.0
#  T_min=200.0
#  sigma_b=7.0/10.0
#  CE = 0.0044
##  CH = 0.0044
#  CTr = 0.004
#  p_pbl = 85000.0
#  p_strato = 10000.0
#  T_virt_surf = 290.0
#  T_min_ref = 220.0
#  H_t = 8.e3
#  q_0 = 0.018                # Maximum specific humidity (default: 0.018)
#  q_t = 1e-12
#  T0E = 310.0
#  T0P = 240.0
#  return ParamStruct(
#    k_a,
#    k_f,
#    k_s,
##    DeltaT_y,
#    DeltaTh_z,
#    T_equator,
##    T_min,
#    sigma_b,
#    CE,
#    CH,
##    CTr,
#    p_pbl,
#    p_strato,
#    T_virt_surf,
#    T_min_ref,
#    H_t,
##    q_0,
#    q_t,
#    T0E,
#    T0P,
#  )
#end 

MPI.Init()
comm = MPI.COMM_WORLD
Proc = MPI.Comm_rank(comm) + 1
ProcNumber = MPI.Comm_size(comm)
print("$Proc: \n")
print("$ProcNumber: \n")

OrdPoly = 4
nz = 45

OrdPolyZ=1
nPanel = 12
NF = 6 * nPanel * nPanel
NumV = 5
NumTr = 2
Parallel = true


# Physical parameters
Phys=CGDycore.PhysParameters()


#ModelParameters
Param = ParamStruct()
#day=3600.0*24.0
#Param=(T0E=310.0,
#       T0P=240.0,
#       B=2.0,
#       K=3.0,
#       LapseRate=0.005,
#       U0=-0.5,
#       PertR=1.0/6.0,
#       Up=0.0, #1.0,
#       PertExpR=0.1,
#       PertLon=pi/9.0,
#       PertLat=2.0*pi/9.0,
#       PertZ=15000.0,
#       NBr=1.e-2,
#       DeltaT=1,
#       ExpDist=5,
#       T0=300,
#       T_init = 315,
#       lapse_rate = -0.008,
#       Deep=false,
#       k_a=1.0/(40.0 * day),
#       k_f=1.0/day,
#       k_s=1/(4 * day),
#       pert = 0.1,
#       uMax = 0.0,
#       vMax = 0.0,
#       DeltaT_y=0,
#       DeltaTh_z=-5.0,
#       T_equator=315.0,
#       T_min=200.0,
#       sigma_b=7.0/10.0,
#       z_D=20.0e3,
##      Moist
#       q_0 = 0.018,                # Maximum specific humidity (default: 0.018)
#       q_t = 1e-12,
#       T_virt_surf = 290.0,
#       T_min_ref = 220.0,
#       H_t = 8.e3,
##      Boundary layer       
#       CE = 0.0044,
#       CH = 0.0044,
#       p_pbl = 85000.0,
#       p_strato = 10000.0,
##      Surface flux       
#       CTr = 0.004,
#       DeltaTS = 29.0,
#       TSMin = 271.0,
#       DeltaLat = 26.0 * pi / 180.0,
#       RelCloud = 1.0 / 1000.0,       
#       )
Model = CGDycore.Model()
# Initial conditions
  Model.Equation="CompressibleMoist"
  Model.NumV=NumV
  Model.NumTr=NumTr
  Model.Problem="HeldSuarezMoistSphere"
  Model.ProfRho="DecayingTemperatureProfile"
  Model.ProfTheta="DecayingTemperatureProfile"
  Model.ProfVel="Const"
  Model.RhoPos=1
  Model.uPos=2
  Model.vPos=3
  Model.wPos=4
  Model.ThPos=5
  Model.RhoVPos = 1
  Model.RhoCPos = 2
  Model.HorLimit = false
  Model.Source = true
  Model.Upwind = true
  Model.Microphysics = true
  Model.RelCloud = 1.0 / 1000.0       
  Model.Damping = false
  Model.StrideDamp=20000.0
  Model.Relax = 1.0/100.0
  Model.Coriolis=true
  Model.CoriolisType="Sphere"
  Model.VerticalDiffusion = true
  Model.SurfaceFlux = true

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
  sigma = 1.0
  lambda = 3.16
  CGDycore.AddStretchICONVerticalGrid!(SubGrid,nz,H,sigma,lambda)
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
  Model.HyperDCurl=7.e15
  Model.HyperDGrad=7.e15
  Model.HyperDDiv=7.e15

# Output
  Output.OrdPrint=CG.OrdPoly
  vtkGrid=CGDycore.vtkCGGrid(CG,CGDycore.TransSphereX,CGDycore.Topo,Global)


  U = zeros(Float64,nz,CG.NumG,Model.NumV+Model.NumTr)
  U[:,:,Model.RhoPos]=CGDycore.Project(CGDycore.fRho,0.0,CG,Global,Param)
  (U[:,:,Model.uPos],U[:,:,Model.vPos])=CGDycore.ProjectVec(CGDycore.fVel,0.0,CG,Global,Param)
  U[:,:,Model.ThPos]=CGDycore.Project(CGDycore.fTheta,0.0,CG,Global,Param).*U[:,:,Model.RhoPos]
  if NumTr>0
    U[:,:,Model.RhoVPos+Model.NumV]=CGDycore.Project(CGDycore.fQv,0.0,CG,Global,Param).*U[:,:,Model.RhoPos]
  end   

# Output
  Output.vtkFileName=string("HeldSuarezMoist",string(Proc),"_")
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
    "Tr1",
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
    dtau = 120
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
  SimDays=1000
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
  Global.Cache.TSurf=CGDycore.ProjectSurf(CGDycore.fTSurf,0.0,CG,Global,Param)

# Print initial conditions
  @show "Print initial conditions"
  Global.Output.vtk=CGDycore.vtkOutput(U,vtkGrid,CG,Global)

  @show "Choose integration method"
  if IntMethod == "Rosenbrock"
    @time begin
      for i=1:nIter
        Δt = @elapsed begin
          CGDycore.RosenbrockSchur!(U,dtau,CGDycore.FcnNHCurlVecI!,CGDycore.JacSchur!,CG,Global,Param);
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