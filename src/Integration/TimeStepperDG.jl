mutable struct CacheStruct{FT<:AbstractFloat,
                           FE<:FiniteElements.DGElement,
                           AT3<:AbstractArray,
                           AT4<:AbstractArray,
                           AT5<:AbstractArray}
  U::AT4
  Vn::AT4
  S::AT4
  k::AT5
  fV::AT3
end

function CacheStruct{FT}(backend,DG::FiniteElements.DGElement,M,nz,NumV,NumAux,NumStages) where 
  {FT<:AbstractFloat, FE<:FiniteElements.DGElement}
  NumG = DG.NumG
  NumI = DG.NumI
  U = KernelAbstractions.zeros(backend,FT,M,nz,NumG,NumV+NumAux)
  fV = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV)
  S = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV)
  k = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV,NumStages)
  @views UI = U[:,:,1:NumI,1:NumV]
  return CacheStruct{FT,
                     typeof(DG),
                     typeof(fV),
                     typeof(U),
                     typeof(k)}(
    U,
    UI,
    S,
    k,
    fV,
  )
end

function TimeStepperDG(U,Fcn,Jac,dt,nIter,nPrint,DG,Exchange,Metric,Trans,Phys,Param,Grid,
  Global,ElemType::Grids.ElementType,VelForm)

  backend = get_backend(U)
  FT = eltype(U)
  TimeStepper = Global.TimeStepper
  Output = Global.Output
  Proc = Global.ParallelCom.Proc
  ProcNumber = Global.ParallelCom.ProcNumber
  IntMethod = TimeStepper.IntMethod
  Table = TimeStepper.Table

  if IntMethod == "Rosenbrock" || IntMethod == "RosenbrockSSP" || IntMethod == "RosenbrockAMD"
    TimeStepper.ROS=RosenbrockStruct{FT}(Table)
  elseif IntMethod == "RungeKutta"
    TimeStepper.RK=RungeKuttaMethod{FT}(Table)
  elseif IntMethod == "IMEX"
    TimeStepper.IMEX=IMEXMethod(Table)
  elseif IntMethod == "MIS"
    TimeStepper.MIS=MISMethod(Table)
  elseif IntMethod == "LinIMEX"
    TimeStepper.LinIMEX=LinIMEXMethod(Table)
  end

  # Simulation period
  time=[0.0]
  dtau = FT(TimeStepper.dtau)
  SimDays = TimeStepper.SimDays
  SimHours = TimeStepper.SimHours
  SimMinutes = TimeStepper.SimMinutes
  SimSeconds = TimeStepper.SimSeconds
  SimTime = TimeStepper.SimTime
  PrintDays = Output.PrintDays
  PrintHours = Output.PrintHours
  PrintMinutes = Output.PrintMinutes
  PrintSeconds = Output.PrintSeconds
  PrintTime = Output.PrintTime
  PrintStartTime = Output.PrintStartTime
  StartAverageDays = Output.StartAverageDays
  SimTime = 24*3600*SimDays+3600*SimHours+60*SimMinutes+SimSeconds+SimTime
  nIter=ceil(SimTime/dtau)
  PrintInt=ceil((24*3600*PrintDays+3600*PrintHours+60*PrintMinutes+PrintSeconds+PrintTime)/dtau)
  StartAverageTime = StartAverageDays * 3600 * 24
  if StartAverageTime >= 0 && StartAverageTime < SimTime
    UAver = KernelAbstractions.zeros(backend,FT,size(U))
    @. UAver = 0
    iAv = 1
  end
  PrintStartInt=0

  NumV = Global.Model.NumV
  NumTr = Global.Model.NumTr
  NumAux = Global.Model.NumAux
  ND = Global.Model.NDEDMF
  nz = Global.Grid.nz
  M = DG.OrdPolyZ + 1


  if IntMethod == "Rosenbrock" || IntMethod == "RosenbrockSSP" || IntMethod == "RosenbrockAMD"
    Cache = CacheStruct{FT}(backend,DG,M,nz,NumV,NumAux,TimeStepper.ROS.nStage)
    CacheU = Cache.U
    @views Aux = CacheU[:,:,:,NumV+1:NumV+NumAux]
    @views Un = CacheU[:,:,:,1:NumV]
    @views UI = U[:,:,1:DG.NumI,:]
    @views UnI = Un[:,:,1:DG.NumI,:]
  elseif IntMethod == "RungeKutta"
    TimeStepper.RK=RungeKuttaMethod{FT}(Table)
  elseif IntMethod == "IMEX"
    TimeStepper.IMEX=IMEXMethod(Table)
  elseif IntMethod == "MIS"
    TimeStepper.MIS=MISMethod(Table)
  elseif IntMethod == "LinIMEX"
    TimeStepper.LinIMEX=LinIMEXMethod(Table)
  end

  dz = Metric.dz

  @show TimeStepper.ROS

  Outputs.unstructured_vtkSphere(U,Trans,DG,Metric,Phys,Global,Proc,ProcNumber;Thermo=Aux)
  JCache = DGSEM.JacDGVert{FT}(backend,M,nz,DG.NumI)
  @time begin
    @inbounds for i = 1 : nIter
      Δt = @elapsed begin
        Rosenbrock!(U,dt,Fcn,Jac,DG,Metric,Phys,Cache,JCache,Exchange,
          Global,Param,VelForm)
        if mod(i,nPrint) == 0
          Outputs.unstructured_vtkSphere(U,Trans,DG,Metric,Phys,Global,Proc,ProcNumber;Thermo=Aux)
        end
      end
      percent = i/nIter*100
      if Proc == 1
        @info "Iteration: $i took $Δt, $percent% complete"
      end
    end
  end  
  Outputs.unstructured_vtkSphere(U,Trans,DG,Metric,Phys,Global,Proc,ProcNumber;Thermo=Aux)

#=
  if IntMethod == "Rosenbrock"
    Ros = Integration.RosenbrockStruct{FTB}(Table)
    DGSEM.Rosenbrock(Ros,U,DGSEM.FcnGPUSplit!,dtau,IterTime,nPrint,DG,Exchange,Metric,
      Trans,Phys,Param,Grid,Global,Grid.Type,VelForm)
  elseif IntMethod == "RosenbrockNonConservative"
    Ros = Integration.RosenbrockStruct{FTB}(Table)
    DGSEM.Rosenbrock(Ros,U,DGSEM.FcnGPUNonConservativeSplit!,dtau,IterTime,nPrint,DG,Exchange,Metric,
      Trans,Phys,Param,Grid,Global,Grid.Type)
  elseif IntMethod == "MIS"
    Ros = Integration.RosenbrockStruct{FTB}(Table)
    Mis = DGSEM.MISStruct{FTB}("MISRK4")
  DGSEM.MIS_Method(Ros,Mis,U,DGSEM.FcnGPUSplitSlow!,DGSEM.FcnGPUSplitFast!,dtauSmall,dtau,IterTime,nPrint,DG,Exchange,Metric,Trans,Phys,Param,Grid,Global)
  elseif IntMethod == "RungeKutta"
    DGSEM.RK3(U,DGSEM.FcnGPUSplit!,dtau,IterTime,nPrint,DG,Exchange,Metric,
      Trans,Phys,Grid,Global)
  elseif IntMethod == "RungeKuttaNonConservative"
    DGSEM.RK3(U,DGSEM.FcnGPUNonConservativeSplit!,dtau,IterTime,nPrint,DG,Exchange,Metric,
      Trans,Phys,Grid,Global)
  end
=#

end

#=
function Rosenbrock(ROS,U,Fcn,dtau,nIter,nPrint,DG,Exchange,Metric,Trans,Phys,Param,Grid,
  Global,ElemType::Grids.ElementType,VelForm)

  ROS = Global.TimeStepper.ROS
  nStage = ROS.nStage
  k = Cache.k
  fV = Cache.fV
  Vn = Cache.Vn
  Global.TimeStepper.dtauStage = dt

  fac = dtau * ROS.gamma
  Jac!(Jac,fac,U,DG,Metric,Phys,Cache,Global,Param)
  @inbounds for iStage = 1 : nStage
    @. UnI = UI
    @inbounds for jStage = 1 : iStage-1
      @views @. UnI = UnI + ROS.a[iStage,jStage] * k[:,:,:,:,jStage]
    end
    @views Fcn(k[:,:,:,:,iStage],Un,DG,Metric,Phys,Cache,Exchange,Global,ElemType,VelForm)
    @inbounds for jStage = 1 : iStage - 1
      fac = ROS.c[iStage,jStage] / dtau
      @views @. k[:,:,:,:,iStage] += fac * k[:,:,:,:,jStage]
    end
    fac = dtau * ROS.gamma
    @views Solve!(k[:,:,:,:,iStage],k[:,:,:,:,iStage],Jac,fac,DG,Metric,Global,VelForm)
  end
  @inbounds for iStage = 1 : nStage
    @views @. UI = UI + ROS.m[iStage] * k[:,:,:,:,iStage]
  end
end  

function RosenbrockSparse(ROS,U,Fcn,dtau,nIter,nPrint,DG,Exchange,Metric,Trans,Phys,Param,Grid,Global)

  backend = get_backend(U)
  FTB = eltype(U)
  Proc = Global.ParallelCom.Proc
  ProcNumber = Global.ParallelCom.ProcNumber
  Model = Global.Model
  NumV = Model.NumV
  NumAux = Model.NumAux
  M = size(U,1)
  nz = size(U,2)
  CacheU = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumG,NumV+NumAux)
  @views Aux = CacheU[:,:,:,NV+1:NV+NAUX]
  @views Un = CacheU[:,:,:,1:NumV]
  CacheS = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV)
  @views UI = U[:,:,1:DG.NumI,:]
  @views UnI = Un[:,:,1:DG.NumI,:]
  nStage = ROS.nStage
  k = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV,nStage)
  kS = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV)

  dSdS,dSdM,dMdS,dMdM = InitJacDG(DG,nz,Param)
  N = size(dSdS,1)
  M = DG.OrdPolyZ + 1
  dz = Metric.dz 


  Outputs.unstructured_vtkSphere(U,Trans,DG,Metric,Phys,Global,Proc,ProcNumber;Aux)
  kLoc = zeros(3*N)
  @inbounds for i = 1 : nIter
    if Proc == 1
      @show "Sparse",i,nPrint
    end  
    fac = (1.0 / (dtau * ROS.gamma))
    Jac = JacDG(U,DG,fac,dSdS,dSdM,dMdS,dMdM,Metric.dz,Phys)  
    @inbounds for iStage = 1 : nStage
      @. UnI = UI
      @inbounds for jStage = 1 : iStage-1
        @views @. UnI = UnI + ROS.a[iStage,jStage] * k[:,:,:,:,jStage]
      end
      @views Fcn(k[:,:,:,:,iStage],Un,DG,Model,Metric,Exchange,Grid,CacheU,CacheS,Phys,Global,Grid.Type)
      @inbounds for jStage = 1 : iStage - 1
        fac = ROS.c[iStage,jStage] / dtau
        @views @. k[:,:,:,:,iStage] += fac * k[:,:,:,:,jStage]
      end
      @inbounds for ID = 1 : DG.NumI
        ikLoc = 0
        @inbounds for iv in [1,4,5]
          @inbounds for iz = 1 : nz
            @inbounds for i = 1 : M
              ikLoc += 1
              kLoc[ikLoc] = k[i,iz,ID,iv,iStage]
            end
          end
        end
        ldiv!(Jac[ID],kLoc)
        ikLoc = 0
        @inbounds for iv in [1,4,5]
          @inbounds for iz = 1 : nz
            @inbounds for i = 1 : M
              ikLoc += 1
              k[i,iz,ID,iv,iStage] = kLoc[ikLoc]
            end
          end
        end
      end
      @views @. k[:,:,:,2:3,iStage] *= (dtau * ROS.gamma)
    end
    @inbounds for iStage = 1 : nStage
      @views @. UI = UI + ROS.m[iStage] * k[:,:,:,:,iStage]
    end
    if mod(i,nPrint) == 0 || i == nIter
      if Proc == 1
        @show "Print",i
      end
      Outputs.unstructured_vtkSphere(U,Trans,DG,Metric,Phys,Global,Proc,ProcNumber;Aux)
    end
  end  
end  


function RosenbrockEul(ROS,U,Fcn,dtau,nIter,nPrint,DG,Exchange,Metric,Trans,Phys,Param,Grid,Global)

  backend = get_backend(U)
  FTB = eltype(U)
  Proc = Global.ParallelCom.Proc
  ProcNumber = Global.ParallelCom.ProcNumber
  Model = Global.Model
  NumV = Model.NumV
  NumAux = Model.NumAux
  nz = size(U,2)
  CacheU = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumG,NumV+NumAux)
  CacheS = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV)
  fU = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV)
  Un = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV)
  @views UI = U[:,:,1:DG.NumI,:]
  nStage = ROS.nStage
  k = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV,nStage)

  dSdS,dSdM,dMdS,dMdM = InitJacDG(DG,nz,Param)
  N = size(dSdS,1)
  M = DG.OrdPolyZ + 1


  Outputs.unstructured_vtkSphere(U,Trans,DG,Metric,Phys,Global,Proc,ProcNumber)
  @inbounds for i = 1 : nIter
    if Proc == 1
      @show "RosEul",i,nPrint
    end  
    Jac = JacDG(U,DG,dSdS,dSdM,dMdS,dMdM,Metric.dz,Phys)  
    for ID = 1 : DG.NumI
       Jac[ID] = (1.0 / dtau) * sparse(I,3*N,3*N) - Jac[ID]
    end  
    Fcn(fU,U,DG,Model,Metric,Exchange,Grid,CacheU,CacheS,Phys,Global,Grid.Type)
    for ID = 1 : DG.NumI
      @views fU[:,:,ID,[1,4,5]] .= reshape(Jac[ID] \ reshape(fU[:,:,ID,[1,4,5]],3*N),M,nz,3)
    end  
    @views @. fU[:,:,:,2:3] = fU[:,:,:,2:3,:] * dtau 
    @. UI = UI + fU
    if mod(i,nPrint) == 0 || i == nIter
      if Proc == 1
        @show "Print",i
      end
      Outputs.unstructured_vtkSphere(U,Trans,DG,Metric,Phys,Global,Proc,ProcNumber)
    end
  end  
end  

function RK3(U,Fcn,dtau,nIter,nPrint,DG,Exchange,Metric,Trans,Phys,Grid,Global)

  backend = get_backend(U)
  FTB = eltype(U)
  Proc = Global.ParallelCom.Proc
  ProcNumber = Global.ParallelCom.ProcNumber
  Model = Global.Model
  NumV = Model.NumV
  NumAux = Model.NumAux
  M = size(U,1)
  nz = size(U,2)
  Cache = CacheStructDG{FTB}(backend,DG.NumG,DG.NumI,M,nz,NumV,NumAux)
  CacheU = Cache.U
  @views Un = CacheU[:,:,:,1:NumV]
  @views UI = U[:,:,1:DG.NumI,:]
  @views UnI = Un[:,:,1:DG.NumI,:]
  @views UNew = CacheU[:,:,:,1:NumV]
  CacheS = Cache.S
  FU = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV)
  @views UI = U[:,:,1:DG.NumI,:]
  @views UNewI = UNew[:,:,1:DG.NumI,:]

  ElemType = Grid.Type
  
  Outputs.unstructured_vtkSphere(U,Trans,DG,Metric,Phys,Global,Proc,ProcNumber)
  @inbounds for i = 1 : nIter
      if Proc == 1
        @show i,nPrint
      end  

    @. UNewI = UI 

    Fcn(FU,UNew,DG,Metric,Phys,Cache,Exchange,Global,ElemType)
    fac = FTB(1/3 * dtau)
    @. UNewI = UI + fac * FU

    Fcn(FU,UNew,DG,Metric,Phys,Cache,Exchange,Global,ElemType)
    fac = FTB(1/2 * dtau)
    @. UNewI = UI + fac * FU

    Fcn(FU,UNew,DG,Metric,Phys,Cache,Exchange,Global,ElemType)
    fac = FTB(dtau)
    @. UI = UI + fac * FU

    if mod(i,nPrint) == 0 || i == nIter
      if Proc == 1
        @show "Print",i
      end  
      Outputs.unstructured_vtkSphere(U,Trans,DG,Metric,Phys,Global,Proc,ProcNumber)
    end
  end
end


function RK1(U,Fcn,dtau,nIter,nPrint,DG,Exchange,Metric,Trans,Phys,Grid,Global)

  backend = get_backend(U)
  FTB = eltype(U)
  Proc = Global.ParallelCom.Proc
  ProcNumber = Global.ParallelCom.ProcNumber
  Model = Global.Model
  NumV = Model.NumV
  NumAux = Model.NumAux
  nz = size(U,1)
  CacheU = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumG,NumV+NumAux)
  CacheS = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV)
  FU = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV)
  @views UI = U[:,:,1:DG.NumI,:]
  
  Outputs.unstructured_vtkSphere(U,Trans,DG,Metric,Phys,Global,Proc,ProcNumber)
  @inbounds for i = 1 : nIter
      if Proc == 1
        @show i,nPrint
      end  

    Fcn(FU,U,DG,Model,Metric,Exchange,Grid,CacheU,CacheS,Phys,Global,Grid.Type)
    fac = FTB(dtau)
    @. UI = UI + fac * FU

    if mod(i,nPrint) == 0 || i == nIter
      if Proc == 1
        @show "Print",i
      end  
      Outputs.unstructured_vtkSphere(U,Trans,DG,Metric,Phys,Global,Proc,ProcNumber)
    end
  end
end
=#
