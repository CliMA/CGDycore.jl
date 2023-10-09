using KernelAbstractions
function TimeStepper!(U,Fcn!,FcnPrepare!,Jac!,Trans,CG,Metric,Phys,Global,Param,DiscType)  
  backend = get_backend(U)
  FT = eltype(U)
  TimeStepper = Global.TimeStepper
  Output = Global.Output
  Proc = Global.ParallelCom.Proc
  ProcNumber = Global.ParallelCom.ProcNumber
  IntMethod = TimeStepper.IntMethod
  Table = TimeStepper.Table

  if IntMethod == "Rosenbrock" || IntMethod == "RosenbrockSSP"
    TimeStepper.ROS=RosenbrockStruct{FT}(backend,Table)  
  elseif IntMethod == "RungeKutta"  
    TimeStepper.RK=RungeKuttaMethod(Table)
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
    UAver = similar(U)
    @. UAver = 0.0
    iAv = 1
  end  
  PrintStartInt=0
  Output.OrdPrint=CG.OrdPoly

  NumV = Global.Model.NumV
  NumTr = Global.Model.NumTr
  nz = Global.Grid.nz
  NumG = CG.NumG
  Cache=CacheStruct{FT}(backend,CG.DoF,Global.Grid.NumFaces,Global.Grid.NumGhostFaces,NumG,nz,
    NumV,NumTr)

  if IntMethod == "Rosenbrock" || IntMethod == "RosenbrockD"
    JCache = JStruct{FT}(backend,NumG,nz,NumTr)
    Cache.k = KernelAbstractions.zeros(backend,FT,size(U[:,:,1:NumV+NumTr])..., TimeStepper.ROS.nStage);
    Cache.fV = KernelAbstractions.zeros(backend,FT,size(U))
    Cache.Vn = similar(U)
  elseif IntMethod == "RosenbrockSSP"
    Global.J = JStruct(NumG,nz,NumTr)
    Cache.k=zeros(size(U[:,:,1:NumV])..., TimeStepper.ROS.nStage);
    Cache.fV=zeros(size(U))
    Cache.VS=zeros(size(U[:,:,NumV+1:end])..., TimeStepper.ROS.nStage+1);
    Cache.fS=zeros(size(U[:,:,NumV+1:end])..., TimeStepper.ROS.nStage);
    Cache.fRhoS=zeros(size(U[:,:,1])..., TimeStepper.ROS.nStage);
    Cache.RhoS=zeros(size(U[:,:,1])..., TimeStepper.ROS.nStage+1);
  elseif IntMethod == "LinIMEX"
    Global.J = JStruct(NumG,nz,NumTr)
    Cache.Ymyn=zeros(size(U[:,:,1:NumV+NumTr])..., TimeStepper.LinIMEX.nStage-1);
    Cache.fV=zeros(size(U))
    Cache.f=zeros(size(U)..., TimeStepper.LinIMEX.nStage)
    Cache.fV=zeros(size(U))
    Cache.Vn=zeros(size(U))
  elseif IntMethod == "RungeKutta"
    Cache.f=KernelAbstractions.zeros(backend,FT,size(U)..., TimeStepper.RK.nStage)
  elseif IntMethod == "MIS"
    Cache.f=zeros(size(U)..., TimeStepper.MIS.nStage)
    Cache.VS=zeros(size(U)..., TimeStepper.MIS.nStage - 1)
    Global.J = JStruct(NumG,nz,NumTr)
    Cache.k=zeros(size(U[:,:,1:NumV+NumTr])..., TimeStepper.ROS.nStage);
    Cache.fV=zeros(size(U))
    Cache.Vn=zeros(size(U))
    Cache.R=zeros(size(U))
    Cache.dZ=zeros(size(U))  
  elseif IntMethod == "IMEX"
    Global.J = JStruct(NumG,nz,NumTr)
    Cache.fV=zeros(size(U))
    Cache.R=zeros(size(U))
    Cache.dZ=zeros(size(U))
    Cache.Y=zeros(size(U[:,:,1:NumV+NumTr])..., TimeStepper.IMEX.nStage);
    Cache.Z=zeros(size(U[:,:,1:NumV+NumTr])..., TimeStepper.IMEX.nStage);
    Cache.Vn=zeros(size(U))  
  end

 # Boundary values
  if Global.Model.SurfaceFlux
    Cache.TSurf=ProjectSurf(fTSurf,0.0,CG,Global,Param)
  end


# Print initial conditions
  unstructured_vtkSphere(U,Trans,CG,Metric,Cache,Global,Proc,ProcNumber)

  if IntMethod == "Rosenbrock"
    @time begin
      for i=1:nIter
        Δt = @elapsed begin
          RosenbrockSchur!(U,dtau,Fcn!,FcnPrepare!,Jac!,CG,Metric,Phys,Cache,JCache,Global,Param,DiscType);
          time[1] += dtau
          if mod(i,PrintInt) == 0 && time[1] >= PrintStartTime
            unstructured_vtkSphere(U,Trans,CG,Metric,Cache,Global,Proc,ProcNumber)
          end
          if time[1] >= StartAverageTime && StartAverageTime >= 0.0
            AverageInTime!(UAver,U,iAv)
            iAv += 1
          end  
        end
        percent = i/nIter*100
        if Proc == 1
          @info "Iteration: $i took $Δt, $percent% complete"
        end  
      end
    end
  elseif IntMethod == "RosenbrockD"
    @time begin
      for i=1:nIter
        Δt = @elapsed begin
          RosenbrockDSchur!(U,dtau,FcnNHCurlVec!,JacSchur!,CG,Global);
          time[1] += dtau
          if mod(i,PrintInt) == 0 && time[1] >= PrintStartTime
            unstructured_vtkSphere(U,Trans,CG,Global,Proc,ProcNumber)
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
          RosenbrockSchurSSP!(U,dtau,FcnNHCurlVec!,JacSchur!,CG,Global);
          time[1] += dtau
          if mod(i,PrintInt) == 0 && time[1] >= PrintStartTime
            unstructured_vtkSphere(U,Trans,CG,Global,Proc,ProcNumber)
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
          LinIMEXSchur!(U,dtau,FcnNHCurlVecI!,JacSchur!,CG,Global,Param);
          time[1] += dtau
          if mod(i,PrintInt) == 0 && time[1] >= PrintStartTime
            unstructured_vtkSphere(U,TransSphereX,CG,Global,Proc,ProcNumber)
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
          IMEXSchur!(U,dtau,FcnNHCurlExp1DVecI!,FcnNHCurlImp1DGlobalVecI!,JacSchur!,CG,Global,Param);
          time[1] += dtau
          if mod(i,PrintInt) == 0 && time[1] >= PrintStartTime
            unstructured_vtkSphere(U,TransSphereX,CG,Global,Proc,ProcNumber)
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
          MISSchur!(U,dtau,dtauFast,FcnNHCurlExp3DVecI!,FcnNHCurlImp3DVecI!,JacSchur!,CG,Global,Param);
          time[1] += dtau
          if mod(i,PrintInt) == 0 && time[1] >= PrintStartTime
            unstructured_vtkSphere(U,Trans,CG,Global,Proc,ProcNumber)
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
          @time RungeKuttaExplicit!(U,dtau,Fcn!,FcnPrepare!,CG,Metric,Phys,Cache,Global,Param,DiscType)
          time[1] += dtau
          if mod(i,PrintInt) == 0 && time[1] >= PrintStartTime
            unstructured_vtkSphere(U,Trans,CG,Metric,Cache,Global,Proc,ProcNumber)
          end
        end
        percent = i/nIter*100
        @info "Iteration: $i took $Δt, $percent% complete"
      end
    end
  else
    error("Bad IntMethod")
  end
  if StartAverageTime >= 0 && StartAverageTime < SimTime
    unstructured_vtkSphere(UAver,Trans,CG,Metric,Cache,Global,Proc,ProcNumber)
  end  
end  

function TimeStepperAdvection!(U,Trans,CG,Metric,Phys,Global,Param)  
  TimeStepper = Global.TimeStepper
  Output = Global.Output
  Proc = Global.ParallelCom.Proc
  ProcNumber = Global.ParallelCom.ProcNumber
  IntMethod = TimeStepper.IntMethod
  Table = TimeStepper.Table

  if IntMethod == "Rosenbrock" || IntMethod == "RosenbrockSSP"
    TimeStepper.ROS=RosenbrockMethod(Table)  
  elseif IntMethod == "RungeKutta"  
    TimeStepper.RK=RungeKuttaMethod(Table)
  elseif IntMethod == "SSPRungeKutta"  
    TimeStepper.SSP=SSPRungeKuttaMethod(Table)
  end

# Simulation period
  time=[0.0]
  dtau = TimeStepper.dtau
  SimDays = TimeStepper.SimDays
  SimHours = TimeStepper.SimHours
  SimMinutes = TimeStepper.SimMinutes
  SimSeconds = TimeStepper.SimSeconds
  SimTime = TimeStepper.SimTime
  PrintDays = Output.PrintDays
  PrintHours = Output.PrintHours
  PrintSeconds = Output.PrintSeconds
  PrintTime = Output.PrintTime
  PrintStartTime = Output.PrintStartTime
  nIter=ceil((24*3600*SimDays+3600*SimHours+60*SimMinutes+SimSeconds+SimTime)/dtau)
  PrintInt=ceil((24*3600*PrintDays+3600*PrintHours+PrintSeconds+PrintTime)/dtau)
  PrintStartInt=0
  Output.OrdPrint=CG.OrdPoly

  NumV = Global.Model.NumV
  NumTr = Global.Model.NumTr
  nz = Global.Grid.nz
  NumG = CG.NumG
  FT = eltype(U)
  backend = get_backend(U)
  Cache=CacheStruct{FT}(backend,CG.OrdPoly+1,Global.Grid.NumFaces,Global.Grid.NumGhostFaces,NumG,nz,
    NumV,NumTr)

  if IntMethod == "Rosenbrock" || IntMethod == "RosenbrockD"
    Global.J = JStruct(NumG,nz,NumTr)
    Cache.k=zeros(size(U[:,:,1:NumV+NumTr])..., TimeStepper.ROS.nStage);
    Cache.fV=zeros(size(U))
    Cache.Vn=zeros(size(U))
  elseif IntMethod == "RungeKutta"
    Cache.f=zeros(size(U)..., TimeStepper.RK.nStage)
  elseif IntMethod == "SSPRungeKutta"
    Cache.fS=zeros(nz,NumG,NumTr,TimeStepper.SSP.nStage)
    Cache.VS=zeros(nz,NumG,NumTr,TimeStepper.SSP.nStage+1)
    Cache.fV=zeros(size(U))
    Cache.RhoS=zeros(size(U[:,:,1])..., TimeStepper.SSP.nStage+1)
    Cache.fRhoS=zeros(size(U[:,:,1])..., TimeStepper.SSP.nStage)
  end


# Print initial conditions
# unstructured_vtkSphere(U,TransSphereX,CG,Global,Proc,ProcNumber)

  if IntMethod == "Rosenbrock"
    @time begin
      for i=1:nIter
        Δt = @elapsed begin
          RosenbrockSchur!(U,dtau,FcnTracer!,JacSchur!,CG,Global,Param);
          time[1] += dtau
          if mod(i,PrintInt) == 0 && time[1] >= PrintStartTime
            unstructured_vtkSphere(U,Trans,CG,Global,Proc,ProcNumber)
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
          RungeKuttaExplicit!(time[1],U,dtau,FcnTracer!,CG,Global,Param)
          time[1] += dtau
          if mod(i,PrintInt) == 0 && time[1] >= PrintStartTime
            unstructured_vtkSphere(U,Trans,CG,Global,Proc,ProcNumber)
          end
        end
        percent = i/nIter*100
        @info "Iteration: $i took $Δt, $percent% complete"
      end
    end
  elseif IntMethod == "SSPRungeKutta"
    @time begin
      for i = 1 : nIter
        Δt = @elapsed begin
          SSPRungeKutta!(time[1],U,dtau,FcnTracer!,CG,Metric,Cache,Global,Param)
          time[1] += dtau
          if mod(i,PrintInt) == 0 && time[1] >= PrintStartTime
            unstructured_vtkSphere(U,Trans,CG,Metric,Global,Proc,ProcNumber)
          end
        end
        percent = i/nIter*100
        @info "Iteration: $i took $Δt, $percent% complete"
      end
    end
  else
    error("Bad IntMethod")
  end
end  
function TimeStepperAdvectionConv!(U,Trans,CG,Metric,Global,Param)  
  TimeStepper = Global.TimeStepper
  Output = Global.Output
  Proc = Global.ParallelCom.Proc
  ProcNumber = Global.ParallelCom.ProcNumber
  IntMethod = TimeStepper.IntMethod
  Table = TimeStepper.Table

  if IntMethod == "Rosenbrock" || IntMethod == "RosenbrockSSP"
    TimeStepper.ROS=RosenbrockMethod(Table)  
  elseif IntMethod == "RungeKutta"  
    TimeStepper.RK=RungeKuttaMethod(Table)
  elseif IntMethod == "SSPRungeKutta"  
    TimeStepper.SSP=SSPRungeKuttaMethod(Table)
  end

# Simulation period
  time=[0.0]
  dtau = TimeStepper.dtau
  SimDays = TimeStepper.SimDays
  SimHours = TimeStepper.SimHours
  SimMinutes = TimeStepper.SimMinutes
  SimSeconds = TimeStepper.SimSeconds
  SimTime = TimeStepper.SimTime
  PrintDays = Output.PrintDays
  PrintHours = Output.PrintHours
  PrintSeconds = Output.PrintSeconds
  PrintTime = Output.PrintTime
  PrintStartTime = Output.PrintStartTime
  nIter=ceil((24*3600*SimDays+3600*SimHours+60*SimMinutes+SimSeconds+SimTime)/dtau)
  PrintInt=ceil((24*3600*PrintDays+3600*PrintHours+PrintSeconds+PrintTime)/dtau)
  PrintStartInt=0
  Output.OrdPrint=CG.OrdPoly

  NumV = Global.Model.NumV
  NumTr = Global.Model.NumTr
  nz = Global.Grid.nz
  NumG = CG.NumG
  Global.Cache=CacheCreate(CG.OrdPoly+1,Global.Grid.NumFaces,Global.Grid.NumGhostFaces,NumG,nz,
    NumV,NumTr)

  if IntMethod == "Rosenbrock" || IntMethod == "RosenbrockD"
    Global.J = JStruct(NumG,nz,NumTr)
    Global.Cache.k=zeros(size(U[:,:,1:NumV+NumTr])..., TimeStepper.ROS.nStage);
    Global.Cache.fV=zeros(size(U))
    Global.Cache.Vn=zeros(size(U))
  elseif IntMethod == "RungeKutta"
    Global.Cache.f=zeros(size(U)..., TimeStepper.RK.nStage)
  elseif IntMethod == "SSPRungeKutta"
    Global.Cache.fS=KernelAbstractions.zeros(backend,FT,nz,NumG,NumTr,TimeStepper.SSP.nStage)
    Global.Cache.VS=KernelAbstractions.zeros(backend,FT,nz,NumG,NumTr,TimeStepper.SSP.nStage+1)
    Global.Cache.fV=KernelAbstractions.zeros(backend,FT,size(U))
    Global.Cache.RhoS=KernelAbstractions.zeros(backend,FT,size(U[:,:,1])..., TimeStepper.SSP.nStage+1)
    Global.Cache.fRhoS=KernelAbstractions.zeros(backend,FT,size(U[:,:,1])..., TimeStepper.SSP.nStage)
  end


# Print initial conditions
  unstructured_vtkSphere(U,TransSphereX,CG,Global,Proc,ProcNumber)

  if IntMethod == "Rosenbrock"
    @time begin
      for i=1:nIter
        Δt = @elapsed begin
          RosenbrockSchur!(U,dtau,FcnTracer!,JacSchur!,CG,Global,Param);
          time[1] += dtau
          if mod(i,PrintInt) == 0 && time[1] >= PrintStartTime
            unstructured_vtkSphere(U,Trans,CG,Global,Proc,ProcNumber)
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
          RungeKuttaExplicit!(time[1],U,dtau,FcnTracer!,CG,Global,Param)
          time[1] += dtau
          if mod(i,PrintInt)==0 && i >= PrintStartInt
            unstructured_vtkSphere(U,Trans,CG,Global,Proc,ProcNumber)
          end
        end
        percent = i/nIter*100
        @info "Iteration: $i took $Δt, $percent% complete"
      end
    end
  elseif IntMethod == "SSPRungeKutta"
    @time begin
      for i=1:nIter
        Δt = @elapsed begin
          SSPRungeKutta!(time[1],U,dtau,FcnTracerConv!,Metric,CG,Global,Param)
          time[1] += dtau
          if (mod(i,PrintInt) == 0 && time[1] >= PrintStartTime) || i == nIter
            unstructured_vtkSphere(U,Trans,CG,Global,Proc,ProcNumber)
          end
        end
        percent = i/nIter*100
        @info "Iteration: $i took $Δt, $percent% complete"
      end
    end
  else
    error("Bad IntMethod")
  end
end  

function TimeStepperGPUAdvection!(U,Trans,CG,Metric,Phys,Global,Param,Profile)  
  TimeStepper = Global.TimeStepper
  Output = Global.Output
  Proc = Global.ParallelCom.Proc
  ProcNumber = Global.ParallelCom.ProcNumber
  IntMethod = TimeStepper.IntMethod
  Table = TimeStepper.Table

  if IntMethod == "SSPRungeKutta"  
    TimeStepper.SSP=SSPRungeKuttaMethod(Table)
  end

# Simulation period
  time=[0.0]
  dtau = TimeStepper.dtau
  SimDays = TimeStepper.SimDays
  SimHours = TimeStepper.SimHours
  SimMinutes = TimeStepper.SimMinutes
  SimSeconds = TimeStepper.SimSeconds
  SimTime = TimeStepper.SimTime
  PrintDays = Output.PrintDays
  PrintHours = Output.PrintHours
  PrintSeconds = Output.PrintSeconds
  PrintTime = Output.PrintTime
  PrintStartTime = Output.PrintStartTime
  nIter=ceil((24*3600*SimDays+3600*SimHours+60*SimMinutes+SimSeconds+SimTime)/dtau)
  PrintInt=ceil((24*3600*PrintDays+3600*PrintHours+PrintSeconds+PrintTime)/dtau)
  PrintStartInt=0
  Output.OrdPrint=CG.OrdPoly

  NumV = Global.Model.NumV
  NumTr = Global.Model.NumTr
  nz = Global.Grid.nz
  NumG = CG.NumG
  FT = eltype(U)
  backend = get_backend(U)
  Cache=CacheStruct{FT}(backend,CG.OrdPoly+1,Global.Grid.NumFaces,Global.Grid.NumGhostFaces,NumG,nz,
    NumV,NumTr)

  if IntMethod == "SSPRungeKutta"
    Cache.fS=KernelAbstractions.zeros(backend,FT,nz,NumG,NumTr,TimeStepper.SSP.nStage)
    Cache.VS=KernelAbstractions.zeros(backend,FT,nz,NumG,NumTr,TimeStepper.SSP.nStage+1)
    Cache.fV=KernelAbstractions.zeros(backend,FT,size(U))
    Cache.RhoS=KernelAbstractions.zeros(backend,FT,size(U[:,:,1])..., TimeStepper.SSP.nStage+1)
    Cache.fRhoS=KernelAbstractions.zeros(backend,FT,size(U[:,:,1])..., TimeStepper.SSP.nStage)
  end


# Print initial conditions
  unstructured_vtkSphere(U,TransSphereX,CG,Metric,Cache,Global,Proc,ProcNumber)

  if IntMethod == "SSPRungeKutta"
    @time begin
      for i = 1 : nIter
        Δt = @elapsed begin
          SSPRungeKutta!(time[1],U,dtau,FcnAdvectionGPU!,CG,Metric,Phys,Cache,Global,Param,Profile)
          time[1] += dtau
          if mod(i,PrintInt) == 0 && time[1] >= PrintStartTime
            unstructured_vtkSphere(U,Trans,CG,Metric,Cache,Global,Proc,ProcNumber)
          end
        end 
        percent = i/nIter*100
        @info "Iteration: $i took $Δt, $percent% complete"
      end
    end
  else
    error("Bad IntMethod")
  end
  unstructured_vtkSphere(U,Trans,CG,Metric,Cache,Global,Proc,ProcNumber)
end  
