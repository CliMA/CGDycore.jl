function TimeStepper!(U,Trans,CG,Global,Param)  
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
  elseif IntMethod == "IMEX"   
    TimeStepper.IMEX=IMEXMethod(Table)
  elseif IntMethod == "MIS"  
    TimeStepper.MIS=MISMethod(Table)
  elseif IntMethod == "LinIMEX"  
    TimeStepper.LinIMEX=LinIMEXMethod(Table)
  end

# Simulation period
  time=[0.0]
  dtau = TimeStepper.dtau
  SimDays = TimeStepper.SimDays
  SimHours = TimeStepper.SimHours
  SimMinutes = TimeStepper.SimMinutes
  SimSeconds = TimeStepper.SimSeconds
  PrintDays = Output.PrintDays
  PrintHours = Output.PrintHours
  PrintSeconds = Output.PrintSeconds
  PrintStartDays = Output.PrintStartDays
  nIter=ceil((24*3600*SimDays+3600*SimHours+60*SimMinutes+SimSeconds)/dtau)
  PrintInt=ceil((24*3600*PrintDays+3600*PrintHours+PrintSeconds)/dtau)
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
  elseif IntMethod == "RosenbrockSSP"
    Global.J = JStruct(NumG,nz,NumTr)
    Global.Cache.k=zeros(size(U[:,:,1:NumV])..., TimeStepper.ROS.nStage);
    Global.Cache.fV=zeros(size(U))
    Global.Cache.VS=zeros(size(U[:,:,NumV+1:end])..., TimeStepper.ROS.nStage+1);
    Global.Cache.fS=zeros(size(U[:,:,NumV+1:end])..., TimeStepper.ROS.nStage);
    Global.Cache.fRhoS=zeros(size(U[:,:,1])..., TimeStepper.ROS.nStage);
    Global.Cache.RhoS=zeros(size(U[:,:,1])..., TimeStepper.ROS.nStage+1);
  elseif IntMethod == "LinIMEX"
    Global.J = JStruct(NumG,nz,NumTr)
    Global.Cache.Ymyn=zeros(size(U[:,:,1:NumV+NumTr])..., TimeStepper.LinIMEX.nStage-1);
    Global.Cache.fV=zeros(size(U))
    Global.Cache.f=zeros(size(U)..., TimeStepper.LinIMEX.nStage)
    Global.Cache.fV=zeros(size(U))
    Global.Cache.Vn=zeros(size(U))
  elseif IntMethod == "RungeKutta"
    Global.Cache.f=zeros(size(U)..., TimeStepper.RK.nStage)
  elseif IntMethod == "MIS"
    Global.Cache.f=zeros(size(U)..., TimeStepper.MIS.nStage)
    Global.Cache.VS=zeros(size(U)..., TimeStepper.MIS.nStage - 1)
    Global.J = JStruct(NumG,nz,NumTr)
    Global.Cache.k=zeros(size(U[:,:,1:NumV+NumTr])..., TimeStepper.ROS.nStage);
    Global.Cache.fV=zeros(size(U))
    Global.Cache.Vn=zeros(size(U))
    Global.Cache.R=zeros(size(U))
    Global.Cache.dZ=zeros(size(U))  
  elseif IntMethod == "IMEX"
    Global.J = JStruct(NumG,nz,NumTr)
    Global.Cache.fV=zeros(size(U))
    Global.Cache.R=zeros(size(U))
    Global.Cache.dZ=zeros(size(U))
    Global.Cache.Y=zeros(size(U[:,:,1:NumV+NumTr])..., TimeStepper.IMEX.nStage);
    Global.Cache.Z=zeros(size(U[:,:,1:NumV+NumTr])..., TimeStepper.IMEX.nStage);
    Global.Cache.Vn=zeros(size(U))  
  end

 # Boundary values
  if Global.Model.SurfaceFlux
    Global.Cache.TSurf=ProjectSurf(fTSurf,0.0,CG,Global,Param)
  end


# Print initial conditions
  @show "Print initial conditions"
  unstructured_vtkSphere(U,TransSphereX,CG,Global,Proc,ProcNumber)

  @show "Choose integration method"
  if IntMethod == "Rosenbrock"
    @time begin
      for i=1:nIter
        Δt = @elapsed begin
          RosenbrockSchur!(U,dtau,Fcn!,JacSchur!,CG,Global,Param);
#         RosenbrockSchur!(U,dtau,FcnNHCurlVecI!,JacSchur!,CG,Global,Param);
          time[1] += dtau
          if mod(i,PrintInt) == 0 && i >= PrintStartInt
            unstructured_vtkSphere(U,Trans,CG,Global,Proc,ProcNumber)
          end
        end
        percent = i/nIter*100
        @show sum(abs.(U))
        @info "Iteration: $i took $Δt, $percent% complete"
      end
    end
  elseif IntMethod == "RosenbrockD"
    @time begin
      for i=1:nIter
        Δt = @elapsed begin
          RosenbrockDSchur!(U,dtau,FcnNHCurlVec!,JacSchur!,CG,Global);
          time[1] += dtau
          if mod(i,PrintInt) == 0 && i >= PrintStartInt
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
          if mod(i,PrintInt) == 0 && i >= PrintStartInt
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
          if mod(i,PrintInt) == 0 && i >= PrintStartInt
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
          if mod(i,PrintInt) == 0 && i >= PrintStartInt
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
          if mod(i,PrintInt) == 0 && i >= PrintStartInt
            unstructured_vtkSphere(U,TransSphereX,CG,Global,Proc,ProcNumber)
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
          RungeKuttaExplicit!(U,dtau,Fcn!,CG,Global,Param)
          time[1] += dtau
          if mod(i,PrintInt)==0 && i >= PrintStartInt
            unstructured_vtkSphere(U,Trans,CG,Global,Proc,ProcNumber)
          end
        end
        percent = i/nIter*100
        @show sum(abs.(U))
        @info "Iteration: $i took $Δt, $percent% complete"
      end
    end
  else
    error("Bad IntMethod")
  end
end  

function TimeStepperAdvection!(U,Trans,CG,Global,Param)  
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
  PrintDays = Output.PrintDays
  PrintHours = Output.PrintHours
  PrintSeconds = Output.PrintSeconds
  PrintStartDays = Output.PrintStartDays
  nIter=ceil((24*3600*SimDays+3600*SimHours+60*SimMinutes+SimSeconds)/dtau)
  PrintInt=ceil((24*3600*PrintDays+3600*PrintHours+PrintSeconds)/dtau)
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
    Global.Cache.fS=zeros(nz,NumG,NumTr,TimeStepper.SSP.nStage)
    Global.Cache.VS=zeros(nz,NumG,NumTr,TimeStepper.SSP.nStage+1)
    Global.Cache.fV=zeros(size(U))
    Global.Cache.RhoS=zeros(size(U[:,:,1])..., TimeStepper.SSP.nStage+1)
    Global.Cache.fRhoS=zeros(size(U[:,:,1])..., TimeStepper.SSP.nStage)
  end


# Print initial conditions
  @show "Print initial conditions"
  unstructured_vtkSphere(U,TransSphereX,CG,Global,Proc,ProcNumber)

  @show "Choose integration method"
  if IntMethod == "Rosenbrock"
    @time begin
      for i=1:nIter
        Δt = @elapsed begin
          RosenbrockSchur!(U,dtau,FcnTracer!,JacSchur!,CG,Global,Param);
          time[1] += dtau
          if mod(i,PrintInt) == 0 && i >= PrintStartInt
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
    @show "Hallo Oswald"
    @show nIter,"Start"
    @time begin
      for i=1:nIter
        Δt = @elapsed begin
          SSPRungeKutta!(time[1],U,dtau,FcnTracer!,CG,Global,Param)
          @show sum(abs.(U))
          time[1] += dtau
          if mod(i,PrintInt)==0 && i >= PrintStartInt
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
function TimeStepperAdvectionConv!(U,Trans,CG,Global,Param)  
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
  PrintDays = Output.PrintDays
  PrintHours = Output.PrintHours
  PrintSeconds = Output.PrintSeconds
  PrintStartDays = Output.PrintStartDays
  nIter=ceil((24*3600*SimDays+3600*SimHours+60*SimMinutes+SimSeconds)/dtau)
  PrintInt=ceil((24*3600*PrintDays+3600*PrintHours+PrintSeconds)/dtau)
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
    Global.Cache.fS=zeros(nz,NumG,NumTr,TimeStepper.SSP.nStage)
    Global.Cache.VS=zeros(nz,NumG,NumTr,TimeStepper.SSP.nStage+1)
    Global.Cache.fV=zeros(size(U))
    Global.Cache.RhoS=zeros(size(U[:,:,1])..., TimeStepper.SSP.nStage+1)
    Global.Cache.fRhoS=zeros(size(U[:,:,1])..., TimeStepper.SSP.nStage)
  end


# Print initial conditions
  @show "Print initial conditions"
  unstructured_vtkSphere(U,TransSphereX,CG,Global,Proc,ProcNumber)

  @show "Choose integration method"
  if IntMethod == "Rosenbrock"
    @time begin
      for i=1:nIter
        Δt = @elapsed begin
          RosenbrockSchur!(U,dtau,FcnTracer!,JacSchur!,CG,Global,Param);
          time[1] += dtau
          if mod(i,PrintInt) == 0 && i >= PrintStartInt
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
    @show "Hallo Oswald"
    @show nIter,"Start"
    @time begin
      for i=1:nIter
        Δt = @elapsed begin
          @show "B",sum(abs.(U))
          SSPRungeKutta!(time[1],U,dtau,FcnTracerConv!,CG,Global,Param)
          @show sum(abs.(U[:,:,1]))
          @show sum(abs.(U[:,:,2]))
          time[1] += dtau
          if mod(i,PrintInt)==0 && i >= PrintStartInt
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
