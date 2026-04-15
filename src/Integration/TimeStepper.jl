function TimeStepper(IntMethod,dt,U,Fcn,Jac,FE,Exchange,Metric,Trans,Phys,Param,Grid,
  Global,ElemType::Grids.ElementType,VelForm)

  backend = get_backend(U)
  FT = eltype(U)
  TimeStepper = Global.TimeStepper
  Output = Global.Output
  Proc = Global.ParallelCom.Proc
  ProcNumber = Global.ParallelCom.ProcNumber
  NumV = Global.Model.NumV
  NumTr = Global.Model.NumTr
  NumAux = Global.Model.NumAux
  ND = Global.Model.NDEDMF
  nz = Global.Grid.nz
  M = FE.Mz

  CacheInt = Cache(backend,FT,IntMethod,FE,M,nz,NumV)
  CacheAux = CacheAuxStruct(backend,FT,FE,M,nz,Global.Model,Grid)
  if IntMethod.JacComp
    JCache = CacheJac(backend,FT,M,nz,Global.Model,FE)
  else
    JCache = nothing  
  end  
  Aux = CacheAux.Aux
  if Global.Model.GPAuxPos > 0
    DGSEM.GeoPot(Aux,FE,Metric,Exchange,Global)
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


  Outputs.unstructured_vtkSphere(U,Trans,FE,Metric,Phys,Global,Proc,ProcNumber;Thermo=Aux)
# TimeIntegration!(IntMethod,U,dt,Fcn,CacheAux,Jac,FE,Metric,Phys,CacheInt,JCache,Exchange,
#   Global,Param,VelForm)
  @time begin
    @inbounds for i = 1 : nIter
      Δt = @elapsed begin
        TimeIntegration!(IntMethod,U,dt,Fcn,CacheAux,Jac,FE,Metric,Phys,CacheInt,JCache,Exchange,
          Global,Param,VelForm)
        if mod(i,PrintInt) == 0
          @. @views U[:,:,FE.BoundaryDoF,3] = FT(0)  
          Outputs.unstructured_vtkSphere(U,Trans,FE,Metric,Phys,Global,Proc,ProcNumber;Thermo=Aux)
        end
      end
      percent = i/nIter*100
      if Proc == 1
        @info "Iteration: $i took $Δt, $percent% complete"
      end
    end
  end  
  @. @views U[:,:,FE.BoundaryDoF,3] = FT(0)  
  Outputs.unstructured_vtkSphere(U,Trans,FE,Metric,Phys,Global,Proc,ProcNumber;Thermo=Aux)

end

