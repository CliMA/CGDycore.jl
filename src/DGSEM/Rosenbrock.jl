function Rosenbrock!(U,dt,Fcn,FcnPrepare,Jac,DG,Metric,Phys,Cache,JCache,Exchange,Global,Param,Equation)

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
  @views Un = CacheU[:,:,:,1:NumV]
  CacheS = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV)
  @views UI = U[:,:,1:DG.NumI,:]
  @views UnI = Un[:,:,1:DG.NumI,:]
  nStage = ROS.nStage
  k = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV,nStage)

  M = DG.OrdPolyZ + 1
  dz = Metric.dz 

  Jac(Jac,dt*ROS.gamma,U,DG,Metric,Phys,Cache,Global,Param,Equation)
  @inbounds for iStage = 1 : nStage
    @. UnI = UI
    @inbounds for jStage = 1 : iStage-1
      @views @. UnI = UnI + ROS.a[iStage,jStage] * k[:,:,:,:,jStage]
    end
    @views Fcn(k[:,:,:,:,iStage],Un,DG,Model,Metric,Exchange,Grid,CacheU,CacheS,Phys,Global,Grid.Type)
    @inbounds for jStage = 1 : iStage - 1
      fac = ROS.c[iStage,jStage] / dt
      @views @. k[:,:,:,:,iStage] += fac * k[:,:,:,:,jStage]
    end
    @views Solve!(Jac,k[:,:,:,:,iStage])
    @views @. k[:,:,:,2:3,iStage] *= (dt * ROS.gamma)
  end
  @inbounds for iStage = 1 : nStage
    @views @. UI = UI + ROS.m[iStage] * k[:,:,:,:,iStage]
  end
end  
