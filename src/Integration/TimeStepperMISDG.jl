# MIS Method

function TimeStepperMIS(ROS,MIS,U,FcnSlow,FcnFast,Jac,dtauSmall,dtau,nIter,nPrint,DG,Exchange,Metric,Trans,
  Phys,Param,Grid,Global,VelForm)
	
  backend = get_backend(U)
  FTB = eltype(U)
  Proc = Global.ParallelCom.Proc
  ProcNumber = Global.ParallelCom.ProcNumber
  Model = Global.Model
  NumV = Model.NumV
  NumAux = Model.NumAux
  M = size(U,1)
  nz = size(U,2)
  CacheU = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumG,NumV)
  Un = CacheU
  CacheAux = CacheAuxStruct(backend,FTB,DG,M,nz,Model,Grid)
  Aux = CacheAux.Aux
  DGSEM.GeoPot(Aux,DG,Metric,Exchange,Global)
  @views UI = U[:,:,1:DG.NumI,:]
  @views UnI = Un[:,:,1:DG.NumI,1:NumV]
  nStage = ROS.nStage
  k = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV,nStage)

  M = DG.OrdPolyZ + 1
  dz = Metric.dz 
  Outputs.unstructured_vtkSphere(U,Trans,DG,Metric,Phys,Global,Proc,ProcNumber;Thermo=Aux)
  JCache = DGSEM.JacDGVert{FTB}(backend,M,nz,DG.NumI)
  SlowStages = MIS.nStage
  dZn = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV)
  Sdu = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV,MIS.nStage)
  Yn = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV,MIS.nStage+1)

  @time begin
    @inbounds for i = 1 : nIter
      time_elapsed = @elapsed begin
        @views @. Yn[:,:,:,:,1] = UI
        @inbounds for iStage = 1 : SlowStages
          @. UnI = Yn[:,:,:,:,iStage]
          @views FcnSlow(Sdu[:,:,:,:,iStage],Un,DG,Metric,Phys,CacheAux,Exchange,Global,VelForm)
          @views @.  Yn[:,:,:,:,iStage+1] = UI
          @views InitialConditionMISDG!(Yn[:,:,:,:,iStage+1], Yn, UI, MIS.alpha, iStage)
          @views SlowTendencyDG!(dZn, Yn, Sdu, UI, MIS.d, MIS.gamma, MIS.beta, dtau, iStage)
          @views RosenbrockMISDG!(ROS,dtauSmall,MIS.d[iStage+1]*dtau,Yn[:,:,:,:,iStage+1],dZn,FcnFast,Jac,DG,
            Exchange,Metric,Phys,Param,Grid,Global,JCache,CacheAux,Un,k,VelForm)
        end 
        @views @. UI = Yn[:,:,:,:,end] 
      end  
      percent = i/nIter*100
      if Proc == 1
        @info "Iteration: $i took $time_elapsed, $percent% complete"
      end
      if mod(i,nPrint) == 0 || i == nIter
        Outputs.unstructured_vtkSphere(UI,Trans,DG,Metric,Phys,Global,Proc,ProcNumber;Thermo=Aux)
      end
    end
  end
end  

function SlowTendencyDG!(dZn, Yn, Sdu, u, d, gamma, beta, dtL, stage)
  @views @. dZn = (beta[stage+1,1] / d[stage+1]) * Sdu[:,:,:,:,1]
  @inbounds for j in 2 : stage
    @views @. dZn += (gamma[stage+1,j] / (dtL * d[stage+1])) * (Yn[:,:,:,:,j] - u) + 
      (beta[stage+1,j] / d[stage+1]) * Sdu[:,:,:,:,j]
  end
end

function InitialConditionMISDG!(Zn0, Yn, u, alpha, stage)
  @inbounds for j in 2 : stage
    @views @. Zn0 += alpha[stage+1, j] * (Yn[:,:,:,:,j] - u)
  end
end

function RosenbrockMISDG!(ROS,dt_small,dt,U,FSlow,Fcn,Jac,DG,Exchange,Metric,Phys,Param,Grid,Global,
  JCache,CacheAux,Un,k,VelForm)
  numit = ceil(Int,dt/dt_small)
  FTB = eltype(U)
  Model = Global.Model
  NumV = Model.NumV
  M = size(U,1)
  nz = size(U,2)
  nStage = ROS.nStage
  M = DG.OrdPolyZ + 1
  dz = Metric.dz
  dtau = dt/numit
  @views UI = U[:,:,1:DG.NumI,:]
  @views UIn = Un[:,:,1:DG.NumI,1:size(U,4)]

  fac = dtau * ROS.gamma
  Jac(U,fac,DG,Metric,Phys,CacheAux,JCache,Global,VelForm)
  @inbounds for i = 1 : numit
    @inbounds for iStage = 1 : nStage
      @. UIn = UI
      @inbounds for jStage = 1 : iStage-1
        @views @. UIn = UIn + ROS.a[iStage,jStage] * k[:,:,:,:,jStage]
      end
      @views Fcn(k[:,:,:,:,iStage],Un,DG,Metric,Phys,CacheAux,Exchange,Global,VelForm)
      @views @. k[:,:,:,:,iStage] += FSlow 
      @inbounds for jStage = 1 : iStage - 1
        fac = ROS.c[iStage,jStage] / dtau
        @views @. k[:,:,:,:,iStage] += fac * k[:,:,:,:,jStage]
       end
       @views Solve!(k[:,:,:,:,iStage],k[:,:,:,:,iStage],JCache,dtau*ROS.gamma,DG,Metric,Global,VelForm)
     end
     @inbounds for iStage = 1 : nStage
       @views @. UI = UI + ROS.m[iStage] * k[:,:,:,:,iStage]
    end
  end 
end  
