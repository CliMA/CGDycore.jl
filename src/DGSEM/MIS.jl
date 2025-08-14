# MIS Method

function MIS_Method(ROS,MIS,U,FcnSlow,FcnFast,dtauSmall,dtau,nIter,nPrint,DG,Exchange,Metric,Trans,Phys,Param,Grid,Global)
	
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
  Outputs.unstructured_vtkSphere(U,Trans,DG,Metric,Phys,Global,Proc,ProcNumber)
  Jac = JacDGVert{FTB}(backend,M,nz,DG.NumI)
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
          @views FcnSlow(Sdu[:,:,:,:,iStage],Un,DG,Model,Metric,Exchange,Grid,CacheU,CacheS,Phys,Global,Grid.Type, Param)
          @views @.  Yn[:,:,:,:,iStage+1] = UI
          @views InitialConditionMIS!(Yn[:,:,:,:,iStage+1], Yn, UI, MIS.alfa, iStage)
          @views SlowTendency!(dZn, Yn, Sdu, UI, MIS.d, MIS.gamma, MIS.beta, dtau, iStage)
          @views RosenbrockMIS!(ROS, dtauSmall, MIS.d[iStage+1] * dtau,Yn[:,:,:,:,iStage+1],dZn,FcnFast,DG,
            Exchange,Metric,Phys,Param,Grid,Global,Jac,CacheU,CacheS,k)
        end 
        @views @. UI = Yn[:,:,:,:,end] 
      end  
      percent = i/nIter*100
      if Proc == 1
        @info "Iteration: $i took $time_elapsed, $percent% complete"
      end
      if mod(i,nPrint) == 0 || i == nIter
        Outputs.unstructured_vtkSphere(UI,Trans,DG,Metric,Phys,Global,Proc,ProcNumber)
      end
    end
  end
end  

function SlowTendency!(dZn, Yn, Sdu, u, d, gamma, beta, dtL, stage)
  @views @. dZn = (beta[stage+1,1] / d[stage+1]) * Sdu[:,:,:,:,1]
  @inbounds for j in 2 : stage
    @views @. dZn += (gamma[stage+1,j] / (dtL * d[stage+1])) * (Yn[:,:,:,:,j] - u) + 
      (beta[stage+1,j] / d[stage+1]) * Sdu[:,:,:,:,j]
  end
end

function InitialConditionMIS!(Zn0, Yn, u, alfa, stage)
  @inbounds for j in 2 : stage
    @views @. Zn0 += alfa[stage+1, j] * (Yn[:,:,:,:,j] - u)
  end
end

function RosenbrockMIS!(ROS,dt_small,dt,U,FSlow,Fcn,DG,Exchange,Metric,Phys,Param,Grid,Global,
  Jac,CacheU,CacheS,k)
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
  @views Un = CacheU[:,:,:,1:NumV] 
  @views UIn = Un[:,:,1:DG.NumI,:]

  fac = (1.0 / (dtau * ROS.gamma))
  FillJacDGVert!(Jac,U,DG,dz,fac,Phys,Param)
  SchurBoundary!(Jac)
  @inbounds for i = 1 : numit
    @inbounds for iStage = 1 : nStage
      @. UIn = UI
      @inbounds for jStage = 1 : iStage-1
        @views @. UIn = UIn + ROS.a[iStage,jStage] * k[:,:,:,:,jStage]
      end
      @views Fcn(k[:,:,:,:,iStage],Un,DG,Model,Metric,Exchange,Grid,CacheU,CacheS,Phys,Global,Grid.Type,Param)
      @views @. k[:,:,:,:,iStage] += FSlow 
      @inbounds for jStage = 1 : iStage - 1
        fac = ROS.c[iStage,jStage] / dtau
        @views @. k[:,:,:,:,iStage] += fac * k[:,:,:,:,jStage]
       end
       @views Solve!(Jac,k[:,:,:,:,iStage])
       @views @. k[:,:,:,2:3,iStage] *= (dtau * ROS.gamma)
     end
     @inbounds for iStage = 1 : nStage
       @views @. UI = UI + ROS.m[iStage] * k[:,:,:,:,iStage]
    end
  end 
end  
