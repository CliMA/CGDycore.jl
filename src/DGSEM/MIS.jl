# MIS Method

function MIS_Method(ROS,MIS,U,FcnSlow,FcnFast,dtau,nIter,nPrint,DG,Exchange,Metric,Trans,Phys,Param,Grid,Global)
	
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
  dt_small = dtau * 20/340.0	
  Outputs.unstructured_vtkSphere(U,Trans,DG,Metric,Phys,Global,Proc,ProcNumber)
  Jac = JacDGVert{FTB}(backend,M,nz,DG.NumI)
  SlowStages = MIS.nStage
  Zn0 = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV)
  yn = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV)
  y = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV)
  dZn = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV)
  Sdu = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV,MIS.nStage)
  Yn = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV,MIS.nStage)
  @time begin
	@inbounds for i = 1 : nIter
		time_elapsed = @elapsed begin

		## MIS Slow
		@inbounds for iStage = 1:SlowStages
		@. Zn0 = UI
	InitialConditionMIS!(Zn0, Yn, UI, MIS.alfa, iStage)
		
		@. dZn = 0.0
	PreparationODEMIS!(dZn, Yn, Sdu, UI, MIS.d, MIS.gamma, MIS.beta, dtau, SlowStages)

		if iStage == 1
		@. Yn[:,:,:,:,iStage] = UI 
		else
	RosenbrockMIS!(ROS, dt_small, MIS.d[iStage] * dtau, iStage,Yn, Zn0, dZn, yn, y, k, FcnFast,DG,Exchange,Metric,Trans,Phys,Param,Grid,Global,Jac,CacheU,CacheS)
		end
		@views @. yn = Yn[:,:,:,:,iStage]
	@views FcnSlow(Sdu[:,:,:,:,iStage],yn,DG,Model,Metric,Exchange,Grid,CacheU,CacheS,Phys,Global,Grid.Type, Param)

		end ## Fine Stages
		
		@. UI = Yn[:,:,:,:,end] 
        if mod(i,nPrint) == 0 || i == nIter
          Outputs.unstructured_vtkSphere(UI,Trans,DG,Metric,Phys,Global,Proc,ProcNumber)
        end
		end
      		percent = i/nIter*100
      		if Proc == 1
        	@info "Iteration: $i took $time_elapsed, $percent% complete"
      		end
		## fine n_Iter
	end
    end

end  

function PreparationODEMIS!(dZn, Yn, Sdu, u, d, gamma, beta, dtL, stage)
	@inbounds for j in 1:(stage-1)
		@views @. dZn = dZn + 1/d[stage] * (1/dtL * gamma[stage,j] * (Yn[:,:,:,:,j] - u) + beta[stage,j] * Sdu[:,:,:,:,j])
	end
end

function InitialConditionMIS!(Zn0, Yn, u, alfa, stage)
	@inbounds for j in 1:(stage-1)
	@views @. Zn0 = Zn0 + alfa[stage, j] * (Yn[:,:,:,:,j] - u)
	end
end


function RosenbrockMIS!(ROS,dt_small, dt, stage,Yn,Zn0, dZn, UI, UIn, k, Fcn,DG,Exchange,Metric,Trans,Phys,Param,Grid,Global,Jac,CacheU,CacheS)
  numit = round(dt/dt_small)
	if numit <= 1.0
	   numit = 1.0
	end
  FTB = eltype(UI)
  Model = Global.Model
  NumV = Model.NumV
  M = size(UI,1)
  nz = size(UI,2)
  nStage = ROS.nStage
  M = DG.OrdPolyZ + 1
  dz = Metric.dz
  dtau = dt/numit
  @. UI = Zn0

        fac = (1.0 / (dtau * ROS.gamma))
	FillJacDGVert!(Jac,UI,DG,dz,fac,Phys,Param)
       SchurBoundary!(Jac)
    @inbounds for i = 1 : numit
        @inbounds for iStage = 1 : nStage
          @. UIn = UI
          @inbounds for jStage = 1 : iStage-1
            @views @. UIn = UIn + ROS.a[iStage,jStage] * k[:,:,:,:,jStage]
          end
          @views Fcn(k[:,:,:,:,iStage],UIn,DG,Model,Metric,Exchange,Grid,CacheU,CacheS,Phys,Global,Grid.Type, Param)
	@. k[:,:,:,:,iStage] += dZn 
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
	@. Yn[:,:,:,:,stage] = UI
end  
