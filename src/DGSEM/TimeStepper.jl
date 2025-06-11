function Rosenbrock(ROS,U,Fcn,dtau,IterTime,nPrint,DG,Exchange,Metric,Trans,Phys,Param,Grid,Global)

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
  CacheF = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV)
  @views UI = U[:,:,1:DG.NumI,:]
  @views UnI = Un[:,:,1:DG.NumI,:]
  nStage = ROS.nStage
  k = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV,nStage)

  dSdS,dSdM,dMdS,dMdM = InitJacDG(DG,nz,Param)
  N = size(dSdS,1)
  M = DG.OrdPolyZ + 1


  Outputs.unstructured_vtkSphere(U,Trans,DG,Metric,Phys,Global,Proc,ProcNumber)
  kLoc = zeros(3*N)
  @inbounds for i = 1 : IterTime
    if Proc == 1
      @show i,nPrint
    end  
    fac = (1.0 / (dtau * ROS.gamma))
    Jac = JacDG(U,DG,fac,dSdS,dSdM,dMdS,dMdM,Metric.dz,Phys)  
    @inbounds for iStage = 1 : nStage
      @. UnI = UI
      @inbounds for jStage = 1 : iStage-1
        @views @. UnI = UnI + ROS.a[iStage,jStage] * k[:,:,:,:,jStage]
      end
      @views Fcn(k[:,:,:,:,iStage],Un,DG,Model,Metric,Exchange,Grid,CacheU,CacheF,Phys,Global,Grid.Type)
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
#       @views ldiv!(Jac[ID],reshape(k[:,:,ID,[1,4,5],iStage],3*N))
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
    if mod(i,nPrint) == 0 || i == IterTime
      if Proc == 1
        @show "Print",i
      end
      Outputs.unstructured_vtkSphere(U,Trans,DG,Metric,Phys,Global,Proc,ProcNumber)
    end
  end  
end  

function RosenbrockEul(ROS,U,Fcn,dtau,IterTime,nPrint,DG,Exchange,Metric,Trans,Phys,Param,Grid,Global)

  backend = get_backend(U)
  FTB = eltype(U)
  Proc = Global.ParallelCom.Proc
  ProcNumber = Global.ParallelCom.ProcNumber
  Model = Global.Model
  NumV = Model.NumV
  NumAux = Model.NumAux
  nz = size(U,2)
  CacheU = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumG,NumV+NumAux)
  CacheF = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV)
  fU = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV)
  Un = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV)
  @views UI = U[:,:,1:DG.NumI,:]
  nStage = ROS.nStage
  k = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV,nStage)

  dSdS,dSdM,dMdS,dMdM = InitJacDG(DG,nz,Param)
  N = size(dSdS,1)
  M = DG.OrdPolyZ + 1


  Outputs.unstructured_vtkSphere(U,Trans,DG,Metric,Phys,Global,Proc,ProcNumber)
  @inbounds for i = 1 : IterTime
    if Proc == 1
      @show "RosEul",i,nPrint
    end  
    Jac = JacDG(U,DG,dSdS,dSdM,dMdS,dMdM,Metric.dz,Phys)  
    for ID = 1 : DG.NumI
       Jac[ID] = (1.0 / dtau) * sparse(I,3*N,3*N) - Jac[ID]
    end  
    Fcn(fU,U,DG,Model,Metric,Exchange,Grid,CacheU,CacheF,Phys,Global,Grid.Type)
    for ID = 1 : DG.NumI
      @views fU[:,:,ID,[1,4,5]] .= reshape(Jac[ID] \ reshape(fU[:,:,ID,[1,4,5]],3*N),M,nz,3)
    end  
    @views @. fU[:,:,:,2:3] = fU[:,:,:,2:3,:] * dtau 
    @. UI = UI + fU
    if mod(i,nPrint) == 0 || i == IterTime
      if Proc == 1
        @show "Print",i
      end
      Outputs.unstructured_vtkSphere(U,Trans,DG,Metric,Phys,Global,Proc,ProcNumber)
    end
  end  
end  

function RK3(U,Fcn,dtau,IterTime,nPrint,DG,Exchange,Metric,Trans,Phys,Grid,Global)

  backend = get_backend(U)
  FTB = eltype(U)
  Proc = Global.ParallelCom.Proc
  ProcNumber = Global.ParallelCom.ProcNumber
  Model = Global.Model
  NumV = Model.NumV
  NumAux = Model.NumAux
  nz = size(U,1)
  CacheU = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumG,NumV+NumAux)
  @views UNew = CacheU[:,:,:,1:NumV]
  CacheF = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV)
  FU = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV)
  @views UI = U[:,:,1:DG.NumI,:]
  @views UNewI = UNew[:,:,1:DG.NumI,:]
  
  Outputs.unstructured_vtkSphere(U,Trans,DG,Metric,Phys,Global,Proc,ProcNumber)
  @inbounds for i = 1 : IterTime
      if Proc == 1
        @show i,nPrint
      end  

    @. UNewI = UI 

    Fcn(FU,UNew,DG,Model,Metric,Exchange,Grid,CacheU,CacheF,Phys,Global,Grid.Type)
    fac = FTB(1/3 * dtau)
    @. UNewI = UI + fac * FU

    Fcn(FU,UNew,DG,Model,Metric,Exchange,Grid,CacheU,CacheF,Phys,Global,Grid.Type)
    fac = FTB(1/2 * dtau)
    @. UNewI = UI + fac * FU

    Fcn(FU,UNew,DG,Model,Metric,Exchange,Grid,CacheU,CacheF,Phys,Global,Grid.Type)
    fac = FTB(dtau)
    @. UI = UI + fac * FU

    if mod(i,nPrint) == 0 || i == IterTime
      if Proc == 1
        @show "Print",i
      end  
      Outputs.unstructured_vtkSphere(U,Trans,DG,Metric,Phys,Global,Proc,ProcNumber)
    end
  end
end


function RK1(U,Fcn,dtau,IterTime,nPrint,DG,Exchange,Metric,Trans,Phys,Grid,Global)

  backend = get_backend(U)
  FTB = eltype(U)
  Proc = Global.ParallelCom.Proc
  ProcNumber = Global.ParallelCom.ProcNumber
  Model = Global.Model
  NumV = Model.NumV
  NumAux = Model.NumAux
  nz = size(U,1)
  CacheU = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumG,NumV+NumAux)
  CacheF = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV)
  FU = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV)
  @views UI = U[:,:,1:DG.NumI,:]
  
  Outputs.unstructured_vtkSphere(U,Trans,DG,Metric,Phys,Global,Proc,ProcNumber)
  @inbounds for i = 1 : IterTime
      if Proc == 1
        @show i,nPrint
      end  

    Fcn(FU,U,DG,Model,Metric,Exchange,Grid,CacheU,CacheF,Phys,Global,Grid.Type)
    fac = FTB(dtau)
    @. UI = UI + fac * FU

    if mod(i,nPrint) == 0 || i == IterTime
      if Proc == 1
        @show "Print",i
      end  
      Outputs.unstructured_vtkSphere(U,Trans,DG,Metric,Phys,Global,Proc,ProcNumber)
    end
  end
end
