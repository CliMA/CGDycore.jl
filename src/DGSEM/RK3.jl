function TestJacobian(ROS,U,Fcn,dtau,IterTime,nPrint,DG,Exchange,Metric,Trans,Phys,Param,Grid,Global)

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
  fT = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV)
  Un = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV)
  @views UI = U[:,:,1:DG.NumI,:]
  nStage = ROS.nStage
  k = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV,nStage)

  dSdS,dSdM,dMdS,dMdM = InitJacDG(DG,nz,Param)
  N = size(dSdS,1)
  M = DG.OrdPolyZ + 1
  Fcn(fU,U,DG,Model,Metric,Exchange,Grid,CacheU,CacheF,Phys,Global,Grid.Type)

  JacDiff = zeros(3*N,3*N)
  ColJacDiff = zeros(M,nz,3)
  iC = 0
  for iv in [1,4,5]
    for iz = 1 : nz
      for k = 1 : M
        iC += 1
        temp = U[k,iz,1,iv]
        U[k,iz,1,iv] = (1 + 1.e-8) * U[k,iz,1,iv] + 1.e-8
        Fcn(fT,U,DG,Model,Metric,Exchange,Grid,CacheU,CacheF,Phys,Global,Grid.Type)
        ColJacDiff .= (fT[:,:,1,[1,4,5]] - fU[:,:,1,[1,4,5]]) ./ (U[k,iz,1,iv] - temp)
        JacDiff[:,iC] .= reshape(ColJacDiff,3*N)
        U[k,iz,1,iv] = temp
      end
    end
  end
  JacD = sparse(JacDiff)

  Jac = JacDG1(U,DG,dSdS,dSdM,dMdS,dMdM,Metric.dz,Phys)  
  return JacD, Jac[1]
end  

function Rosenbrock(ROS,U,Fcn,dtau,IterTime,nPrint,DG,Exchange,Metric,Trans,Phys,Param,Grid,Global)

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
      @show i,nPrint
    end  
    Jac = JacDG(U,DG,dSdS,dSdM,dMdS,dMdM,Metric.dz,Phys)  
    for ID = 1 : DG.NumI
      for i = 1 : Jac[ID].n  
        Jac[ID] .= (1.0 / (dtau * ROS.gamma)) * sparse(I,3*N,3*N) - Jac[ID]
      end
    end  
    @. Un = UI
    @inbounds for iStage = 1 : nStage
      @. UI = Un
      @inbounds for jStage = 1 : iStage-1
        @views @. UI = UI + ROS.a[iStage,jStage] * k[:,:,:,:,jStage]
      end
      Fcn(fU,U,DG,Model,Metric,Exchange,Grid,CacheU,CacheF,Phys,Global,Grid.Type)
      @inbounds for jStage = 1 : iStage - 1
        fac = ROS.c[iStage,jStage] / dtau
        @views @. fU = fU + fac * k[:,:,:,:,jStage]
      end
      for ID = 1 : DG.NumI
        @views k[:,:,ID,[1,4,5],iStage] .= reshape(Jac[ID] \ reshape(fU[:,:,ID,[1,4,5]],3*N),M,nz,3)
      end  
      @views @. k[:,:,:,2:3,iStage] = fU[:,:,:,2:3,:] * (dtau * ROS.gamma)
    end
    @. UI = Un
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
  CacheF = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV)
  FU = KernelAbstractions.zeros(backend,FTB,size(U,1),size(U,2),DG.NumI,NumV)
  UNew = similar(U)
  @views UI = U[:,:,1:DG.NumI,:]
  @views UNewI = UNew[:,:,1:DG.NumI,:]
  
  Outputs.unstructured_vtkSphere(U,Trans,DG,Metric,Phys,Global,Proc,ProcNumber)
  @inbounds for i = 1 : IterTime
      if Proc == 1
        @show i,nPrint
      end  

    Fcn(FU,U,DG,Model,Metric,Exchange,Grid,CacheU,CacheF,Phys,Global,Grid.Type)
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
