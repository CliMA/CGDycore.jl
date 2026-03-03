mutable struct CacheMISStruct{FT<:AbstractFloat,
                           AT4<:AbstractArray,
                           AT5<:AbstractArray}
  Vn::AT4
  dZn::AT4
  fV::AT4
  Sdu::AT5
  Yn::AT5
# ROS
  k::AT5
end

function Cache(backend,FT,IntMethod::MISMethod,FE,M,nz,NumV) 
  NumG = FE.NumG
  NumI = FE.NumI
  nStage = IntMethod.nStage
  Vn = KernelAbstractions.zeros(backend,FT,M,nz,NumG,NumV)

  dZn = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV)
  fV = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV)
  Sdu = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV,nStage)
  Yn = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV,nStage+1)

  k = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV,IntMethod.FastMethod.nStage)
  return CacheMISStruct{FT,
                     typeof(Vn),
                     typeof(k)}(
    Vn,
    dZn,
    fV,
    Sdu,
    Yn,
    k,
  )
end

function TimeIntegration!(MIS::MISMethod,V,dt,Fcn,Aux,Jac,FE,Metric,Phys,Cache,JCache,Exchange,
  Global,Param,Type)

  FcnSlow,FcnFast = Fcn
  dtSlow,dtFast = dt
  Fast = MIS.FastMethod

  nStage = MIS.nStage
  k = Cache.k
  Vn = Cache.Vn
  Sdu = Cache.Sdu
  dZn = Cache.dZn
  Yn = Cache.Yn
  @views VnI = Vn[:,:,1:size(V,3),1:size(V,4)]
  Global.TimeStepper.dtauStage = dtSlow

  @views @. Yn[:,:,:,:,1] = V
  @inbounds for iStage = 1 : nStage
    @. VnI = Yn[:,:,:,:,iStage]
    @views FcnSlow(Sdu[:,:,:,:,iStage],Vn,FE,Metric,Phys,Aux,Exchange,Global,Type)
    @views @.  Yn[:,:,:,:,iStage+1] = V
    @views InitialConditionMIS!(Yn[:,:,:,:,iStage+1], Yn, V, MIS.alpha, iStage)
    @views SlowTendency!(dZn, Yn, Sdu, V, MIS.d, MIS.gamma, MIS.beta, dtSlow, iStage)
    @views TimeStepperFast!(Fast,dtFast,MIS.d[iStage+1]*dtSlow,Yn[:,:,:,:,iStage+1],dZn,FcnFast,Jac,FE,
       Exchange,Metric,Phys,Param,Global,Cache,JCache,Aux,Vn,k,Type)
  end
  @views @. V = Yn[:,:,:,:,end]
end  

function SlowTendency!(dZn, Yn, Sdu, u, d, gamma, beta, dtL, stage)
  @views @. dZn = (beta[stage+1,1] / d[stage+1]) * Sdu[:,:,:,:,1]
  @inbounds for j in 2 : stage
    @views @. dZn += (gamma[stage+1,j] / (dtL * d[stage+1])) * (Yn[:,:,:,:,j] - u) +
      (beta[stage+1,j] / d[stage+1]) * Sdu[:,:,:,:,j]
  end
end

function InitialConditionMIS!(Zn0, Yn, u, alpha, stage)
  @inbounds for j in 2 : stage
    @views @. Zn0 += alpha[stage+1, j] * (Yn[:,:,:,:,j] - u)
  end
end

function TimeStepperFast!(IntMethod,dt,tEnd,U,FSlow,Fcn,Jac,FE,Exchange,Metric,Phys,Param,Global,
  Cache,JCache,CacheAux,Un,k,VelForm)
  numit = ceil(Int,tEnd/dt)
  FTB = eltype(U)
  dtau = tEnd/numit
  fac = dtau * IntMethod.gammaD
  Jac(U,fac,FE,Metric,Phys,CacheAux,JCache,Global,VelForm)
  @inbounds for i = 1 : numit
    TimeIntegrationFast!(IntMethod,U,dtau,Fcn,FSlow,CacheAux,FE,Metric,Phys,Cache,JCache,Exchange,
      Global,Param,VelForm)
  end
end

function TimeStepperFast!(IntMethod,dt,tEnd,U,D,FSlow,Fcn,Jac,FE,Exchange,Metric,Phys,Param,Global,
  Cache,JCache,CacheAux,Un,k,VelForm)
  numit = ceil(Int,tEnd/dt)
  FTB = eltype(U)
  dtau = tEnd/numit
  fac = dtau * IntMethod.gammaD
  Jac(U,fac,FE,Metric,Phys,CacheAux,JCache,Global,VelForm)
  @. U -= D
  @inbounds for i = 1 : numit
    TimeIntegrationFast!(IntMethod,U,dtau,Fcn,FSlow,CacheAux,FE,Metric,Phys,Cache,JCache,Exchange,
      Global,Param,VelForm)
  end
  @. U += D
end
