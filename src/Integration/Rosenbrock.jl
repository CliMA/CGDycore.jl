mutable struct CacheROSStruct{FT<:AbstractFloat,
                           AT4<:AbstractArray,
                           AT5<:AbstractArray}
  Vn::AT4
  k::AT5
  fV::AT4
end

function Cache(backend,FT,IntMethod::RosenbrockMethod,FE,M,nz,NumV) 
  NumG = FE.NumG
  NumI = FE.NumI
  nStage = IntMethod.nStage
  Vn = KernelAbstractions.zeros(backend,FT,M,nz,NumG,NumV)
  fV = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV)
  k = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV,nStage)
  return CacheROSStruct{FT,
                     typeof(Vn),
                     typeof(k)}(
    Vn,
    k,
    fV,
  )
end


NVTX.@annotate function TimeIntegration!(ROS::RosenbrockMethod,V,dt,Fcn,Aux,Jac,FE,Metric,Phys,Cache,JCache,Exchange,
  Global,Param,DiscType)

  nStage = ROS.nStage
  k = Cache.k
  fV = Cache.fV
  Vn = Cache.Vn
  @views VnI = Vn[:,:,1:size(fV,3),1:size(fV,4)]
  dtau, = dt
  FcnFull, = Fcn
  Global.TimeStepper.dtauStage = dtau  

  JCache.CompTri = true
  @inbounds for iStage = 1 : nStage
    @. VnI = V
    @inbounds for jStage = 1 : iStage-1
      @views @. VnI = VnI + ROS.a[iStage,jStage] * k[:,:,:,:,jStage]
    end
    FcnFull(fV,Vn,FE,Metric,Phys,Aux,Exchange,Global,DiscType)
    if iStage == 1
      Jac(V,dtau*ROS.gammaD,FE,Metric,Phys,Aux,JCache,Global,DiscType)
    end  
    @inbounds for jStage = 1 : iStage - 1
      fac = ROS.c[iStage,jStage] / dtau
      @views @. fV = fV + fac * k[:,:,:,:,jStage]
    end
    @views Solve!(k[:,:,:,:,iStage],fV,JCache,dtau*ROS.gammaD,FE,Metric,Global,DiscType)
  end
  @inbounds for iStage = 1 : nStage
    @views @. V = V + ROS.m[iStage] * k[:,:,:,:,iStage]
  end
end

NVTX.@annotate function TimeIntegrationFast!(ROS::RosenbrockMethod,V,dt,Fcn,FSlow,Aux,FE,Metric,Phys,Cache,JCache,Exchange,
  Global,Param,DiscType)

  nStage = ROS.nStage
  k = Cache.k
  fV = Cache.fV
  Vn = Cache.Vn
  @views VnI = Vn[:,:,1:size(fV,3),1:size(fV,4)]
  Global.TimeStepper.dtauStage = dt

  JCache.CompTri = true
  @inbounds for iStage = 1 : nStage
    @. VnI = V
    @inbounds for jStage = 1 : iStage-1
      @views @. VnI = VnI + ROS.a[iStage,jStage] * k[:,:,:,:,jStage]
    end
    Fcn(fV,Vn,FE,Metric,Phys,Aux,Exchange,Global,DiscType)
    @. fV += FSlow
    @inbounds for jStage = 1 : iStage - 1
      fac = ROS.c[iStage,jStage] / dt
      @views @. fV = fV + fac * k[:,:,:,:,jStage]
    end
    @views Solve!(k[:,:,:,:,iStage],fV,JCache,dt*ROS.gammaD,FE,Metric,Global,DiscType)
  end
  @inbounds for iStage = 1 : nStage
    @views @. V = V + ROS.m[iStage] * k[:,:,:,:,iStage]
  end
end

