mutable struct CacheRKEStruct{FT<:AbstractFloat,
                           AT4<:AbstractArray,
                           AT5<:AbstractArray}
  Vn::AT4
  f::AT5
end

function Cache(backend,FT,IntMethod::RungeKuttaExMethod,FE,M,nz,NumV)
  NumG = FE.NumG
  NumI = FE.NumI
  nStage = IntMethod.nStage
  Vn = KernelAbstractions.zeros(backend,FT,M,nz,NumG,NumV)
  f = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV,nStage)
  return CacheRKEStruct{FT,
                     typeof(Vn),
                     typeof(f)}(
    Vn,
    f,
  )
end

function TimeIntegration!(RK::RungeKuttaExMethod,V,dt,Fcn,Aux,Jac,FE,Metric,Phys,Cache,JCache,
  Exchange,Global,Param,DiscType)

  f=Cache.f
  Vn=Cache.Vn
  @views VnI = Vn[:,:,1:size(f,3),1:size(f,4)]
  dtau, = dt
  FcnFull, = Fcn

  @inbounds for iStage=1:RK.nStage
    @. VnI = V
    @inbounds for jStage=1:iStage-1
      @views @. VnI = VnI + dt * RK.A[iStage,jStage] * f[:,:,:,:,jStage]
    end
    @views FcnFull(f[:,:,:,:,iStage],Vn,FE,Metric,Phys,Aux,Exchange,Global,DiscType)
  end
  @inbounds for iStage=1:RK.nStage
    @views @. V = V + dt * RK.b[iStage] * f[:,:,:,:,iStage]
  end
end

function RungeKuttaExplicit!(time,V,dt,Fcn,CG,Global,Param,DiscType)
  RK=Global.TimeStepper.RK
  f=Global.Cache.f
  Vn=Global.Cache.Vn

  @. Vn = V
  @inbounds for iStage=1:RK.nStage
    @. V = Vn
    @inbounds for jStage=1:iStage-1
      @views @. V = V + dt * RK.ARKE[iStage,jStage] * f[:,:,:,jStage]
    end
    @views Fcn(f[:,:,:,iStage],V,time + RK.cRKE[iStage] * dt,CG,Global,Param)
  end
  @. V = Vn
  @inbounds for iStage=1:RK.nStage
    @views @. V = V + dt * RK.bRKE[iStage] * f[:,:,:,iStage]
  end
end

function RungeKuttaExplicitLS!(time,V,dt,Fcn,CG,Global,Param)
  RK=Global.RK
  @views f=Global.Cache.f[:,:,:,1]
  Vn=Global.Cache.Vn

  @. Vn = V
  @inbounds for iStage=1:RK.nStage-1
    Global.Model.Param.time = time + RK.cRKE[iStage] * dt
    Fcn(f,V,CG,Global,Param)
    @. V = Vn + dt * RK.ARKE[iStage+1,iStage] * f
  end
  Fcn(f,V,CG,Global)
  @. V = Vn + dt * f
end

