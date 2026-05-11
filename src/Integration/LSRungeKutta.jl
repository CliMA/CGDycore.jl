mutable struct CacheLSRKStruct{FT<:AbstractFloat,
                           AT4<:AbstractArray}
  f::AT4
  fTemp::AT4
end

function Cache(backend,FT,IntMethod::LSRungeKuttaMethod,FE,M,nz,NumV)
  NumG = FE.NumG
  NumI = FE.NumI
  nStage = IntMethod.nStage
  f = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV)
  fTemp = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV)
  return CacheLSRKStruct{FT,
                     typeof(f)}(
    f,
    fTemp,
  )
end
function TimeIntegration!(LSRK::LSRungeKuttaMethod,V,dt,Fcn,Aux,Jac,FE,Metric,Phys,Cache,JCache,Exchange,
  Global,Param,DiscType)

  f = Cache.f
  fTemp = Cache.fTemp

  dtau, = dt
  FcnFull, = Fcn

  FcnFull(f,V,FE,Metric,Phys,Aux,Exchange,Global,DiscType)
  @. fTemp = f
  @. V += (dtau * LSRK.b[1]) * f
  @inbounds for iStage = 2 : LSRK.nStage
    FcnFull(f,V,FE,Metric,Phys,Aux,Exchange,Global,DiscType)
    @. fTemp = f - LSRK.a[iStage] * fTemp
    @. V += (dtau * LSRK.b[iStage]) * fTemp
  end
end

function TimeIntegrationFast!(LSRK::LSRungeKuttaMethod,V,dt,Fcn,FSlow,Aux,FE,Metric,Phys,Cache,JCache,Exchange,
  Global,Param,DiscType)

  f = Cache.CacheFast.f
  fTemp = Cache.CacheFast.fTemp
  Vn = Cache.Vn
  @views VnI = Vn[:,:,1:size(V,3),1:size(V,4)]

  @. VnI = V
  Fcn(f,Vn,FE,Metric,Phys,Aux,Exchange,Global,DiscType)
  @. fTemp = f + FSlow
  @. VnI += (dt * LSRK.b[1]) * fTemp
  @inbounds for iStage = 2 : LSRK.nStage
    Fcn(f,Vn,FE,Metric,Phys,Aux,Exchange,Global,DiscType)
    @. fTemp = f + FSlow - LSRK.a[iStage] * fTemp
    @. VnI += (dt * LSRK.b[iStage]) * fTemp
  end
  @. V = VnI
end

