mutable struct CacheLinIMEXStruct{FT<:AbstractFloat,
                           AT4<:AbstractArray,
                           AT5<:AbstractArray}
  Vn::AT4
  f::AT5
  fV::AT4
  Ymyn::AT5
end

function Cache(backend,FT,IntMethod::LinIMEXMethod,FE,M,nz,NumV)
  NumG = FE.NumG
  NumI = FE.NumI
  nStage = IntMethod.nStage
  Vn = KernelAbstractions.zeros(backend,FT,M,nz,NumG,NumV)
  fV = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV)
  Ymyn = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV,nStage)
  f = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV,nStage)
  return CacheLinIMEXStruct{FT,
                     typeof(Vn),
                     typeof(f)}(
    Vn,
    f,
    fV,
    Ymyn,
  )
end

#function LinIMEXSchur!(V,dt,Fcn,Jac,CG,Global,Param)
function TimeIntegration!(LinIMEX::LinIMEXMethod,V,dt,Fcn,Aux,Jac,FE,Metric,Phys,Cache,JCache,Exchange,
  Global,Param,DiscType)

  nStage = LinIMEX.nStage
  gammaD = LinIMEX.AI[2,2]
  f = Cache.f
  fV = Cache.fV
  Vn = Cache.Vn
  Ymyn = Cache.Ymyn
  @views VnI = Vn[:,:,1:size(fV,3),1:size(fV,4)]
  dtau, = dt
  FcnFull, = Fcn

  @. VnI = V
  Jac(V,dtau*gammaD,FE,Metric,Phys,Aux,JCache,Global,DiscType)

  @. Ymyn[:,:,:,:,1] = 0
  @inbounds for iStage = 2 : nStage
    @views FcnFull(f[:,:,:,:,iStage],Vn,FE,Metric,Phys,Aux,Exchange,Global,DiscType)

    @views @. fV = LinIMEX.D[iStage,1] * f[:,:,:,:,1]
    @inbounds for jStage = 2 : iStage - 1
      @views @. fV += LinIMEX.D[iStage,jStage] * f[:,:,:,:,jStage] +
        (LinIMEX.E[iStage,jStage] / dtau) * Ymyn[:,:,:,:,jStage] 
    end
    @views Solve!(Ymyn[:,:,:,:,iStage],fV,JCache,dtau*gammaD,FE,Metric,Global,DiscType)
    @views @. VnI = Ymyn[:,:,:,:,iStage] + V
  end

  if LinIMEX.d[nStage] != 0.0
    @views FcnFull(f[:,:,:,:,nStage],Vn,FE,Metric,Phys,Aux,Exchange,Global,DiscType)
    @views @. V += dtau * LinIMEX.d[nStage] * f[:,:,:,:,nStage]
  end  
  @inbounds for iStage = 1 : nStage - 1
    @views @. V += (dtau * LinIMEX.d[iStage]) * f[:,:,:,:,iStage] +
      LinIMEX.e[iStage+1] * Ymyn[:,:,:,:,iStage+1]
  end
end
