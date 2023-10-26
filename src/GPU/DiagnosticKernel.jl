@kernel function PressureKernel!(p,@Const(U),Phys)
  Iz,IC = @index(Global, NTuple)

  NumG = @uniform @ndrange()[2]

  if IC <= NumG
    @inbounds p[Iz,IC] = PressureGPU(view(U,Iz,IC,:),Phys)
  end
end

@inline function PressureGPU(U,Phys)
  @inbounds Phys.p0 * fast_powGPU(Phys.Rd * U[5] / Phys.p0, 1 / (1 - Phys.kappa))
end

# we may be hitting a slow path:
# https://stackoverflow.com/questions/14687665/very-slow-stdpow-for-bases-very-close-to-1
fast_powGPU(x::FT, y::FT) where {FT <: AbstractFloat} = exp(y * log(x))

@kernel function uStarCoefficientKernel!(uStar,@Const(U),@Const(dXdxI),@Const(nS),@Const(Glob))
  ID,IF = @index(Global, NTuple)

  NumF = @uniform @ndrange()[2]

  if IF <= NumF
    ind = Glob[ID,IF]
    @inbounds @views uStar[ID,IF] = uStarCoefficientGPU(U[ind,2],U[ind,3],U[ind,4],dXdxI[:,ID,IF],nS[:,ID,IF])
  end
end

@inline function uStarCoefficientGPU(v1,v2,w,dXdxI,nS)
# Computation norm_v_a
# |v_a| = |v - n(n*v)| = sqrt(v*v -(n*v)^2)
  @inbounds wS = -(dXdxI[3,1]* v1 + dXdxI[3,2] * v2) / dXdxI[3,3]
  w = 0.5 * (wS + w)
  @inbounds nU = nS[1] * v1 + nS[2] * v2 + nS[3] * w
  @inbounds sqrt((v1 - nS[1] * nU) * (v1 - nS[1] * nU) + 
    (v2 - nS[2] * nU) * (v2 - nS[2] * nU) +
    (w - nS[3] * nU) * (w - nS[3] * nU))
end

@kernel function EddyCoefficientKernel!(K,@Const(uStar),@Const(p),@Const(dz),@Const(Glob),Param)
  ID,Iz,IF = @index(Global, NTuple)

  NumF = @uniform @ndrange()[2]

  if Iz <= Nz && IF <= NumF
    ind = Glob[ID,IF]
    @inbounds @views K[ID,iZ,IF] = EddyCoefficientGPU(uStar[ID,IF],p[Iz,ind],dz,Param)
  end
end

@inline function EddyCoefficientGPU(uStar,p,dz,Param)
  K = Param.CE * uStar * dz / 2
  if p < Param.p_pbl
    dpR = (Param.p_pbl - p) / Param.p_strato  
    K = K * exp(-dpR * dpR)  
  end
end

function FcnPrepareGPU!(U,FE,Metric,Phys,Cache,Exchange,Global,Param,DiscType)
  backend = get_backend(U)
  dXdxI = Metric.dXdxI
  nS = Metric.nS
  FT = eltype(U)
  N = size(FE.DS,1)
  Nz = size(U,1)
  NumG = size(U,2)
  Glob = FE.Glob
  NumF = size(Glob,2)

  NG = min(div(512,Nz),NumG)
  group = (Nz, NG)
  ndrange = (Nz, NumG)
  NF = min(div(512,N*N),NumF)
  groupS = (N * N, NF)
  ndrangeS = (N * N, NumF)
  groupK = (N * N, NF)
  ndrangeK = (N * N, NumF)
  @views p = Cache.AuxG[:,:,1]
  uStar = Cache.uStar

  KPressureKernel! = PressureKernel!(backend,group)
  KPressureKernel!(p,U,Phys,ndrange=ndrange)

  if Global.Model.VerticalDiffusion
    KuStarCoefficientKernel! = uStarCoefficientKernel!(backend,groupS)
    KEddyCoefficientKernel! = EddyCoefficientKernel!(backend,groupK)
    @views KuStarCoefficientKernel!(uStar,U[1,:,:],dXdxI[3,:,1,:,1,:],nS,Glob,ndrange=ndrangeS)
  end   

  KernelAbstractions.synchronize(backend)
end