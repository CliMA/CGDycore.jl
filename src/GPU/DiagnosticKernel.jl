@kernel function PressureKernel!(Pressure,p,@Const(U))
  Iz,IC = @index(Global, NTuple)

  NumG = @uniform @ndrange()[2]

  if IC <= NumG
    @inbounds p[Iz,IC] = Pressure(view(U,Iz,IC,:))
  end
end

#@inline function PressureGPU(U,Phys)
#  @inbounds Phys.p0 * fast_powGPU(Phys.Rd * U[5] / Phys.p0, 1 / (1 - Phys.kappa))
#end


@kernel function uStarCoefficientKernel!(uStar,@Const(U),@Const(dXdxI),@Const(nS),@Const(Glob))
  ID,IF = @index(Global, NTuple)

  NumF = @uniform @ndrange()[2]

  if IF <= NumF
    ind = Glob[ID,IF]
    @inbounds uStar[ID,IF] = uStarCoefficientGPU(U[1,ind,2],U[1,ind,3],U[1,ind,4],view(dXdxI,3,:,1,ID,1,IF),view(nS,:,ID,IF))
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

@kernel function EddyCoefficientKernel!(Eddy,K,@Const(uStar),@Const(p),@Const(dz),@Const(Glob))
  ID,Iz,IF = @index(Global, NTuple)

  Nz = @uniform @ndrange()[2]
  NumF = @uniform @ndrange()[3]

  if Iz <= Nz && IF <= NumF
    @inbounds ind = Glob[ID,IF]
    @inbounds @views K[ID,Iz,IF] = Eddy(uStar[ID,IF],p[Iz,ind],dz[1,ind])
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
  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU

  NG = min(div(NumberThreadGPU,Nz),NumG)
  group = (Nz, NG)
  ndrange = (Nz, NumG)
  NF = min(div(NumberThreadGPU,N*N),NumF)
  groupS = (N * N, NF)
  ndrangeS = (N * N, NumF)
  NzG = min(div(NumberThreadGPU,N*N),Nz)
  groupK = (N * N, NzG, 1)
  ndrangeK = (N * N, Nz, NumF)
  @views p = Cache.AuxG[:,:,1]
  uStar = Cache.uStar
  KV = Cache.KV
  dz = Metric.dz
  Eddy = Global.Model.Eddy
  Pressure = Global.Model.Pressure

  KPressureKernel! = PressureKernel!(backend,group)
  KPressureKernel!(Pressure,p,U,ndrange=ndrange)

  if Global.Model.VerticalDiffusion
    KuStarCoefficientKernel! = uStarCoefficientKernel!(backend,groupS)
    KEddyCoefficientKernel! = EddyCoefficientKernel!(backend,groupK)
    KuStarCoefficientKernel!(uStar,U,dXdxI,nS,Glob,ndrange=ndrangeS)
    KEddyCoefficientKernel!(Eddy,KV,uStar,p,dz,Glob,ndrange=ndrangeK)
  end   

  KernelAbstractions.synchronize(backend)
end

