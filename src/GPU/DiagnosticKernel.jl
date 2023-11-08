#@kernel function SurfFlux()
#  RhoD = Rho - RhoV
#      Rm = Phys.Rd * RhoD + Phys.Rv * RhoV
#      Cpml = Phys.Cpd * RhoD + Phys.Cpv * RhoV
#      T = p / Rm
#      p_vs = fpvs(TSurf,Phys.  )
#
#      RhoVSurface = p_vs / (Phys.Rv * TSurf)
#      LatFlux = -CE * uStar * (RhoV - RhoVSurface)
#      SensFlux = -CH * uStar *  (T - TSurf)
#
#end

@kernel function PressureKernel!(Pressure,p,@Const(U))
  Iz,IC = @index(Global, NTuple)

  NumG = @uniform @ndrange()[2]

  if IC <= NumG
    @inbounds p[Iz,IC] = Pressure(view(U,Iz,IC,:))
  end
end

@kernel function PressureCKernel!(Pressure,p,@Const(U),@Const(Glob),@Const(JJ),@Const(M))
  ID,Iz,IF = @index(Global, NTuple)

  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  if Iz <= Nz && IF <= NF
    @inbounds ind = Glob[ID,IF]  
    @inbounds @atomic p[Iz,ind] += Pressure(view(U,Iz,ind,:)) * JJ[ID,Iz,IF] / M[Iz,ind]
  end
end

#@inline function PressureGPU(U,Phys)
#  @inbounds Phys.p0 * fast_powGPU(Phys.Rd * U[5] / Phys.p0, 1 / (1 - Phys.kappa))
#end

#=
@kernel function uStarCoefficientKernel!(uStar,@Const(U),@Const(dXdxI),@Const(nS),@Const(Glob))
  ID,IF = @index(Global, NTuple)

  NumF = @uniform @ndrange()[2]

  if IF <= NumF
    ind = Glob[ID,IF]
    @inbounds uStar[ID,IF] = uStarCoefficientGPU(U[1,ind,2],U[1,ind,3],U[1,ind,4],view(dXdxI,3,:,1,ID,1,IF),view(nS,:,ID,IF))
  end
end

@kernel function SurfaceKoeffKernel!(uStar,CT,CH,@Const(U),@Const(dXdxI),@Const(nS),@Const(Glob),Surface)
  ID,IF = @index(Global, NTuple)

  NumF = @uniform @ndrange()[2]

  if IF <= NumF
    ind = Glob[ID,IF]
    @inbounds uStar[ID,IF], CT[ID,IF], CH[ID,IF] = Surface(view(U,1,ind,:),view(dXdxI,3,:,1,ID,1,IF),view(nS,:,ID,IF))
  end
end
=#

@kernel function SurfaceKernel!(TSurf,RhoVSurf,uStar,CT,CH,@Const(U),@Const(p),@Const(X),
  @Const(dXdxI),@Const(nS),@Const(Glob),SurfaceValues,SurfaceData)

  ID,IF = @index(Global, NTuple)

  NumF = @uniform @ndrange()[2]

  if IF <= NumF
    @inbounds ind = Glob[ID,IF]
    @inbounds x1 = X[ID,1,1,1,IF]
    @inbounds x2 = X[ID,1,2,1,IF]
    @inbounds x3 = X[ID,1,3,1,IF] 
    xS = SVector{3}(x1, x2 ,x3)
    @inbounds TSurf[ID,IF], RhoVSurf[ID,IF]  = SurfaceValues(xS,view(U,1,ind,:),p[1,ind])
    @inbounds uStar[ID,IF], CT[ID,IF], CH[ID,IF] = SurfaceData(view(U,1,ind,:),p[1,ind],
     view(dXdxI,3,:,1,ID,1,IF),view(nS,ID,:,IF))
  end
end

@kernel function EddyCoefficientKernel!(Eddy,K,@Const(Rho),@Const(uStar),@Const(p),@Const(dz),@Const(Glob))
  ID,Iz,IF = @index(Global, NTuple)

  Nz = @uniform @ndrange()[2]
  NumF = @uniform @ndrange()[3]

  if Iz <= Nz && IF <= NumF
    @inbounds ind = Glob[ID,IF]
    @inbounds @views K[ID,Iz,IF] = Eddy(uStar[ID,IF],p[Iz,ind],dz[1,ind]) * Rho[Iz,ind]
  end
end

@kernel function EddyCoefficientCKernel!(Eddy,K,@Const(Rho),@Const(uStar),@Const(p),@Const(dz),@Const(Glob))
  ID,Iz,IF = @index(Global, NTuple)

  Nz = @uniform @ndrange()[2]
  NumF = @uniform @ndrange()[3]

  if Iz <= Nz && IF <= NumF
    @inbounds ind = Glob[ID,IF]
    @inbounds p = Pressure(view(U,Iz,ind,:)) 
    @inbounds K[Iz,ind] += Eddy(uStar[ID,IF],p,dz[1,ind]) * Rho[Iz,ind] * JJ[ID,Iz,IF] / M[Iz,ind]
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
  X = Metric.X
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
  TSurf = Cache.TSurf
  RhoVSurf = Cache.RhoVSurf
  uStar = Cache.uStar
  KV = Cache.KV
  CT = Cache.CT
  CH = Cache.CH
  @views Rho = U[:,:,1]
# LatFlux = Cache.LatFlux
# SensFlux = Cache.SensFlux
  dz = Metric.dz
  Eddy = Global.Model.Eddy
  Pressure = Global.Model.Pressure

  KPressureKernel! = PressureKernel!(backend,group)
  KPressureKernel!(Pressure,p,U,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)

  if Global.Model.SurfaceFlux
    NFG = min(div(NumberThreadGPU,N*N),NumF)
    groupS = (N * N, NFG)
    ndrangeS = (N * N, NumF)
    KSurfaceKernel! = SurfaceKernel!(backend,groupS) 
    KSurfaceKernel!(TSurf,RhoVSurf,uStar,CT,CH,U,p,X,dXdxI,nS,Glob,Global.Model.SurfaceValues,
      Global.Model.SurfaceData,ndrange=ndrangeS)
    KernelAbstractions.synchronize(backend)
  end  
  if Global.Model.VerticalDiffusion 
    if Global.Model.SurfaceFlux  
      KEddyCoefficientKernel! = EddyCoefficientKernel!(backend,groupK)
      KEddyCoefficientKernel!(Eddy,KV,Rho,uStar,p,dz,Glob,ndrange=ndrangeK)
      KernelAbstractions.synchronize(backend)
    else    
      NFG = min(div(NumberThreadGPU,N*N),NumF)
      groupS = (N * N, NFG)
      ndrangeS = (N * N, NumF)
      KSurfaceKernel! = SurfaceKernel!(backend,groupS) 
      SurfaceValues = Global.Model.SurfaceValues
      SurfaceData = Global.Model.SurfaceData
      @show "vor KSurfaceKernel"
      KSurfaceKernel!(TSurf,RhoVSurf,uStar,CT,CH,U,p,X,dXdxI,nS,Glob,SurfaceValues,
        SurfaceData,ndrange=ndrangeS)
      @show "nack KSurfaceKernel"
      KernelAbstractions.synchronize(backend)
      KEddyCoefficientKernel! = EddyCoefficientKernel!(backend,groupK)
      KEddyCoefficientKernel!(Eddy,KV,Rho,uStar,p,dz,Glob,ndrange=ndrangeK)
      KernelAbstractions.synchronize(backend)
    end
  end
end

