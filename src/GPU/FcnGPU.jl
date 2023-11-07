function FcnAdvectionGPU!(F,U,time,FE,Metric,Phys,Cache,Global,Param,Profile)

  backend = get_backend(F)
  FT = eltype(F)
  Glob = FE.Glob
  DS = FE.DS
  DW = FE.DW
  M = FE.M
  dXdxI = Metric.dXdxI
  X = Metric.X
  J = Metric.J
  N = FE.OrdPoly+1
  Nz = size(F,1)
  NF = size(Glob,3)
  Koeff = Global.Model.HyperDDiv
  Temp1 = Cache.Temp1


# State vector
  @views Rho = U[:,:,1]
  @views u = U[:,:,2]
  @views v = U[:,:,3]
  @views w = U[:,:,4]
# Cache
  @views CacheF = Temp1[:,:,1:5]
# Ranges
  NzG = min(div(1024,N*N),Nz)
  group = (N, N, NzG, 1)
  ndrange = (N, N, Nz, NF)

  KuvwFunCKernel! = uvwFunCKernel!(backend, group)
  KDivRhoGradKernel! = DivRhoGradKernel!(backend, group)
  KHyperViscKernel! = HyperViscKernel!(backend, group)
  KDivRhoTrUpwind3Kernel! = DivRhoTrUpwind3Kernel!(backend, group)

  KuvwFunCKernel!(Profile,u,v,w,time,Glob,X,Param,Phys,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)

  @. CacheF = 0
  KHyperViscKernel!(CacheF,U,DS,DW,dXdxI,J,M,Glob,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)
  @. F = 0
  KDivRhoTrUpwind3Kernel!(F,U,CacheF,DS,DW,dXdxI,J,M,Glob,Koeff,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)
end

function FcnGPU!(F,U,FE,Metric,Phys,Cache,Exchange,Global,Param,DiscType)

  backend = get_backend(F)
  FT = eltype(F)
  @. F = 0
  Glob = FE.Glob
  DS = FE.DS
  DW = FE.DW
  M = FE.M
  dXdxI = Metric.dXdxI
  X = Metric.X
  J = Metric.J
  lat = Metric.lat  
  dz = Metric.dz  
  DoF = FE.DoF
  N = size(FE.DS,1)
  Nz = size(F,1)
  NDoF = size(F,2)
  NF = size(Glob,2)
  NumV  = Global.Model.NumV 
  NumTr  = Global.Model.NumTr 
  Koeff = Global.Model.HyperDDiv
  Temp1 = Cache.Temp1
  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU
  Force = Global.Model.Force
  MicrophysicsSource = Global.Model.MicrophysicsSource

  KoeffCurl = Global.Model.HyperDCurl
  KoeffGrad = Global.Model.HyperDGrad
  KoeffDiv = Global.Model.HyperDDiv

# State vector
  @views Rho = U[:,:,1]
  @views u = U[:,:,2]
  @views v = U[:,:,3]
  @views w = U[:,:,4]
  @views RhoTr = U[:,:,5]
# Tendency
  @views FRho = F[:,:,1]
  @views FRhoTr = F[:,:,5]
# Cache
  @views CacheF = Temp1[:,:,1:6]
  @views CacheFF = Temp1[:,:,1:6+NumTr]
  @views p = Cache.AuxG[:,:,1]
  KV = Cache.KV
  TSurf = Cache.TSurf
  RhoVSurf = Cache.RhoVSurf
  uStar = Cache.uStar
  CT = Cache.CT
  CH = Cache.CH
# Ranges
  NzG = min(div(NumberThreadGPU,N*N),Nz)
  group = (N, N, NzG, 1)
  ndrange = (N, N, Nz, NF)
  groupTr = group
  ndrangeTr = ndrange
  NDoFG = min(div(NumberThreadGPU,Nz),NDoF)
  groupG = (Nz, NDoFG)  
  ndrangeG = (Nz, NDoF)  

  KRhoGradKinKernel! = RhoGradKinKernel!(backend,group)
  KGradKernel! = GradKernel!(backend,group)
  KDivRhoGradKernel! = DivRhoGradKernel!(backend, group)
  KHyperViscKernel! = HyperViscKernel!(backend, group)
  KHyperViscKoeffKernel! = HyperViscKoeffKernel!(backend, group)
  KDivRhoThUpwind3Kernel! = DivRhoThUpwind3Kernel!(backend, group)
  KMomentumCoriolisKernel! = MomentumCoriolisKernel!(backend, group)
  KHyperViscTracerKernel! = HyperViscTracerKernel!(backend, groupTr)
  KHyperViscTracerKoeffKernel! = HyperViscTracerKoeffKernel!(backend, groupTr)
  KDivRhoTrUpwind3Kernel! = DivRhoTrUpwind3Kernel!(backend, groupTr)
# KMomentumKernel! = MomentumKernel!(backend, group)

  @. Temp1 = FT(0)
  @views MRho = CacheF[:,:,6]
  KHyperViscKernel!(CacheF,MRho,U,DS,DW,dXdxI,J,M,Glob,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)
  for iT = 1 : NumTr
    @views CacheTr = Temp1[:,:,iT + 6]  
    KHyperViscTracerKernel!(CacheTr,U[:,:,iT+NumV],Rho,DS,DW,dXdxI,J,M,Glob,ndrange=ndrangeTr)
    KernelAbstractions.synchronize(backend)
  end  

  Parallels.ExchangeData3DSendGPU(CacheFF,Exchange)
  Parallels.ExchangeData3DRecvGPU!(CacheFF,Exchange)
  KernelAbstractions.synchronize(backend)

  @. F = 0
  KHyperViscKoeffKernel!(F,U,CacheF,DS,DW,dXdxI,J,M,Glob,KoeffCurl,KoeffGrad,KoeffDiv,ndrange=ndrange)
  for iT = 1 : NumTr
    @views CacheTr = Temp1[:,:,iT + 6]  
    @views KHyperViscTracerKoeffKernel!(F[:,:,iT+NumV],CacheTr,Rho,DS,DW,dXdxI,J,M,Glob,
      KoeffDiv,ndrange=ndrangeTr)
    KernelAbstractions.synchronize(backend)
  end  


  KGradKernel!(F,U,p,DS,dXdxI,J,M,MRho,Glob,Phys,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)

  KMomentumCoriolisKernel!(F,U,DS,dXdxI,J,X,MRho,M,Glob,Phys,ndrange=ndrange)
# KMomentumKernel!(F,U,DS,dXdxI,MRho,M,Glob,Phys,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)

  KDivRhoThUpwind3Kernel!(F,U,DS,dXdxI,J,M,Glob,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)

  for iT = 1 : NumTr
    @views KDivRhoTrUpwind3Kernel!(F[:,:,iT+NumV],U[:,:,iT+NumV],U,DS,
      dXdxI,J,M,Glob,ndrange=ndrangeTr)
    KernelAbstractions.synchronize(backend)
  end  

  if Global.Model.VerticalDiffusion
    KVerticalDiffusionScalarKernel! = VerticalDiffusionScalarKernel!(backend,groupTr)  
    @views KVerticalDiffusionScalarKernel!(F[:,:,5],U[:,:,5],Rho,KV,
      dXdxI,J,M,Glob,ndrange=ndrangeTr)
    KernelAbstractions.synchronize(backend)
    for iT = 1 : NumTr
      @views KVerticalDiffusionScalarKernel!(F[:,:,iT+NumV],U[:,:,iT+NumV],Rho,KV,
        dXdxI,J,M,Glob,ndrange=ndrangeTr)
      KernelAbstractions.synchronize(backend)
    end  
  end    

  if Global.Model.SurfaceFlux
    NFG = min(div(NumberThreadGPU,N*N),NF)
    groupS = (N * N, NFG)
    ndrangeS = (N * N, NF)  
    KSurfaceFluxScalarsKernel = SurfaceFluxScalarsKernel(backend,groupS)  
    KSurfaceFluxScalarsKernel(F,U,p,TSurf,RhoVSurf,uStar,CT,CH,dXdxI,Glob,M,Phys,ndrange=ndrangeS)
    KernelAbstractions.synchronize(backend)
  end  
      
  Parallels.ExchangeData3DSendGPU(F,Exchange)
  Parallels.ExchangeData3DRecvGPU!(F,Exchange)
  KernelAbstractions.synchronize(backend)

# if Global.Model.VerticalDiffusion
#   KVerticalDiffusionCScalarKernel! = VerticalDiffusionScalarKernel!(backend, groupG)
#   @views KVerticalDiffusionCScalarKernel!(F[:,:,5],U[:,:,5],Rho,KV,
#     dz,ndrange=ndrangeG)
#   for iT = 1 : NumTr
#     @views KVerticalDiffusionCScalarKernel!(F[:,:,iT+NumV],U[:,:,iT+NumV],Rho,KV,
#       dz,ndrange=ndrangeG)
#   end  
#   KernelAbstractions.synchronize(backend)
# end    


  if Global.Model.Forcing
    KForceKernel! = ForceKernel!(backend, groupG)
    KForceKernel!(Force,F,U,p,lat,ndrange=ndrangeG)  
    KernelAbstractions.synchronize(backend)
  end  

  if Global.Model.Microphysics
    KMicrophysicsKernel! = MicrophysicsKernel!(backend, groupG)
    KMicrophysicsKernel!(MicrophysicsSource,F,U,p,ndrange=ndrangeG)
    KernelAbstractions.synchronize(backend)
  end

end


