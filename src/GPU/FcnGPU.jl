function FcnAdvectionGPU!(F,U,time,FE,Metric,Phys,Cache,Exchange,Global,Param,Profile)

  backend = get_backend(F)
  FT = eltype(F)
  dtau = Global.TimeStepper.dtauStage
  Glob = FE.Glob
  DS = FE.DS
  DW = FE.DW
  M = FE.M
  DoF = FE.DoF
  Stencil = FE.Stencil
  dXdxI = Metric.dXdxI
  X = Metric.X
  J = Metric.J
  JC = Metric.JC
  JCW = Metric.JCW
  N = FE.OrdPoly+1
  ww = FE.w
  Nz = size(F,1)
  NF = size(Glob,2)
  NDoF = size(Glob,1)
  KoeffDiv = Global.Model.HyperDDiv
  NumV  = Global.Model.NumV
  NumTr  = Global.Model.NumTr
  Temp1 = Cache.Temp1
  @views qMin = Cache.qMin[:,:,1:NumTr]
  @views qMax = Cache.qMax[:,:,1:NumTr]


# State vector
  @views Rho = U[:,:,1]
  @views u = U[:,:,2]
  @views v = U[:,:,3]
  @views w = U[:,:,4]
# Cache
  @views CacheF = Temp1[:,:,1:5]
  @views CacheTr = Temp1[:,:,1]
# Ranges
  NzG = min(div(1024,N*N),Nz)
  group = (N, N, NzG, 1)
  ndrange = (N, N, Nz, NF)
  groupC = (N* N, NzG, 1)
  ndrangeC = (N* N, Nz, NF)
  groupLim = (10, 1)
  ndrangeLim = (Nz, NF)
  NFG = min(div(512,Nz),NF)
  groupL = (Nz, NFG, 1)
  ndrangeL = (Nz, NF, NumTr)

  KLimitKernel! = LimitKernel!(backend, groupL)
  KuvwFunCKernel! = uvwFunCKernel!(backend, groupC)
  KDivRhoKernel! = DivRhoKernel!(backend, group)
  KHyperViscTracerKernel! = HyperViscTracerKernel!(backend, group)
  KHyperViscTracerKoeffKernel! = HyperViscTracerKoeffKernel!(backend, group)
  KDivRhoTrUpwind3Kernel! = DivRhoTrUpwind3Kernel!(backend, group)
  KDivRhoTrUpwind3LimKernel! = DivRhoTrUpwind3LimKernel!(backend, group)
  KDivRhoTrViscUpwind3LimKernel! = DivRhoTrViscUpwind3LimKernel!(backend, group)

  if Global.Model.HorLimit
    @views KLimitKernel!(DoF,qMin,qMax,U[:,:,NumV+1:NumV+NumTr],Rho,Glob,ndrange=ndrangeL)
    KernelAbstractions.synchronize(backend)
    Parallels.ExchangeDataFSendGPU(qMin,qMax,Exchange)
    Parallels.ExchangeDataFRecvGPU!(qMin,qMax,Exchange)  
  end


# Velocity 
  KuvwFunCKernel!(Profile,u,v,w,time,Glob,X,Param,Phys,ndrange=ndrangeC)
  KernelAbstractions.synchronize(backend)

# Hyperviscosity Part 1
  if ~Global.Model.HorLimit
     KHyperViscTracerKernel!(CacheTr,U[:,:,1+NumV],Rho,DS,DW,dXdxI,J,M,Glob,ndrange=ndrange)
     KernelAbstractions.synchronize(backend)

#   Data exchange  
    Parallels.ExchangeData3DSendGPU(CacheTr,Exchange)
    Parallels.ExchangeData3DRecvGPU!(CacheTr,Exchange)
  end  

  F .= FT(0)

  KDivRhoKernel!(F,U,DS,dXdxI,J,M,Glob,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)  

  if Global.Model.HorLimit
    @views KDivRhoTrUpwind3LimKernel!(F[:,:,1+NumV],U[:,:,1+NumV],U,DS,
      dXdxI,J,M,Glob,dtau,ww,qMin[:,:,1],qMax[:,:,1],Stencil,ndrange=ndrange)
    KernelAbstractions.synchronize(backend)  
  else
    @views KHyperViscTracerKoeffKernel!(F[:,:,1+NumV],CacheTr,Rho,DS,DW,dXdxI,J,M,Glob,
      KoeffDiv,ndrange=ndrange)
    KernelAbstractions.synchronize(backend)

    @views KDivRhoTrUpwind3Kernel!(F[:,:,1+NumV],U[:,:,1+NumV],U,DS,
      dXdxI,J,M,Glob,ndrange=ndrange)
    KernelAbstractions.synchronize(backend)  
  end

# Data exchange  
  @views Parallels.ExchangeData3DSendGPU(F[:,:,1:1+NumV],Exchange)
  @views Parallels.ExchangeData3DRecvGPU!(F[:,:,1:1+NumV],Exchange)
end

function FcnGPU!(F,U,FE,Metric,Phys,Cache,Exchange,Global,Param,Equation::Models.CompressibleShallow)

  backend = get_backend(F)
  FT = eltype(F)
  Glob = FE.Glob
  DS = FE.DS
  DW = FE.DW
  M = FE.M
  MW = FE.MW
  dXdxI = Metric.dXdxI
  nS = Metric.nS
  X = Metric.X
  J = Metric.J
  NF = Global.Grid.NumFaces
  NBF = Global.Grid.NumBoundaryFaces
  @views dXdxI_B = dXdxI[:,:,:,:,:,1:NBF]
  @views nS_B = nS[:,:,1:NBF]
  @views J_B = J[:,:,:,1:NBF]
  @views X_B = X[:,:,:,:,1:NBF]
  @views Glob_B = Glob[:,1:NBF]

  @views dXdxI_I = dXdxI[:,:,:,:,:,NBF+1:NF]
  @views J_I = J[:,:,:,NBF+1:NF]
  @views X_I = X[:,:,:,:,NBF+1:NF]
  @views nS_I = nS[:,:,NBF+1:NF]
  @views Glob_I = Glob[:,NBF+1:NF]
  lat = Metric.lat  
  dz = Metric.dz  
  zP = Metric.zP  
  DoF = FE.DoF
  N = size(FE.DS,1)
  Nz = size(F,1)
  NDoF = size(F,2)
  NumV  = Global.Model.NumV 
  NumTr  = Global.Model.NumTr 
  Koeff = Global.Model.HyperDDiv
  Temp1 = Cache.Temp1
  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU
  Force = Global.Model.Force
  Damp = Global.Model.Damp
  MicrophysicsSource = Global.Model.MicrophysicsSource
  CoriolisFun = Global.Model.CoriolisFun
  GravitationFun = Global.Model.GravitationFun

  KoeffCurl = Global.Model.HyperDCurl
  KoeffGrad = Global.Model.HyperDGrad
  KoeffDiv = Global.Model.HyperDDiv
  KoeffDivW = Global.Model.HyperDDivW

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
  @views CacheFF = Temp1[:,:,1:6+NumTr+1]
  @views Cachew = Temp1[:,:,6 + 1 + NumTr]  
  @views p = Cache.AuxG[:,:,1]
# Ranges
  NzG = min(div(NumberThreadGPU,N*N),Nz)
  group = (N, N, NzG, 1)
  ndrange = (N, N, Nz, NF)
  ndrangeB = (N, N, Nz, NBF)
  ndrangeI = (N, N, Nz, NF-NBF)
  groupTr = group
  ndrangeTr = ndrange
  NDoFG = min(div(NumberThreadGPU,Nz),NDoF)
  groupG = (Nz, NDoFG)  
  ndrangeG = (Nz, NDoF)  
  NzG = min(div(NumberThreadGPU,N*N),Nz-1)
  groupw = (N, N, NzG, 1)
  ndrangewB = (Nz-1, NBF)  
  ndrangewI = (Nz-1, NF-NBF)  

  KRhoGradKinKernel! = RhoGradKinKernel!(backend,group)
  KGradKernel! = GradKernel!(backend,group)
  KHyperViscKernel! = HyperViscKernel!(backend, group)
  KHyperViscKoeffKernel! = HyperViscKoeffKernel!(backend, group)
  KDivRhoThUpwind3Kernel! = DivRhoThUpwind3Kernel!(backend, group)
  KMomentumCoriolisKernel! = MomentumVectorInvariantCoriolisKernel!(backend, group)
  KHyperViscTracerKernel! = HyperViscTracerKernel!(backend, groupTr)
  KHyperViscTracerKoeffKernel! = HyperViscTracerKoeffKernel!(backend, groupTr)
  KDivRhoTrUpwind3Kernel! = DivRhoTrUpwind3Kernel!(backend, groupTr)


####
# First phase  
####
  Temp1 .= FT(0)
  @views MRho = CacheF[:,:,6]
  KHyperViscKernel!(CacheF,MRho,U,DS,DW,dXdxI,J,M,Glob,ndrange=ndrangeB)
  KernelAbstractions.synchronize(backend)
  for iT = 1 : NumTr
    @views CacheTr = Temp1[:,:,iT + 6]  
    @views KHyperViscTracerKernel!(CacheTr,U[:,:,iT+NumV],Rho,DS,DW,dXdxI,J,M,Glob,ndrange=ndrangeB)
    KernelAbstractions.synchronize(backend)
  end  
  if KoeffDivW > 0
    KHyperViscWKernel! = HyperViscWKernel!(backend, groupTr)
    @views KHyperViscWKernel!(Cachew,U[:,:,4],DS,DW,dXdxI,J,MW,Glob,ndrange=ndrangeB)
    KernelAbstractions.synchronize(backend)
  end  
  Parallels.ExchangeData3DSendGPU(CacheFF,Exchange)

  KHyperViscKernel!(CacheF,MRho,U,DS,DW,dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeI)
  KernelAbstractions.synchronize(backend)
  for iT = 1 : NumTr
    @views CacheTr = Temp1[:,:,iT + 6]  
    @views KHyperViscTracerKernel!(CacheTr,U[:,:,iT+NumV],Rho,DS,DW,dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeI)
    KernelAbstractions.synchronize(backend)
  end  
  if KoeffDivW > 0
    @views KHyperViscWKernel!(Cachew,U[:,:,4],DS,DW,dXdxI_I,J_I,MW,Glob_I,ndrange=ndrangeI)
    KernelAbstractions.synchronize(backend)
  end  

  Parallels.ExchangeData3DRecvGPU!(CacheFF,Exchange)
  KernelAbstractions.synchronize(backend)

####
# Second phase  
####

  F .= FT(0)
  KHyperViscKoeffKernel!(F,U,CacheF,DS,DW,dXdxI,J,M,Glob,KoeffCurl,KoeffGrad,KoeffDiv,ndrange=ndrangeB)
  KernelAbstractions.synchronize(backend)
  for iT = 1 : NumTr
    @views CacheTr = Temp1[:,:,iT + 6]  
    @views KHyperViscTracerKoeffKernel!(F[:,:,iT+NumV],CacheTr,Rho,DS,DW,dXdxI,J,M,Glob,
      KoeffDiv,ndrange=ndrangeB)
    KernelAbstractions.synchronize(backend)
  end  
  if KoeffDivW > 0
    KHyperViscWKoeffKernel! = HyperViscWKoeffKernel!(backend, groupTr)
    @views KHyperViscWKoeffKernel!(F[:,:,4],Cachew,DS,DW,dXdxI,J,MW,Glob,KoeffDivW,ndrange=ndrangeB)
    KernelAbstractions.synchronize(backend)
  end  
  KGradKernel!(F,U,p,DS,dXdxI,J,X,M,MRho,Glob,GravitationFun,ndrange=ndrangeB)
  KernelAbstractions.synchronize(backend)
  KMomentumCoriolisKernel!(F,U,DS,dXdxI,J,X,MRho,M,Glob,CoriolisFun,ndrange=ndrangeB)
  KernelAbstractions.synchronize(backend)
  KRhoGradKinKernel!(F,U,DS,dXdxI,J,M,MRho,Glob,ndrange=ndrangeB)
  KernelAbstractions.synchronize(backend)
  KDivRhoThUpwind3Kernel!(F,U,DS,dXdxI,J,M,Glob,ndrange=ndrangeB)
  KernelAbstractions.synchronize(backend)
  for iT = 1 : NumTr
    @views KDivRhoTrUpwind3Kernel!(F[:,:,iT+NumV],U[:,:,iT+NumV],U,DS,
      dXdxI,J,M,Glob,ndrange=ndrangeB)
    KernelAbstractions.synchronize(backend)
  end  
  if Global.Model.VerticalDiffusionMom
    KVerticalDiffusionMomentumKernel! = VerticalDiffusionMomentumKernel!(backend,group)
    KV = Cache.KV
    KVerticalDiffusionMomentumKernel!(F,U,KV,dXdxI,J,M,Glob,ndrange=ndrangeB)
  end    
  if Global.Model.VerticalDiffusion
    KVerticalDiffusionScalarKernel! = VerticalDiffusionScalarKernel!(backend,groupTr)  
    KV = Cache.KV
    @views KVerticalDiffusionScalarKernel!(F[:,:,5],U[:,:,5],Rho,KV,
      dXdxI,J,M,Glob,ndrange=ndrangeB)
    KernelAbstractions.synchronize(backend)
    for iT = 1 : NumTr
      @views KVerticalDiffusionScalarKernel!(F[:,:,iT+NumV],U[:,:,iT+NumV],Rho,KV,
        dXdxI,J,M,Glob,ndrange=ndrangeB)
      KernelAbstractions.synchronize(backend)
    end  
  end    
  if Global.Model.SurfaceFluxMom
    NFG = min(div(NumberThreadGPU,N*N),NF)
    groupS = (N * N, NFG)
    KSurfaceFluxMomentumKernel! = SurfaceFluxMomentumKernel!(backend,groupS)  
    ndrangeSB = (N * N,NBF)
    CM = Global.SurfaceData.CM
    KSurfaceFluxMomentumKernel!(F,U,dXdxI,nS,CM,M,Glob,ndrange=ndrangeSB)
  end  
  if Global.Model.SurfaceFlux
    NFG = min(div(NumberThreadGPU,N*N),NF)
    groupS = (N * N, NFG)
    KSurfaceFluxScalarsKernel = SurfaceFluxScalarsKernel(backend,groupS)  
    ndrangeSB = (N * N,NBF)
    CT = Global.SurfaceData.CT
    CH = Global.SurfaceData.CH
    uStar = Global.SurfaceData.uStar
    TSurf = Global.SurfaceData.TS
    RhoVSurf = Global.SurfaceData.RhoVS
    KSurfaceFluxScalarsKernel(F,U,p,TSurf,RhoVSurf,uStar,CT,CH,dXdxI,Glob,M,Phys,ndrange=ndrangeSB)
    KernelAbstractions.synchronize(backend)
  end  

  Parallels.ExchangeData3DSendGPU(F,Exchange)


  KHyperViscKoeffKernel!(F,U,CacheF,DS,DW,dXdxI_I,J_I,M,Glob_I,KoeffCurl,KoeffGrad,KoeffDiv,ndrange=ndrangeI)
  KernelAbstractions.synchronize(backend)
  for iT = 1 : NumTr
    @views CacheTr = Temp1[:,:,iT + 6]  
    @views KHyperViscTracerKoeffKernel!(F[:,:,iT+NumV],CacheTr,Rho,DS,DW,dXdxI_I,J_I,M,Glob_I,
      KoeffDiv,ndrange=ndrangeI)
    KernelAbstractions.synchronize(backend)
  end  
  if KoeffDivW > 0
    @views KHyperViscWKoeffKernel!(F[:,:,4],Cachew,DS,DW,dXdxI_I,J_I,MW,Glob,KoeffDivW,ndrange=ndrangeI)
    KernelAbstractions.synchronize(backend)
  end  
  KGradKernel!(F,U,p,DS,dXdxI_I,J_I,X_I,M,MRho,Glob_I,GravitationFun,ndrange=ndrangeI)
  KernelAbstractions.synchronize(backend)
  KMomentumCoriolisKernel!(F,U,DS,dXdxI_I,J_I,X_I,MRho,M,Glob_I,CoriolisFun,ndrange=ndrangeI)
  KernelAbstractions.synchronize(backend)
  KRhoGradKinKernel!(F,U,DS,dXdxI_I,J_I,M,MRho,Glob_I,ndrange=ndrangeI)
  KernelAbstractions.synchronize(backend)
  KDivRhoThUpwind3Kernel!(F,U,DS,dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeI)
  KernelAbstractions.synchronize(backend)
  for iT = 1 : NumTr
    @views KDivRhoTrUpwind3Kernel!(F[:,:,iT+NumV],U[:,:,iT+NumV],U,DS,
      dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeI)
    KernelAbstractions.synchronize(backend)
  end  
  if Global.Model.VerticalDiffusionMom
    KV = Cache.KV
    KVerticalDiffusionMomentumKernel!(F,U,KV,dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeI)
  end    
  if Global.Model.VerticalDiffusion
    @views KV_I = Cache.KV[:,:,NBF+1:NF]
    @views KVerticalDiffusionScalarKernel!(F[:,:,5],U[:,:,5],Rho,KV_I,
      dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeI)
    KernelAbstractions.synchronize(backend)
    for iT = 1 : NumTr
      @views KVerticalDiffusionScalarKernel!(F[:,:,iT+NumV],U[:,:,iT+NumV],Rho,KV_I,
        dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeI)
      KernelAbstractions.synchronize(backend)
    end  
  end    
  if Global.Model.SurfaceFluxMom
    ndrangeSI = (N * N,NF-NBF)  
    @views CM_I = Global.SurfaceData.CM[:,NBF+1:NF]
    KSurfaceFluxMomentumKernel!(F,U,dXdxI_I,nS_I,CM_I,M,Glob_I,ndrange=ndrangeSI)
  end  
  if Global.Model.SurfaceFlux
    ndrangeSI = (N * N,NF-NBF)  
    @views CT_I = Global.SurfaceData.CT[:,NBF+1:NF]
    @views CH_I = Global.SurfaceData.CH[:,NBF+1:NF]
    @views uStar_I = Global.SurfaceData.uStar[:,NBF+1:NF]
    @views TSurf_I = Global.SurfaceData.TS[:,NBF+1:NF]
    @views RhoVSurf_I = Global.SurfaceData.RhoVS[:,NBF+1:NF]
    KSurfaceFluxScalarsKernel(F,U,p,TSurf_I,RhoVSurf_I,uStar_I,CT_I,CH_I,dXdxI_I,Glob_I,M,Phys,ndrange=ndrangeSI)
    KernelAbstractions.synchronize(backend)
  end  
      
  Parallels.ExchangeData3DRecvGPU!(F,Exchange)
  KernelAbstractions.synchronize(backend)

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

  if Global.Model.Damping
    KDampKernel! = DampKernel!(backend, groupG)
    KDampKernel!(Damp,F,U,zP,ndrange=ndrangeG)
    KernelAbstractions.synchronize(backend)
  end  

end

function FcnGPU!(F,U,FE,Metric,Phys,Cache,Exchange,Global,Param,Equation::Models.CompressibleDeep)

  backend = get_backend(F)
  FT = eltype(F)
  Glob = FE.Glob
  DS = FE.DS
  DW = FE.DW
  M = FE.M
  MW = FE.MW
  dXdx = Metric.dXdx
  dXdxI = Metric.dXdxI
  X = Metric.X
  J = Metric.J
  NF = Global.Grid.NumFaces
  NBF = Global.Grid.NumBoundaryFaces
  @views dXdxI_B = dXdxI[:,:,:,:,:,1:NBF]
  @views J_B = J[:,:,:,1:NBF]
  @views X_B = X[:,:,:,:,1:NBF]
  @views Glob_B = Glob[:,1:NBF]

  @views dXdx_I = dXdx[:,:,:,:,:,NBF+1:NF]
  @views dXdxI_I = dXdxI[:,:,:,:,:,NBF+1:NF]
  @views J_I = J[:,:,:,NBF+1:NF]
  @views X_I = X[:,:,:,:,NBF+1:NF]
  @views Glob_I = Glob[:,NBF+1:NF]
  lat = Metric.lat  
  dz = Metric.dz  
  zP = Metric.zP  
  DoF = FE.DoF
  N = size(FE.DS,1)
  Nz = size(F,1)
  NDoF = size(F,2)
  NumV  = Global.Model.NumV 
  NumTr  = Global.Model.NumTr 
  Koeff = Global.Model.HyperDDiv
  Temp1 = Cache.Temp1
  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU
  Force = Global.Model.Force
  Damp = Global.Model.Damp
  MicrophysicsSource = Global.Model.MicrophysicsSource

  KoeffCurl = Global.Model.HyperDCurl
  KoeffGrad = Global.Model.HyperDGrad
  KoeffDiv = Global.Model.HyperDDiv
  KoeffDivW = Global.Model.HyperDDivW

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
  @views CacheFF = Temp1[:,:,1:6+NumTr+1]
  @views Cachew = Temp1[:,:,6 + 1 + NumTr]  
  @views p = Cache.AuxG[:,:,1]
  KV = Cache.KV
  TSurf = Cache.TSurf
  RhoVSurf = Cache.RhoVSurf
  uStar = Cache.uStar
  CT = Cache.CT
  CH = Cache.CH
  @views KV_I = Cache.KV[:,:,NBF+1:NF]
  @views TSurf_I = Cache.TSurf[:,NBF+1:NF]
  @views RhoVSurf_I = Cache.RhoVSurf[:,NBF+1:NF]
  @views uStar_I = Cache.uStar[:,NBF+1:NF]
  @views CH_I = Cache.CH[:,NBF+1:NF]
# Ranges
  NzG = min(div(NumberThreadGPU,N*N),Nz)
  group = (N, N, NzG, 1)
  ndrange = (N, N, Nz, NF)
  ndrangeB = (N, N, Nz, NBF)
  ndrangeI = (N, N, Nz, NF-NBF)
  groupTr = group
  ndrangeTr = ndrange
  NDoFG = min(div(NumberThreadGPU,Nz),NDoF)
  groupG = (Nz, NDoFG)  
  ndrangeG = (Nz, NDoF)  
  NzG = min(div(NumberThreadGPU,N*N),Nz-1)
  groupw = (N, N, NzG, 1)
  ndrangewB = (Nz-1, NBF)  
  ndrangewI = (Nz-1, NF-NBF)  

  KGradDeepKernel! = GradDeepKernel!(backend,group)
  KHyperViscKernel! = HyperViscKernel!(backend, group)
  KHyperViscKoeffKernel! = HyperViscKoeffKernel!(backend, group)
  KDivRhoThUpwind3Kernel! = DivRhoThUpwind3Kernel!(backend, group)
  KMomentumDeepCoriolisKernel! = MomentumDeepCoriolisKernel!(backend, group)
  KHyperViscTracerKernel! = HyperViscTracerKernel!(backend, groupTr)
  KHyperViscTracerKoeffKernel! = HyperViscTracerKoeffKernel!(backend, groupTr)
  KHyperViscWKernel! = HyperViscWKernel!(backend, groupTr)
  KHyperViscWKoeffKernel! = HyperViscWKoeffKernel!(backend, groupTr)
  KDivRhoTrUpwind3Kernel! = DivRhoTrUpwind3Kernel!(backend, groupTr)
  if Global.Model.SurfaceFlux
    NFG = min(div(NumberThreadGPU,N*N),NF)
    groupS = (N * N, NFG)
    KSurfaceFluxScalarsKernel = SurfaceFluxScalarsKernel(backend,groupS)  
  end  
  if Global.Model.VerticalDiffusion
    KVerticalDiffusionScalarKernel! = VerticalDiffusionScalarKernel!(backend,groupTr)  
  end  


####
# First phase  
####
  Temp1 .= FT(0)
  @views MRho = CacheF[:,:,6]
  KHyperViscKernel!(CacheF,MRho,U,DS,DW,dXdxI,J,M,Glob,ndrange=ndrangeB)
  KernelAbstractions.synchronize(backend)
  for iT = 1 : NumTr
    @views CacheTr = Temp1[:,:,iT + 6]  
    @views KHyperViscTracerKernel!(CacheTr,U[:,:,iT+NumV],Rho,DS,DW,dXdxI,J,M,Glob,ndrange=ndrangeB)
    KernelAbstractions.synchronize(backend)
  end  
  @views KHyperViscWKernel!(Cachew,U[:,:,4],DS,DW,dXdxI,J,MW,Glob,ndrange=ndrangeB)
  KernelAbstractions.synchronize(backend)
  Parallels.ExchangeData3DSendGPU(CacheFF,Exchange)

  KHyperViscKernel!(CacheF,MRho,U,DS,DW,dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeI)
  KernelAbstractions.synchronize(backend)
  for iT = 1 : NumTr
    @views CacheTr = Temp1[:,:,iT + 6]  
    @views KHyperViscTracerKernel!(CacheTr,U[:,:,iT+NumV],Rho,DS,DW,dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeI)
    KernelAbstractions.synchronize(backend)
  end  
  @views KHyperViscWKernel!(Cachew,U[:,:,4],DS,DW,dXdxI_I,J_I,MW,Glob_I,ndrange=ndrangeI)
  KernelAbstractions.synchronize(backend)

  Parallels.ExchangeData3DRecvGPU!(CacheFF,Exchange)
  KernelAbstractions.synchronize(backend)

####
# Second phase  
####

  F .= FT(0)
  KHyperViscKoeffKernel!(F,U,CacheF,DS,DW,dXdxI,J,M,Glob,KoeffCurl,KoeffGrad,KoeffDiv,ndrange=ndrangeB)
  KernelAbstractions.synchronize(backend)
  for iT = 1 : NumTr
    @views CacheTr = Temp1[:,:,iT + 6]  
    @views KHyperViscTracerKoeffKernel!(F[:,:,iT+NumV],CacheTr,Rho,DS,DW,dXdxI,J,M,Glob,
      KoeffDiv,ndrange=ndrangeB)
    KernelAbstractions.synchronize(backend)
  end  
  @views KHyperViscWKoeffKernel!(F[:,:,4],Cachew,DS,DW,dXdxI,J,MW,Glob,KoeffDivW,ndrange=ndrangeB)
  KernelAbstractions.synchronize(backend)
  KGradDeepKernel!(F,U,p,DS,dXdxI,J,X,M,MRho,Glob,Phys,ndrange=ndrangeB)
  KernelAbstractions.synchronize(backend)
  KMomentumDeepCoriolisKernel!(F,U,DS,dXdxI,J,X,MRho,M,Glob,Phys,ndrange=ndrangeB)
  KernelAbstractions.synchronize(backend)
  KRhoGradKinKernel!(F,U,DS,dXdxI,J,M,MRho,Glob)
  KernelAbstractions.synchronize(backend)
  KDivRhoThUpwind3Kernel!(F,U,DS,dXdxI,J,M,Glob,ndrange=ndrangeB)
  KernelAbstractions.synchronize(backend)
  for iT = 1 : NumTr
    @views KDivRhoTrUpwind3Kernel!(F[:,:,iT+NumV],U[:,:,iT+NumV],U,DS,
      dXdxI,J,M,Glob,ndrange=ndrangeB)
    KernelAbstractions.synchronize(backend)
  end  
  if Global.Model.VerticalDiffusion
    @views KVerticalDiffusionScalarKernel!(F[:,:,5],U[:,:,5],Rho,KV,
      dXdxI,J,M,Glob,ndrange=ndrangeB)
    KernelAbstractions.synchronize(backend)
    for iT = 1 : NumTr
      @views KVerticalDiffusionScalarKernel!(F[:,:,iT+NumV],U[:,:,iT+NumV],Rho,KV,
        dXdxI,J,M,Glob,ndrange=ndrangeB)
      KernelAbstractions.synchronize(backend)
    end  
  end    
  if Global.Model.SurfaceFlux
    ndrangeSB = (N * N,NBF)
    KSurfaceFluxScalarsKernel(F,U,p,TSurf,RhoVSurf,uStar,CT,CH,dXdxI,Glob,M,Phys,ndrange=ndrangeSB)
    KernelAbstractions.synchronize(backend)
  end  

  Parallels.ExchangeData3DSendGPU(F,Exchange)


  KHyperViscKoeffKernel!(F,U,CacheF,DS,DW,dXdxI_I,J_I,M,Glob_I,KoeffCurl,KoeffGrad,KoeffDiv,ndrange=ndrangeI)
  KernelAbstractions.synchronize(backend)
  for iT = 1 : NumTr
    @views CacheTr = Temp1[:,:,iT + 6]  
    @views KHyperViscTracerKoeffKernel!(F[:,:,iT+NumV],CacheTr,Rho,DS,DW,dXdxI_I,J_I,M,Glob_I,
      KoeffDiv,ndrange=ndrangeI)
    KernelAbstractions.synchronize(backend)
  end  
  @views KHyperViscWKoeffKernel!(F[:,:,4],Cachew,DS,DW,dXdxI_I,J_I,MW,Glob,KoeffDivW,ndrange=ndrangeI)
  KernelAbstractions.synchronize(backend)
  KGradDeepKernel!(F,U,p,DS,dXdxI_I,J_I,X_I,M,MRho,Glob_I,Phys,ndrange=ndrangeI)
  KernelAbstractions.synchronize(backend)
  KMomentumDeepCoriolisKernel!(F,U,DS,dXdxI_I,J_I,X_I,MRho,M,Glob_I,Phys,ndrange=ndrangeI)
  KernelAbstractions.synchronize(backend)
  KRhoGradKinKernel!(F,U,DS,dXdxI_I,J_I,M,MRho,Glob)
  KernelAbstractions.synchronize(backend)

  KDivRhoThUpwind3Kernel!(F,U,DS,dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeI)
  KernelAbstractions.synchronize(backend)
  for iT = 1 : NumTr
    @views KDivRhoTrUpwind3Kernel!(F[:,:,iT+NumV],U[:,:,iT+NumV],U,DS,
      dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeI)
    KernelAbstractions.synchronize(backend)
  end  
  if Global.Model.VerticalDiffusion
    @views KVerticalDiffusionScalarKernel!(F[:,:,5],U[:,:,5],Rho,KV_I,
      dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeI)
    KernelAbstractions.synchronize(backend)
    for iT = 1 : NumTr
      @views KVerticalDiffusionScalarKernel!(F[:,:,iT+NumV],U[:,:,iT+NumV],Rho,KV_I,
        dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeI)
      KernelAbstractions.synchronize(backend)
    end  
  end    
  if Global.Model.SurfaceFlux
    ndrangeSB = (N * N,NBF)  
    ndrangeSI = (N * N,NF-NBF)  
    KSurfaceFluxScalarsKernel(F,U,p,TSurf_I,RhoVSurf_I,uStar_I,CT_I,CH_I,dXdxI_I,Glob_I,M,Phys,ndrange=ndrangeSI)
    KernelAbstractions.synchronize(backend)
  end  
      
  Parallels.ExchangeData3DRecvGPU!(F,Exchange)
  KernelAbstractions.synchronize(backend)

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

  if Global.Model.Damping
    KDampKernel! = DampKernel!(backend, groupG)
    KDampKernel!(Damp,F,U,zP,ndrange=ndrangeG)
    KernelAbstractions.synchronize(backend)
  end  

end

function FcnGPUAMD!(F,U,FE,Metric,Phys,Cache,Exchange,Global,Param,Equation::Models.CompressibleShallow)

  backend = get_backend(F)
  FT = eltype(F)
  Glob = FE.Glob
  DS = FE.DS
  DW = FE.DW
  M = FE.M
  MW = FE.MW
  dXdxI = Metric.dXdxI
  X = Metric.X
  J = Metric.J
  NF = Global.Grid.NumFaces
  NBF = Global.Grid.NumBoundaryFaces
  @views dXdxI_I = dXdxI[:,:,:,:,:,NBF+1:NF]
  @views J_I = J[:,:,:,NBF+1:NF]
  @views X_I = X[:,:,:,:,NBF+1:NF]
  @views Glob_I = Glob[:,NBF+1:NF]
  lat = Metric.lat  
  dz = Metric.dz  
  DoF = FE.DoF
  N = size(FE.DS,1)
  Nz = size(F,1)
  NDoF = size(F,2)
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
  @views MRho = CacheF[:,:,6]
  copyto!(MRho,FT(1))
# Ranges
  NzG = min(div(NumberThreadGPU,N*N),Nz)
  group = (N, N, NzG, 1)
  ndrange = (N, N, Nz, NF)
  ndrangeB = (N, N, Nz, NBF)
  ndrangeI = (N, N, Nz, NF-NBF)
  groupTr = group
  ndrangeTr = ndrange
  NDoFG = min(div(NumberThreadGPU,Nz),NDoF)
  groupG = (Nz, NDoFG)  
  ndrangeG = (Nz, NDoF)  
  CoriolisFun = Global.Model.CoriolisFun

  KMomentumCoriolisKernel! = MomentumVectorInvariantCoriolisKernel!(backend, group)

####
# Second phase  
####

  @show "FcnGPUAMD"
  for i = 1 : 100
    F .= FT(0)
    KMomentumCoriolisKernel!(F,U,DS,dXdxI,J,X,MRho,M,Glob,CoriolisFun,ndrange=ndrangeB)
    KMomentumCoriolisKernel!(F,U,DS,dXdxI_I,J_I,X_I,MRho,M,Glob_I,CoriolisFun,ndrange=ndrangeI)
  end

end


