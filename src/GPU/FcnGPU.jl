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
  N = FE.OrdPoly+1
  ww = FE.w
  Nz = size(F,1)
  NF = size(Glob,2)
  NDoF = size(Glob,1)
  KoeffDiv = Global.Model.HyperDDiv
  NumV  = Global.Model.NumV
  NumTr  = Global.Model.NumTr
  Temp1 = Cache.Temp1
  @views q = Cache.q[:,:,1:2*NumTr]

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
  NF = Global.Grid.NumFaces
  NBF = Global.Grid.NumFacesB
  ndrangeB = (N, N, Nz, NBF)
  ndrangeI = (N, N, Nz, NF-NBF)

  @views dXdxI_I = dXdxI[:,:,:,:,:,NBF+1:NF]
  @views J_I = J[:,:,:,NBF+1:NF]
  @views X_I = X[:,:,:,:,NBF+1:NF]
  @views Glob_I = Glob[:,NBF+1:NF]
  @views Stencil_I = Stencil[NBF+1:NF,:]

  KLimitKernel! = LimitKernel!(backend, groupL)
  KuvwFunCKernel! = uvwFunCKernel!(backend, groupC)
  KDivRhoKernel! = DivRhoKernel!(backend, group)
  KHyperViscTracerKernel! = HyperViscTracerKernel!(backend, group)
  KHyperViscTracerKoeffKernel! = HyperViscTracerKoeffKernel!(backend, group)
# KDivRhoTrUpwind3Kernel! = DivRhoTrUpwind3Kernel!(backend, group)
  KAdvectionTrUpwind3Kernel! = AdvectionTrUpwind3Kernel!(backend, group)
  KDivRhoTrUpwind3LimKernel! = DivRhoTrUpwind3LimKernel!(backend, group)
  KDivRhoTrViscUpwind3LimKernel! = DivRhoTrViscUpwind3LimKernel!(backend, group)

  if Global.Model.HorLimit
    @views KLimitKernel!(DoF,q[:,:,nT],q[:,:,nT+1:2*nT],U[:,:,NumV+1:NumV+NumTr],Rho,Glob,ndrange=ndrangeL)
    KernelAbstractions.synchronize(backend)
    Parallels.ExchangeDataFSendGPU(q,Exchange)
    Parallels.ExchangeDataFRecvGPU!(q,Exchange)  
  end


# Velocity 
  KuvwFunCKernel!(Profile,u,v,w,time,Glob,X,Param,Phys,ndrange=ndrangeC)
  KernelAbstractions.synchronize(backend)

# Hyperviscosity Part 1
  if ~Global.Model.HorLimit
    @. CacheTr = FT(0)  
    KHyperViscTracerKernel!(CacheTr,U[:,:,1+NumV],Rho,DS,DW,dXdxI,J,M,Glob,ndrange=ndrangeB)
    KernelAbstractions.synchronize(backend)

#   Data exchange  
    Parallels.ExchangeData3DSendGPU(CacheTr,Exchange)

    KHyperViscTracerKernel!(CacheTr,U[:,:,1+NumV],Rho,DS,DW,dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeI)
    KernelAbstractions.synchronize(backend)
    Parallels.ExchangeData3DRecvGPU!(CacheTr,Exchange)
  end  

  F .= FT(0)

  KDivRhoKernel!(F,U,DS,dXdxI,J,M,Glob,ndrange=ndrangeB)
  KernelAbstractions.synchronize(backend)  

  if Global.Model.HorLimit
    @views KDivRhoTrUpwind3LimKernel!(F[:,:,1+NumV],U[:,:,1+NumV],U,DS,
      dXdxI,J,M,Glob,dtau,ww,q[:,:,1],q[:,:,nT+1],Stencil,ndrange=ndrangeB)
    KernelAbstractions.synchronize(backend)  
  else
    @views KHyperViscTracerKoeffKernel!(F[:,:,1+NumV],CacheTr,Rho,DS,DW,dXdxI,J,M,Glob,
      KoeffDiv,ndrange=ndrangeB)
    KernelAbstractions.synchronize(backend)

#   @views KDivRhoTrUpwind3Kernel!(F[:,:,1+NumV],U[:,:,1+NumV],U,DS,
#     dXdxI,J,M,Glob,ndrange=ndrangeB)
    @views KAdvectionTrUpwind3Kernel!(F[:,:,1+NumV],U[:,:,1+NumV],U,U[:,:,4],DS,dXdxI,
      J,M,Glob,ndrange=ndrangeB)
    KernelAbstractions.synchronize(backend)  
  end

# Data exchange  
  @views Parallels.ExchangeData3DSendGPU(F[:,:,1:1+NumV],Exchange)

  KDivRhoKernel!(F,U,DS,dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeI)
  KernelAbstractions.synchronize(backend)

  if Global.Model.HorLimit
    @views KDivRhoTrUpwind3LimKernel!(F[:,:,1+NumV],U[:,:,1+NumV],U,DS,
      dXdxI_I,J_I,M,Glob_I,dtau,ww,q[:,:,1],q[:,:,nT+1],Stencil_I,ndrange=ndrangeI)
    KernelAbstractions.synchronize(backend)
  else
    @views KHyperViscTracerKoeffKernel!(F[:,:,1+NumV],CacheTr,Rho,DS,DW,dXdxI_I,J_I,M,Glob_I,
      KoeffDiv,ndrange=ndrangeI)
    KernelAbstractions.synchronize(backend)

#   @views KDivRhoTrUpwind3Kernel!(F[:,:,1+NumV],U[:,:,1+NumV],U,DS,
#     dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeI)
    @views KAdvectionTrUpwind3Kernel!(F[:,:,1+NumV],U[:,:,1+NumV],U,U[:,:,4],DS,
      dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeI)
    KernelAbstractions.synchronize(backend)
  end

# Data exchange  
  @views Parallels.ExchangeData3DRecvGPU!(F[:,:,1:1+NumV],Exchange)
end

function FcnGPU!(F,U,FE,Metric,Phys,Cache,Exchange,Global,Param,Equation::Models.CompressibleShallow)

  backend = get_backend(F)
  FT = eltype(F)
  State = Global.Model.State
  dtau = Global.TimeStepper.dtauStage
  Glob = FE.Glob
  DS = FE.DS
  DW = FE.DW
  M = FE.M
  Stencil = FE.Stencil
  BoundaryDoF  = FE.BoundaryDoF 
  dXdxI = Metric.dXdxI
  nS = Metric.nS
  nSS = Metric.nSS
  X = Metric.X
  J = Metric.J
  N = FE.OrdPoly+1
  ww = FE.w
  NF = Global.Grid.NumFaces
  NBF = Global.Grid.NumFacesB
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
  @views Stencil_I = Stencil[NBF+1:NF,:]
  xS = Metric.xS  
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
  Proc = Global.ParallelCom.Proc
  Force = Global.Model.Force
  Damp = Global.Model.Damp
  EDMF = Global.Model.EDMF
  ND = Global.Model.NDEDMF
  MicrophysicsSource = Global.Model.MicrophysicsSource
  CoriolisFun = Global.Model.CoriolisFun
  GravitationFun = Global.Model.GravitationFun
  HorLimit = Global.Model.HorLimit

  KoeffCurl = Global.Model.HyperDCurl
  KoeffGrad = Global.Model.HyperDGrad
  KoeffDiv = Global.Model.HyperDDiv
  KoeffDivW = Global.Model.HyperDDivW

# Position  
  RhoPos = Global.Model.RhoPos
  uPos = Global.Model.uPos
  vPos = Global.Model.vPos
  wPos = Global.Model.wPos
  ThPos = Global.Model.ThPos
  TkePos = Global.Model.TkePos
# State vector
  @views Rho = U[:,:,RhoPos]
  @views u = U[:,:,uPos]
  @views v = U[:,:,vPos]
  @views w = U[:,:,wPos]
  @views Th = U[:,:,ThPos]
  TrPos = NumV
  if TkePos > 0
    @views Tke = U[:,:,TkePos]  
    TrPos += 1
  end  
  @views UTr = U[:,:,TrPos+1:TrPos+NumTr]
  if EDMF
    aRhoEDMFPos = NumV + NumTr  
    if TkePos > 0
      aRhoEDMFPos += 1
    end  
    wEDMFPos = aRhoEDMFPos + ND
    ThEDMFPos = wEDMFPos + ND
    TrEDMFPos = ThEDMFPos + ND
    @views aRhoEDMF = U[:,:,aRhoEDMFPos+ND-1] 
    @views wEDMF = U[:,:,wEDMFPos:wEDMFPos+ND-1]  
    @views ThEDMF = U[:,:,ThEDMFPos:ThEDMFPos+ND-1]  
    @views TrEDMF = U[:,:,TrEDMFPos:end]  
    @views FaRhoEDMF = F[:,:,aRhoEDMFPos+ND-1] 
    @views FwEDMF = F[:,:,wEDMFPos:wEDMFPos+ND-1]  
    @views FThEDMF = F[:,:,ThEDMFPos:ThEDMFPos+ND-1]  
    @views FTrEDMF = F[:,:,TrEDMFPos:end]  
    RhoEDMF = Cache.RhoEDMF
  end  
# Tendency
  @views FRho = F[:,:,1]
  if TkePos > 0
    @views FTke = F[:,:,TkePos]  
  end  
  @views FTr = F[:,:,TrPos+1:TrPos+NumTr]
# Cache
# Need clearer cache distribution for different setups
#   1...4 Horizontal momentum   
#       5 Thermodynamic variable
#      +1 Vertical momentum   
#      +1 Turbulent kinetic energy
#  +NumTr Tracer
#  +ND*(
#       +1 Thermodynamic variable
#       +1 Vertical velocity
#       +NumTr Tracer)
  LenTemp1 = 5
  if KoeffDivW > 0
    LenTemp1 += 1  
    @views Cachew = Temp1[:,:,LenTemp1]
  end  
  if TkePos > 0
    LenTemp1 += 1
    @views CacheTke = Temp1[:,:,LenTemp1]
  end  
  if ~HorLimit && NumTr > 0
    @views CacheTr = Temp1[:,:,LenTemp1+1:LenTemp1+NumTr]
    LenTemp1 += NumTr
  end  
  if EDMF
    @views CachewEDMF = Temp1[:,:,LenTemp1+1:LenTemp1+ND]
    LenTemp1 += ND 
    @views CacheThEDMF = Temp1[:,:,LenTemp1+1:LenTemp1+ND]
    LenTemp1 += ND 
    @views CacheTrEDMF = Temp1[:,:,LenTemp1+1:LenTemp1+ND*NumTr]
    LenTemp1 += ND * NumTr
  end  
      
  @views CacheF = Temp1[:,:,1:5]
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
  groupS = max(div(NDoF,NumberThreadGPU),1)
  ndrangeS = (NDoF)  
  NzG = min(div(NumberThreadGPU,N*N),Nz-1)
  groupw = (N, N, NzG, 1)
  ndrangewB = (Nz-1, NBF)  
  ndrangewI = (Nz-1, NF-NBF)  
  NFG = min(div(NumberThreadGPU,Nz),NF)
  groupL = (Nz, NFG, 1)
  ndrangeL = (Nz, NF, NumTr)

  KRhoGradKinKernel! = RhoGradKinKernel!(backend,group)
  KGradFullKernel! = GradFullKernel!(backend,group)
  KGradKernel! = GradKernel!(backend,group)
  KHyperViscKernel! = HyperViscKernel!(backend, group)
  KHyperViscKoeffKernel! = HyperViscKoeffKernel!(backend, group)
  if State == "Dry" || State == "ShallowWater" || State == "Moist"
    KDivRhoThUpwind3Kernel! = DivRhoThUpwind3Kernel!(backend, group)
  elseif State == "DryEnergy" || State == "MoistEnergy"
    KDivRhoKEUpwind3Kernel! = DivRhoKEUpwind3Kernel!(backend, group)
  end
  KDivRhoThUpwind3Kernel! = DivRhoThUpwind3Kernel!(backend, group)
  KMomentumCoriolisKernel! = MomentumVectorInvariantCoriolisKernel!(backend, group)
  KHyperViscTracerKernel! = HyperViscTracerKernel!(backend, groupTr)
  KHyperViscTracerKoeffKernel! = HyperViscTracerKoeffKernel!(backend, groupTr)
  KDivRhoTrUpwind3Kernel! = DivRhoTrUpwind3Kernel!(backend, groupTr)
  KDivRhoTrUpwind3New2Kernel! = DivRhoTrUpwind3New2Kernel!(backend, groupTr)
  KDivRhoTrUpwind3LimKernel! = DivRhoTrUpwind3LimKernel!(backend, groupTr)
  KLimitKernel! = LimitKernel!(backend, groupL)

# BoundaryValues
  @. @views U[:,BoundaryDoF,vPos] = FT(0.0)

  if HorLimit
    @views KLimitKernel!(DoF,q,UTr,Rho,Glob,ndrange=ndrangeL)
    Parallels.ExchangeDataFSendGPU(q,Exchange)
  end


####
# First phase  
####
  Temp1 .= FT(0)
  KHyperViscKernel!(CacheF,U,DS,DW,dXdxI,J,M,Glob,ndrange=ndrangeB)
  if ~HorLimit
    for iT = 1 : NumTr
      @views KHyperViscTracerKernel!(CacheTr[:,:,iT],UTr[:,:,iT],Rho,DS,DW,dXdxI,J,M,Glob,ndrange=ndrangeB)
    end  
  end  
  if TkePos > 0
    @views KHyperViscTracerKernel!(CacheTke,Tke,Rho,DS,DW,dXdxI,J,M,Glob,ndrange=ndrangeB)
  end  
      
  if KoeffDivW > 0
    KHyperViscWKernel! = HyperViscWKernel!(backend, groupTr)
    @views KHyperViscWKernel!(Cachew,U[:,:,4],DS,DW,dXdxI,J,M,Glob,ndrange=ndrangeB)
  end  
  
  if EDMF
    ndrangeBEDMF = (N, N, Nz, NBF, ND)
    KHyperViscWEDMFKernel! = HyperViscWEDMFKernel!(backend, groupTr)
    @views KHyperViscWEDMFKernel!(CachewEDMF,wEDMF,aRhoEDMF,DS,DW,dXdxI,J,M,Glob,ndrange=ndrangeBEDMF)
    KHyperViscTracerEDMFKernel! = HyperViscTracerEDMFKernel!(backend, groupTr)
    @views KHyperViscTracerEDMFKernel!(CacheThEDMF,ThEDMF,DS,DW,dXdxI,J,M,Glob,ndrange=ndrangeBEDMF)
    if NumTr > 0
      ndrangeBEDMFTr = (N, N, Nz, NBF, ND*NumTr)
      @views KHyperViscTracerEDMFKernel!(CacheTrEDMF,TrEDMF,DS,DW,dXdxI,J,M,Glob,ndrange=ndrangeBEDMFTr)
    end
  end    

  if HorLimit
    Parallels.ExchangeDataFRecvGPU!(q,Exchange)  
  end  
  @views Parallels.ExchangeData3DSendGPU(Temp1[:,:,1:LenTemp1],Exchange)

  KHyperViscKernel!(CacheF,U,DS,DW,dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeI)
  if ~HorLimit
    for iT = 1 : NumTr
      @views KHyperViscTracerKernel!(CacheTr[:,:,iT],UTr[:,:,iT],Rho,DS,DW,dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeI)
    end  
  end  
  if TkePos > 0
    @views KHyperViscTracerKernel!(CacheTke,Tke,Rho,DS,DW,dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeI)
  end  
  if KoeffDivW > 0
    @views KHyperViscWKernel!(Cachew,U[:,:,4],DS,DW,dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeI)
  end  
  if EDMF
    ndrangeIEDMF = (N, N, Nz, NF-NBF, ND)
    @views KHyperViscWEDMFKernel!(CachewEDMF,wEDMF,aRhoEDMF,DS,DW,dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeIEDMF)
    @views KHyperViscTracerEDMFKernel!(CacheThEDMF,ThEDMF,DS,DW,dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeIEDMF)
    if NumTr > 0
      ndrangeIEDMFTr = (N, N, Nz, NF-NBF, ND*NumTr)
      @views KHyperViscTracerEDMFKernel!(CacheTrEDMF,TrEDMF,DS,DW,dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeIEDMFTr)
    end
  end    

  @views Parallels.ExchangeData3DRecvGPU!(Temp1[:,:,1:LenTemp1],Exchange)

####
# Second phase  
####
  F .= FT(0)
  KHyperViscKoeffKernel!(F,U,CacheF,DS,DW,dXdxI,J,M,Glob,KoeffCurl,KoeffGrad,KoeffDiv,ndrange=ndrangeB)
  if ~HorLimit
    for iT = 1 : NumTr
      @views KHyperViscTracerKoeffKernel!(FTr[:,:,iT],CacheTr[:,:,iT],Rho,DS,DW,dXdxI,J,M,Glob,
        KoeffDiv,ndrange=ndrangeB)
#     @views KDivRhoTrUpwind3Kernel!(FTr[:,:,iT],UTr[:,:,iT],U,DS,
#       dXdxI,J,M,Glob,ndrange=ndrangeB)
    end  
  else  
    for iT = 1 : NumTr
      @views KDivRhoTrUpwind3LimKernel!(FTr[:,:,iT],UTr[:,:,iT],U,DS,
         dXdxI,J,M,Glob,dtau,ww,q[:,:,iT],q[:,:,NumTr+iT],Stencil,ndrange=ndrangeB)
    end  
  end  
  if TkePos > 0
    @views KHyperViscTracerKoeffKernel!(FTke,CacheTke,Rho,DS,DW,dXdxI,J,M,Glob,
      KoeffDiv,ndrange=ndrangeB)
    @views KDivRhoTrUpwind3Kernel!(FTke,Tke,U,DS, dXdxI,J,M,Glob,ndrange=ndrangeB)
  end  
  if KoeffDivW > 0
    KHyperViscWKoeffKernel! = HyperViscWKoeffKernel!(backend, groupTr)
    @views KHyperViscWKoeffKernel!(F[:,:,4],Cachew,DS,DW,dXdxI,J,M,Glob,KoeffDivW,ndrange=ndrangeB)
  end  
  if EDMF
    KHyperViscWKoeffEDMFKernel! = HyperViscWKoeffEDMFKernel!(backend, groupTr)
    @views KHyperViscWKoeffEDMFKernel!(FwEDMF,CachewEDMF,DS,DW,dXdxI,J,M,Glob,KoeffDivW,ndrange=ndrangeBEDMF)
    @views KHyperViscTracerKoeffEDMFKernel!(FThEDMF,CacheThEDMF,DS,DW,dXdxI,J,M,Glob,KoeffDiv,ndrange=ndrangeBEDMF)
      @views KHyperViscTracerKoeffEDMFKernel!(FTrEDMF,CacheTrEDMF,DS,DW,dXdxI,J,M,Glob,KoeffDiv,ndrange=ndrangeBEDMFTr)
  end  
  KMomentumCoriolisKernel!(F,U,DS,dXdxI,J,X,M,Glob,CoriolisFun,ndrange=ndrangeB)
  KGradFullKernel!(F,U,p,DS,dXdxI,X,J,M,Glob,GravitationFun,ndrange=ndrangeB)
  if State == "Dry" || State == "ShallowWater" || State == "Moist"
#   KDivRhoThUpwind3Kernel!(F,U,DS,dXdxI,J,M,Glob,ndrange=ndrangeB)
    KDivRhoTrUpwind3New2Kernel!(F,NumV,NumTr,U,DS,dXdxI,J,M,Glob,ndrange=ndrangeB)
  elseif State == "DryEnergy" || State == "MoistEnergy"
    KDivRhoKEUpwind3Kernel!(F,U,p,DS,dXdxI,J,M,Glob,ndrange=ndrangeB)
  end
  if EDMF
    KMomentumCoriolisDraftKernel! = MomentumVectorInvariantCoriolisDraftKernel!(backend,group)  
    KMomentumCoriolisDraftKernel!(F,U,wEDMF,aRhoEDMF,DS,dXdxI,J,X,M,Glob,CoriolisFun,ndrange=ndrangeBEDMF)
    KRhoGradKinEDMFKernel! = RhoGradKinEDMFKernel!(backend,group)
    KRhoGradKinEDMFKernel!(F,U,wEDMF,aRhoEDMF,DS,dXdxI,J,M,Glob,ndrange=ndrangeBEDMF)
    KAdvectionTrUpwind3Kernel! = AdvectionTrUpwind3Kernel!(backend,groupTr)
    KAdvectionTrUpwind3Kernel!(FThEDMF,ThEDMF,U,wEDMF,DS,dXdxI,J,M,Glob,ndrange=ndrangeBEDMF)
    if NumTr > 0
      KAdvectionTrUpwind3Kernel!(FTrEDMF,TrEDMF,U,wEDMF,DS,dXdxI,J,M,Glob,ndrange=ndrangeBEDMFTr)
    end  
  end    

  Parallels.ExchangeData3DSendGPU(F,Exchange)

  KHyperViscKoeffKernel!(F,U,CacheF,DS,DW,dXdxI_I,J_I,M,Glob_I,KoeffCurl,KoeffGrad,KoeffDiv,ndrange=ndrangeI)
  if ~HorLimit
    for iT = 1 : NumTr
      @views KHyperViscTracerKoeffKernel!(FTr[:,:,iT],CacheTr[:,:,iT],Rho,DS,DW,dXdxI_I,J_I,M,Glob_I,
        KoeffDiv,ndrange=ndrangeI)
    end  
#   for iT = 1 : NumTr
#     @views KDivRhoTrUpwind3Kernel!(FTr[:,:,iT],UTr[:,:,iT],U,DS,
#       dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeI)
#   end  
  else  
    for iT = 1 : NumTr
      @views KDivRhoTrUpwind3LimKernel!(FTr[:,:,iT],UTr[:,:,iT],U,DS,
        dXdxI_I,J_I,M,Glob_I,dtau,ww,q[:,:,iT],q[:,:,NumTr+iT],Stencil_I,ndrange=ndrangeI)
    end  
  end  
  if TkePos > 0
    @views KHyperViscTracerKoeffKernel!(FTke,CacheTke,Rho,DS,DW,dXdxI_I,J_I,M,Glob_I,
      KoeffDiv,ndrange=ndrangeI)
    KDivRhoTrUpwind3Kernel!(FTke,Tke,U,DS,dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeI)
  end  
  if KoeffDivW > 0
    @views KHyperViscWKoeffKernel!(F[:,:,4],Cachew,DS,DW,dXdxI_I,J_I,M,Glob_I,KoeffDivW,ndrange=ndrangeI)
  end  
  if EDMF
    @views KHyperViscWKoeffEDMFKernel!(FwEDMF,CachewEDMF,DS,DW,dXdxI_I,J_I,M,Glob_I,KoeffDivW,ndrange=ndrangeIEDMF)
    @views KHyperViscTracerKoeffEDMFKernel!(FThEDMF,CacheThEDMF,DS,DW,dXdxI_I,J_I,M,Glob_I,
      KoeffDiv,ndrange=ndrangeIEDMF)
    if NumTr > 0
      @views KHyperViscTracerKoeffEDMFKernel!(FTrEDMF,CacheTrEDMF,DS,DW,dXdxI_I,J_I,M,Glob_I,
        KoeffDiv,ndrange=ndrangeIEDMFTr)
    end  
  end  

  KMomentumCoriolisKernel!(F,U,DS,dXdxI_I,J_I,X_I,M,Glob_I,CoriolisFun,ndrange=ndrangeI)
  KGradFullKernel!(F,U,p,DS,dXdxI_I,X_I,J_I,M,Glob_I,GravitationFun,ndrange=ndrangeI)

  if State == "Dry" || State == "ShallowWater" || State == "Moist"
#   KDivRhoThUpwind3Kernel!(F,U,DS,dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeI)
    KDivRhoTrUpwind3New2Kernel!(F,NumV,NumTr,U,DS,dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeI)
  elseif State == "DryEnergy"
    KDivRhoKEUpwind3Kernel!(F,U,p,DS,dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeI)
  end

  if EDMF
    KMomentumCoriolisDraftKernel!(F,U,wEDMF,aRhoEDMF,DS,dXdxI_I,J_I,X_I,M,Glob_I,CoriolisFun,ndrange=ndrangeIEDMF)
    KRhoGradKinEDMFKernel!(F,U,wEDMF,aRhoEDMF,DS,dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeIEDMF)
    KAdvectionTrUpwind3Kernel!(FThEDMF,ThEDMF,U,wEDMF,DS,dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeIEDMF)
    if NumTr > 0
      KAdvectionTrUpwind3Kernel!(FTrEDMF,TrEDMF,U,wEDMF,DS,dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeIEDMFTr)
    end  
  end    
      
  Parallels.ExchangeData3DRecvGPU!(F,Exchange)

  if Global.Model.Turbulence
    KV = Cache.KV
    KTkeSourceKernel! = TkeSourceKernel!(backend, groupG)
    KTkeSourceKernel!(Global.Model.TurbulenceSource,FTke,KV,U,dz,ndrange=ndrangeG)
  end 
  if Global.Model.VerticalDiffusionMom
    KVerticalDiffusionMomentumKernel! = VerticalDiffusionMomentumKernel!(backend,groupG)
    KV = Cache.KV
    KVerticalDiffusionMomentumKernel!(F,U,KV,dz,ndrange=ndrangeG)
  end 
  if Global.Model.SurfaceFluxMom
    CM = Global.SurfaceData.CM
    KSurfaceFluxMomentumKernel! = SurfaceFluxMomentumKernel!(backend,groupS)
    KSurfaceFluxMomentumKernel!(F,U,nSS,CM,dz,ndrange=ndrangeS)
  end  
  if Global.Model.VerticalDiffusion
    KVerticalDiffusionScalarKernel! = VerticalDiffusionScalarKernel!(backend,groupG)
    KV = Cache.KV
    @views FTh = F[:,:,5]
    @views Th = U[:,:,5]
    KVerticalDiffusionScalarKernel!(FTh,Th,Rho,KV,dz,ndrange=ndrangeG)
    if TkePos > 0
      KVerticalDiffusionScalarKernel!(FTke,Tke,Rho,KV,dz,ndrange=ndrangeG)
    end  
    for iT = 1 : NumTr
      @views KVerticalDiffusionScalarKernel!(FTr[:,:,iT],UTr[:,:,iT],Rho,KV,dz,ndrange=ndrangeG)
    end
  end
  if Global.Model.SurfaceFlux
    KSurfaceFluxScalarsKernel! = SurfaceFluxScalarsKernel!(backend,groupS)
    CT = Global.SurfaceData.CT
    CH = Global.SurfaceData.CH
    uStar = Global.SurfaceData.uStar
    TSurf = Global.SurfaceData.TS
    RhoVSurf = Global.SurfaceData.RhoVS
    SurfaceFluxRhs! = Global.Model.SurfaceFluxRhs
    KSurfaceFluxScalarsKernel!(SurfaceFluxRhs!,F,U,p,TSurf,RhoVSurf,uStar,CT,CH,dz,ndrange=ndrangeS)
  end
  if Global.Model.Forcing
    KForceKernel! = ForceKernel!(backend, groupG)
    KForceKernel!(Force,F,U,p,xS,ndrange=ndrangeG)  
  end  

  if Global.Model.Microphysics
    KMicrophysicsKernel! = MicrophysicsKernel!(backend, groupG)
    KMicrophysicsKernel!(MicrophysicsSource,F,U,p,ndrange=ndrangeG)
  end

  if Global.Model.Damping
    KDampKernel! = DampKernel!(backend, groupG)
    KDampKernel!(Damp,F,U,zP,ndrange=ndrangeG)
  end  

end

function FcnGPU!(F,U,FE,Metric,Phys,Cache,Exchange,Global,Param,Equation::Models.CompressibleDeep)

  backend = get_backend(F)
  FT = eltype(F)
  Glob = FE.Glob
  DS = FE.DS
  DW = FE.DW
  M = FE.M
  dXdxI = Metric.dXdxI
  X = Metric.X
  J = Metric.J
  NF = Global.Grid.NumFaces
  NBF = Global.Grid.NumFacesB
  @views dXdxI_B = dXdxI[:,:,:,:,:,1:NBF]
  @views J_B = J[:,:,:,1:NBF]
  @views X_B = X[:,:,:,:,1:NBF]
  @views Glob_B = Glob[:,1:NBF]

  @views dXdx_I = dXdx[:,:,:,:,:,NBF+1:NF]
  @views dXdxI_I = dXdxI[:,:,:,:,:,NBF+1:NF]
  @views J_I = J[:,:,:,NBF+1:NF]
  @views X_I = X[:,:,:,:,NBF+1:NF]
  @views Glob_I = Glob[:,NBF+1:NF]
  xS = Metric.xS  
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
  @views CacheF = Temp1[:,:,1:5]
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
  KHyperViscKernel!(CacheF,MRho,U,DS,DW,dXdxI,J,M,Glob,ndrange=ndrangeB)
  KernelAbstractions.synchronize(backend)
  for iT = 1 : NumTr
    @views CacheTr = Temp1[:,:,iT + 6]  
    @views KHyperViscTracerKernel!(CacheTr,U[:,:,iT+NumV],Rho,DS,DW,dXdxI,J,M,Glob,ndrange=ndrangeB)
    KernelAbstractions.synchronize(backend)
  end  
  @views KHyperViscWKernel!(Cachew,U[:,:,4],DS,DW,dXdxI,J,M,Glob,ndrange=ndrangeB)
  KernelAbstractions.synchronize(backend)
  Parallels.ExchangeData3DSendGPU(CacheFF,Exchange)

  KHyperViscKernel!(CacheF,MRho,U,DS,DW,dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeI)
  KernelAbstractions.synchronize(backend)
  for iT = 1 : NumTr
    @views CacheTr = Temp1[:,:,iT + 6]  
    @views KHyperViscTracerKernel!(CacheTr,U[:,:,iT+NumV],Rho,DS,DW,dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeI)
    KernelAbstractions.synchronize(backend)
  end  
  @views KHyperViscWKernel!(Cachew,U[:,:,4],DS,DW,dXdxI_I,J_I,M,Glob_I,ndrange=ndrangeI)
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
  @views KHyperViscWKoeffKernel!(F[:,:,4],Cachew,DS,DW,dXdxI,J,M,Glob,KoeffDivW,ndrange=ndrangeB)
  KernelAbstractions.synchronize(backend)
  KGradDeepKernel!(F,U,p,DS,dXdxI,J,X,M,MRho,Glob,Phys,ndrange=ndrangeB)
  KernelAbstractions.synchronize(backend)
  KMomentumDeepCoriolisKernel!(F,U,DS,dXdxI,J,X,MRho,M,Glob,Phys,ndrange=ndrangeB)
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
  @views KHyperViscWKoeffKernel!(F[:,:,4],Cachew,DS,DW,dXdxI_I,J_I,M,Glob,KoeffDivW,ndrange=ndrangeI)
  KernelAbstractions.synchronize(backend)
  KGradDeepKernel!(F,U,p,DS,dXdxI_I,J_I,X_I,M,MRho,Glob_I,Phys,ndrange=ndrangeI)
  KernelAbstractions.synchronize(backend)
  KMomentumDeepCoriolisKernel!(F,U,DS,dXdxI_I,J_I,X_I,MRho,M,Glob_I,Phys,ndrange=ndrangeI)
  KernelAbstractions.synchronize(backend)
  KRhoGradKinKernel!(F,U,DS,dXdxI_I,J_I,M,MRho,Glob,ndrange=ndrangeI)
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
    @show maximum(xS[2,:])
    @show minimum(xS[2,:])
    KForceKernel!(Force,F,U,p,xS,ndrange=ndrangeG)  
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

