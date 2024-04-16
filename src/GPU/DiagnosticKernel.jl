@kernel inbounds = true function StrainRate!(Strain,@Const(U),@Const(dz))
  Iz,IC = @index(Global, NTuple)

  NumG = @uniform @ndrange()[2]

  if IC <= NumG
    SU = eltype(Strain)(0)
    SL = eltype(Strain)(0)
    if Iz < Nz
      dudzU = (U[Iz+1,IC,2] - U[Iz,IC,2]) / (dz[Iz+1] - dz[Iz])
      dvdzU = (U[Iz+1,IC,3] - U[Iz,IC,3]) / (dz[Iz+1] - dz[Iz])
      SU = sqrt(dudzU * dudzU + dvdzU * dvdzU)
    end
    if Iz > 1
      dudzL = (U[Iz,IC,2] - U[Iz-1,IC,2]) / (dz[Iz] - dz[Iz-1])
      dvdzL = (U[Iz,IC,3] - U[Iz-1,IC,3]) / (dz[Iz] - dz[Iz-1])
      SL = sqrt(dudzL * dudzL + dvdzL * dvdzL)
    end   
    Strain[Iz,IC] = eltype(Strain)(0.5) * (SU + SL)
  end
end


@kernel inbounds = true function TkeSourceKernel!(Source,FTke,KV,@Const(U),@Const(dz))
  Iz,IC = @index(Global, NTuple)

  Nz = @uniform @ndrange()[1]
  NumG = @uniform @ndrange()[2]

  if IC <= NumG && Iz < Nz
    @views FLoc, KLoc = Source(U[Iz+1,IC,:],U[Iz,IC,:],dz[Iz+1,IC],dz[Iz,IC])  
    @atomic :monotonic FTke[Iz,IC] += eltype(FTke)(.5) * FLoc
    @atomic :monotonic FTke[Iz+1,IC] += eltype(FTke)(.5) * FLoc
    KV[Iz,IC] = KLoc
  end  
end


@kernel inbounds = true function EddyCoefficientTKEKernel!(Eddy,K,@Const(Rho),@Const(Tke),@Const(Glob))
  Iz,IC = @index(Global, NTuple)

  Nz = @uniform @ndrange()[1]
  NumG = @uniform @ndrange()[2]

  if IC <= NumG
    @views K[Iz,IC] = Eddy(Tke[Iz,IC],  )
  end
end

@kernel inbounds = true function PressureKernel!(Pressure,p,@Const(U))
  Iz,IC = @index(Global, NTuple)

  NumG = @uniform @ndrange()[2]

  if IC <= NumG
    p[Iz,IC] = Pressure(view(U,Iz,IC,:))
  end
end

@kernel inbounds = true function PressureCKernel!(Pressure,p,@Const(U),@Const(Glob),@Const(JJ),@Const(M))
  ID,Iz,IF = @index(Global, NTuple)

  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  if Iz <= Nz && IF <= NF
    ind = Glob[ID,IF]  
    @atomic :monotonic p[Iz,ind] += Pressure(view(U,Iz,ind,:)) * JJ[ID,Iz,IF] / M[Iz,ind]
  end
end

@kernel inbounds = true function SurfaceKernel!(Surface,@Const(U),@Const(p),@Const(X),
  @Const(dXdxI),@Const(nS),@Const(Glob),SurfaceValues,SurfaceData)

  ID,IF = @index(Global, NTuple)

  NumF = @uniform @ndrange()[2]

  if IF <= NumF
    ind = Glob[ID,IF]
    xS = SVector{3}(X[ID,1,1,1,IF], X[ID,1,2,1,IF], X[ID,1,3,1,IF])
    dz = sqrt(X[ID,2,1,1,IF]^2 + X[ID,2,2,1,IF]^2 + X[ID,2,3,1,IF]^2) -
      sqrt(X[ID,1,1,1,IF]^2 + X[ID,1,2,1,IF]^2 + X[ID,1,3,1,IF]^2)
    SurfaceValues(xS,view(U,1,ind,:),p[1,ind],Surface[ID,IF])
    SurfaceData(dz,view(U,1,ind,:),p[1,ind],
     view(dXdxI,3,:,1,ID,1,IF),view(nS,ID,:,IF),Surface[ID,IF])
  end
end

@kernel inbounds = true function EddyCoefficientKernel!(Eddy,K,@Const(U),@Const(uStar),@Const(p),@Const(dz),@Const(Glob))
  Iz,IC = @index(Global, NTuple)

  Nz = @uniform @ndrange()[1]
  NumG = @uniform @ndrange()[2]

  if IC <= NumG
    LenScale = 100.0  
    @views K[Iz,IC] = Eddy(U[Iz,IC,:],uStar[IC],p[Iz,IC],dz[1,IC],LenScale) 
  end
end

@kernel inbounds = true function EddyCoefficientCKernel!(Eddy,K,@Const(Rho),@Const(uStar),@Const(p),@Const(dz),@Const(Glob))
  ID,Iz,IF = @index(Global, NTuple)

  Nz = @uniform @ndrange()[2]
  NumF = @uniform @ndrange()[3]

  if Iz <= Nz && IF <= NumF
    ind = Glob[ID,IF]
    p = Pressure(view(U,Iz,ind,:)) 
    K[Iz,ind] += Eddy(uStar[ID,IF],p,dz[1,ind]) * Rho[Iz,ind] * JJ[ID,Iz,IF] / M[Iz,ind]
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
  xS = Metric.xS
  nS = Metric.nS
  nSS = Metric.nSS
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
  NDoFG = min(div(NumberThreadGPU,Nz),NumG)
  groupG = (Nz, NDoFG)
  ndrangeG = (Nz, NumG)
  @views p = Cache.AuxG[:,:,1]
  @views KV = Cache.KV
  @views Rho = U[:,:,1]
  dz = Metric.dz
  Eddy = Global.Model.Eddy
  SurfaceData = Global.SurfaceData
  uStar = SurfaceData.uStar
  LandUseData = Global.LandUseData
  Model = Global.Model
  Pressure = Global.Model.Pressure

  KPressureKernel! = PressureKernel!(backend,group)
  KPressureKernel!(Pressure,p,U,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)

  if Global.Model.SurfaceFlux || Global.Model.SurfaceFluxMom
    Surfaces.SurfaceData!(U,p,xS,Glob,SurfaceData,Model,NumberThreadGPU)  
    Surfaces.SurfaceFluxData!(U,p,dz,nSS,SurfaceData,LandUseData,Model,NumberThreadGPU)  
  end  
  if !Global.Model.Turbulence
    if Global.Model.VerticalDiffusion || Global.Model.VerticalDiffusionMom
      if Global.Model.SurfaceFlux || Global.Model.SurfaceFluxMom 
        KEddyCoefficientKernel! = EddyCoefficientKernel!(backend,groupG)
        KEddyCoefficientKernel!(Eddy,KV,U,uStar,p,dz,Glob,ndrange=ndrangeG)
        KernelAbstractions.synchronize(backend)
      else    
        NFG = min(div(NumberThreadGPU,N*N),NumF)
        Surfaces.SurfaceData!(U,p,xS,Glob,SurfaceData,Model,NumberThreadGPU)  
        Surfaces.SurfaceFluxData!(U,p,dz,nSS,SurfaceData,LandUseData,Model,NumberThreadGPU)  
        KEddyCoefficientKernel! = EddyCoefficientKernel!(backend,groupG)
        @show "EddyCoefficientKernel",groupK,ndrangeK
        KEddyCoefficientKernel!(Eddy,KV,U,uStar,p,dz,Glob,ndrange=ndrangeG)
        KernelAbstractions.synchronize(backend)
      end
    end  
  end
end

