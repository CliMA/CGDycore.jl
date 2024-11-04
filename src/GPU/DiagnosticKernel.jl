#=
@kernel inbounds = true function KineticKernel!(KE,@Const(U),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob))

# gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  KinF = @localmem eltype(F) (N,N,ColumnTilesDim)

  ID = I + (J - 1) * N
  ind = Glob[ID,IF]

  if Iz <= Nz
    KinF[I,J,iz] = eltype(F)(0.5) * (U[Iz,ind,2]^2 + U[Iz,ind,3]^2)
    if iz == 1 && Iz == 1
      wCol = -(dXdxI[3,1,1,ID,1,IF] * U[Iz,ind,2] +
        dXdxI[3,2,1,ID,1,IF] * U[Iz,ind,3]) / dXdxI[3,3,1,ID,1,IF]
      KinF[I,J,iz] +=  eltype(F)(0.5) * wCol^2
    elseif iz == 1
      KinF[I,J,iz] +=  eltype(F)(0.5) * U[Iz-1,ind,4]^2
    else
      KinF[I,J,1,iz] +=  eltype(F)(0.5) * U[Iz-1,ind,4]^2
    end
    if iz == ColumnTilesDim && Iz < Nz
      KinF[I,J,iz] +=  eltype(F)(0.5) * U[Iz,ind,4]^2
    end
  end
  if Iz <= Nz && IF <= NF
    ind = Glob[ID,IF]  
    @atomic :monotonic KE[Iz,ind] += KinF[I,J,iz] * JJ[ID,Iz,IF] / M[Iz,ind]
  end
end  
=#

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
    if Iz == 1
      @atomic :monotonic FTke[Iz,IC] += eltype(FTke)(.5) * FLoc  
    end  
    if Iz == Nz
      @atomic :monotonic FTke[Iz+1,IC] += eltype(FTke)(.5) * FLoc  
    end  
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

@kernel inbounds = true function PressureKernel!(Pressure,p,T,PotT,@Const(U),@Const(nS),@Const(zP))
  Iz,IC = @index(Global, NTuple)

  Nz = @uniform @ndrange()[1]
  NumG = @uniform @ndrange()[2]

  if IC <= NumG
    if Iz == 1
      wL = wCol = -(nS[1,Iz] * U[Iz,IC,2] + nS[1,Iz] * U[Iz,IC,3]) / nS[3,Iz]
    else
      wL = U[Iz,IC,4]
    end
    if Iz == Nz
      wR = eltype(U)(0)
    else
      wR = U[Iz+1,IC,4] 
    end  
    p[Iz,IC], T[Iz,IC], PotT[Iz,IC] = Pressure(view(U,Iz,IC,:),wL,wR,zP[Iz,IC])
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
  @views T = Cache.AuxG[:,:,2]
  @views PotT = Cache.AuxG[:,:,3]
  @views KV = Cache.KV
  @views Rho = U[:,:,1]
  dz = Metric.dz
  zP = Metric.zP
  Eddy = Global.Model.Eddy
  SurfaceData = Global.SurfaceData
  uStar = SurfaceData.uStar
  LandUseData = Global.LandUseData
  Model = Global.Model
  Pressure = Global.Model.Pressure

  KPressureKernel! = PressureKernel!(backend,group)
  KPressureKernel!(Pressure,p,T,PotT,U,nSS,zP,ndrange=ndrange)

  if Global.Model.SurfaceFlux || Global.Model.SurfaceFluxMom
    Surfaces.SurfaceData!(U,p,xS,Glob,SurfaceData,Model,NumberThreadGPU)  
    Surfaces.SurfaceFluxData!(U,p,T,PotT,dz,nSS,SurfaceData,LandUseData,Model,NumberThreadGPU)  
  end  
  if !Global.Model.Turbulence
    if Global.Model.VerticalDiffusion || Global.Model.VerticalDiffusionMom
      if Global.Model.SurfaceFlux || Global.Model.SurfaceFluxMom 
        KEddyCoefficientKernel! = EddyCoefficientKernel!(backend,groupG)
        KEddyCoefficientKernel!(Eddy,KV,U,uStar,p,dz,Glob,ndrange=ndrangeG)
      else    
        NFG = min(div(NumberThreadGPU,N*N),NumF)
        Surfaces.SurfaceData!(U,p,xS,Glob,SurfaceData,Model,NumberThreadGPU)  
        Surfaces.SurfaceFluxData!(U,p,dz,nSS,SurfaceData,LandUseData,Model,NumberThreadGPU)  
        KEddyCoefficientKernel! = EddyCoefficientKernel!(backend,groupG)
        @show "EddyCoefficientKernel",groupK,ndrangeK
        KEddyCoefficientKernel!(Eddy,KV,U,uStar,p,dz,Glob,ndrange=ndrangeG)
      end
    end  
  end
end

