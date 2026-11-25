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
    KinF[I,J,iz] = eltype(F)(0.5) * (U[1,Iz,ind,2]^2 + U[1,Iz,ind,3]^2)
    if iz == 1 && Iz == 1
      wCol = -(dXdxI[3,1,1,ID,1,IF] * U[1,Iz,ind,2] +
        dXdxI[3,2,1,ID,1,IF] * U[1,Iz,ind,3]) / dXdxI[3,3,1,ID,1,IF]
      KinF[I,J,iz] +=  eltype(F)(0.5) * wCol^2
    elseif iz == 1
      KinF[I,J,iz] +=  eltype(F)(0.5) * U[1,Iz-1,ind,4]^2
    else
      KinF[I,J,1,iz] +=  eltype(F)(0.5) * U[1,Iz-1,ind,4]^2
    end
    if iz == ColumnTilesDim && Iz < Nz
      KinF[I,J,iz] +=  eltype(F)(0.5) * U[1,Iz,ind,4]^2
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
      dudzU = (U[1,Iz+1,IC,2] - U[1,Iz,IC,2]) / (dz[Iz+1] - dz[Iz])
      dvdzU = (U[1,Iz+1,IC,3] - U[1,Iz,IC,3]) / (dz[Iz+1] - dz[Iz])
      SU = sqrt(dudzU * dudzU + dvdzU * dvdzU)
    end
    if Iz > 1
      dudzL = (U[1,Iz,IC,2] - U[1,Iz-1,IC,2]) / (dz[Iz] - dz[Iz-1])
      dvdzL = (U[1,Iz,IC,3] - U[1,Iz-1,IC,3]) / (dz[Iz] - dz[Iz-1])
      SL = sqrt(dudzL * dudzL + dvdzL * dvdzL)
    end   
    Strain[Iz,IC] = eltype(Strain)(0.5) * (SU + SL)
  end
end


# TKE closure along a column (one thread per column). Writes the source of rho*e, the eddy
# diffusivity rho*K at the interfaces Iz+1/2 (KV[1,Iz,IC], Iz < Nz) and the derivative of the
# sinks for the implicit solve (DTke). The surface production is u*^3/(kappa z1) plus the
# buoyancy flux of the surface heat flux; without surface data the lowest interface value is used.
# virtual potential temperature (unsaturated), buoyancy of the TKE closure
@inline function ThetaV(U,Iz,IC,RhoVPos,RhoCPos,Eps1)
  FT = eltype(U)
  Rho = U[1,Iz,IC,1]
  qv = RhoVPos > 0 ? U[1,Iz,IC,RhoVPos] / Rho : FT(0)
  qc = RhoCPos > 0 ? U[1,Iz,IC,RhoCPos] / Rho : FT(0)
  return U[1,Iz,IC,5] / Rho * (FT(1) + Eps1 * qv - qc)
end

@kernel inbounds = true function TkeColumnKernel!(Turb,FTke,KV,DTke,@Const(U),@Const(p),@Const(Surf),
  @Const(dz),@Const(zP),Grav,Rd,Cpd,Eps1,TkePos,RhoVPos,RhoCPos,SurfaceProd::Bool,SurfaceHeat::Bool)
  IC, = @index(Global, NTuple)

  Nz = @uniform size(U,2)
  NumG = @uniform @ndrange()[1]

  if IC <= NumG
    FT = eltype(FTke)
    zS = zP[1,IC] - FT(0.5) * dz[1,IC]
    RhoB = U[1,1,IC,1]
    uB = U[1,1,IC,2]
    vB = U[1,1,IC,3]
    ThB = ThetaV(U,1,IC,RhoVPos,RhoCPos,Eps1)
    eB = max(U[1,1,IC,TkePos] / RhoB, Turb.EMin)
    # production at the lower face of the current layer, N^2 below
    PB = FT(0)
    PT = FT(0)
    if SurfaceProd
      z1 = FT(0.5) * dz[1,IC]
      uStar = Surf[Surfaces.uStarPos,IC]
      PB = uStar * uStar * uStar / (Turb.Karman * z1)
      if SurfaceHeat
        # surface buoyancy flux g (w'theta'/T + Eps1 w'q_v'), theta_1 - theta_S as in the fluxes
        TB = p[1,1,IC] / (Rd * RhoB)
        pS = p[1,1,IC] + RhoB * Grav * z1
        TS = Surf[Surfaces.TSurfPos,IC] * (p[1,1,IC] / pS)^(Rd / Cpd)
        SensFlux = -Surf[Surfaces.CTPos,IC] * Surf[Surfaces.UaPos,IC] * (TB - TS)
        PB += Grav * SensFlux / TB
        if RhoVPos > 0
          LatFlux = -Surf[Surfaces.CHPos,IC] * Surf[Surfaces.UaPos,IC] *
            (U[1,1,IC,RhoVPos] - Surf[Surfaces.RhoVSurfPos,IC])
          PB += Grav * Eps1 * LatFlux / RhoB
        end
      end
    end
    N2B = FT(0)
    for Iz = 1 : Nz
      if Iz < Nz
        RhoT = U[1,Iz+1,IC,1]
        uT = U[1,Iz+1,IC,2]
        vT = U[1,Iz+1,IC,3]
        ThT = ThetaV(U,Iz+1,IC,RhoVPos,RhoCPos,Eps1)
        eT = max(U[1,Iz+1,IC,TkePos] / RhoT, Turb.EMin)
        dzF = FT(0.5) * (dz[Iz,IC] + dz[Iz+1,IC])
        S2 = ((uT - uB) * (uT - uB) + (vT - vB) * (vT - vB)) / (dzF * dzF)
        N2T = FT(2) * Grav * (ThT - ThB) / (dzF * (ThT + ThB))
        eF = FT(0.5) * (eB + eT)
        zF = FT(0.5) * (zP[Iz,IC] + zP[Iz+1,IC]) - zS
        lF = Models.MixingLength(Turb,zF,eF,N2T)
        KF = Models.EddyDiffusivity(Turb,lF,eF)
        KV[1,Iz,IC] = FT(0.5) * (RhoB + RhoT) * KF
        PT = KF * (S2 - N2T)
      else
        N2T = N2B
      end
      if Iz == 1 && !SurfaceProd
        PB = PT
      end
      if Iz == Nz
        PT = PB
      end
      # dissipation in the layer
      zC = zP[Iz,IC] - zS
      N2C = Iz == 1 ? N2T : FT(0.5) * (N2B + N2T)
      lC = Models.MixingLength(Turb,zC,eB,N2C)
      FTke[Iz,IC] += RhoB * (FT(0.5) * (PB + PT) - Models.Dissipation(Turb,lC,eB))
      DTke[Iz,IC] = Models.DissipationDerivative(Turb,lC,eB,N2C)
      if Iz < Nz
        PB = PT
        N2B = N2T
        RhoB = RhoT
        uB = uT
        vB = vT
        ThB = ThT
        eB = eT
      end
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

@kernel inbounds = true function PressureKernel!(Pressure!,Aux,@Const(U),@Const(nS),@Const(zP))
  K,Iz,IC = @index(Global, NTuple)

  Nz = @uniform @ndrange()[2]
  NumG = @uniform @ndrange()[3]

  if IC <= NumG
    if Iz == 1
      wL = -(nS[1,IC] * U[K,Iz,IC,2] + nS[2,IC] * U[K,Iz,IC,3]) / nS[3,IC]
    else
      wL = U[K,Iz,IC,4]
    end
    if Iz == Nz
      wR = eltype(U)(0)
    else
      wR = U[K,Iz+1,IC,4] 
    end  
    # start of the saturation adjustment, no temperature before the first evaluation
    TStart = Aux[K,Iz,IC,2] > eltype(U)(0) ? Aux[K,Iz,IC,2] : eltype(U)(300)
    Pressure!(view(Aux,K,Iz,IC,:),view(U,K,Iz,IC,:),wL,wR,zP[Iz,IC];T=TStart)
  end
end


@kernel inbounds = true function PressureKernelDG!(Pressure,p,T,PotT,@Const(U),@Const(nS),@Const(zP))
  Iz,K,IC = @index(Global, NTuple)

  Nz = @uniform @ndrange()[1]
  NumG = @uniform @ndrange()[3]

  if IC <= NumG
    if Iz == 1 && K == 1
      wL = wCol = -(nS[1,Iz] * U[1,Iz,K,IC,2] + nS[1,Iz] * U[1,Iz,K,IC,3]) / nS[3,Iz]
    else
      wL = U[1,Iz,K,IC,4]
    end
    if Iz == Nz
      wR = eltype(U)(0)
    else
      wR = U[1,Iz+1,k,IC,4] 
    end  
    p[Iz,K,IC], T[Iz,K,IC], PotT[Iz,K,IC] = Pressure(view(U,1,Iz,K,IC,:),wL,wR,zP[Iz,IC])
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

# rho*K at the interfaces Iz+1/2 (Iz < Nz), from the pressure of the interface
@kernel inbounds = true function EddyCoefficientKernel!(Eddy,KV,@Const(U),@Const(Surf),@Const(p),
  @Const(dz),@Const(Glob))
  K,Iz,IC = @index(Global, NTuple)

  Nz = @uniform @ndrange()[2]
  NumG = @uniform @ndrange()[3]

  if IC <= NumG && Iz < Nz
    LenScale = eltype(KV)(100)
    pF = eltype(KV)(0.5) * (p[K,Iz,IC] + p[K,Iz+1,IC])
    @views KV[K,Iz,IC] = Eddy(U[1,Iz,IC,:],Surf[:,IC],pF,dz[1,IC],LenScale) 
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

function FcnPrepare!(U,FE,Metric,Phys,Cache,Exchange,Global)
  backend = get_backend(U)
  dXdxI = Metric.dXdxI
  X = Metric.X
  xS = Metric.xS
  nS = Metric.nS
  nSS = Metric.nSS
  FT = eltype(U)
  M = size(U,1)
  Nz = size(U,2)
  NumG = size(U,3)
  Glob = FE.Glob
  NumF = size(Glob,2)
  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU

  NG = min(div(NumberThreadGPU,M*Nz),NumG)
  group = (M, Nz, NG)
  ndrange = (M, Nz, NumG)
  NDoFG = min(div(NumberThreadGPU,M*Nz),NumG)
  groupG = (M, Nz, NDoFG)
  ndrangeG = (M, Nz, NumG)
  Aux = Cache.Aux
  @views p = Aux[:,:,:,1]
  @views T = Aux[:,:,:,2]
  @views PotT = Aux[:,:,:,3]
  @views RhoP = Aux[:,:,:,5:end]
  KV = Cache.KV
  @views Rho = U[:,:,:,1]
  dz = Metric.dz
  zP = Metric.zP
  Eddy = Global.Model.Eddy
  SurfaceData = Global.SurfaceData
  LandUseData = Global.LandUseData
  Model = Global.Model
  Pressure = Global.Model.Pressure

  KPressureKernel! = PressureKernel!(backend,group)
  KPressureKernel!(Pressure,Aux,U,nSS,zP,ndrange=ndrange)

  if Global.Model.SurfaceFlux || Global.Model.SurfaceFluxMom
    Surfaces.SurfaceData!(U,p,xS,Glob,SurfaceData.Data,Model,NumberThreadGPU)  
    Surfaces.SurfaceFluxData!(U,p,T,PotT,dz,nSS,SurfaceData.Data,LandUseData,Model,NumberThreadGPU)  
  end  
  if !Global.Model.Turbulence
    if Global.Model.VerticalDiffusion || Global.Model.VerticalDiffusionMom
      if Global.Model.SurfaceFlux || Global.Model.SurfaceFluxMom 
        KEddyCoefficientKernel! = EddyCoefficientKernel!(backend,groupG)
        KEddyCoefficientKernel!(Eddy,KV,U,SurfaceData.Data,p,dz,Glob,ndrange=ndrangeG)
      else    
        Surfaces.SurfaceData!(U,p,xS,Glob,SurfaceData.Data,Model,NumberThreadGPU)  
        Surfaces.SurfaceFluxData!(U,p,T,PotT,dz,nSS,SurfaceData.Data,LandUseData,Model,NumberThreadGPU)  
        KEddyCoefficientKernel! = EddyCoefficientKernel!(backend,groupG)
        KEddyCoefficientKernel!(Eddy,KV,U,SurfaceData.Data,p,dz,Glob,ndrange=ndrangeG)
      end
    end  
  end
end

