abstract type SurfaceValues end

Base.@kwdef struct HeldSuarezMoistSurface <: SurfaceValues end

function (::HeldSuarezMoistSurface)(Phys,Param,uPos,vPos,wPos,TS,RhoVS,CM,CT,CH,uStar)
  @inline function SurfaceValues(SD,xS,U,p)
    FT = eltype(xS)
    Lon = xS[1]
    Lat = xS[2]
    TSurf = Param.DeltaTS * exp(-FT(0.5) * Lat^2 / Param.DeltaLat^2) + Param.TSMin
    p_vs = Thermodynamics.fpvs(TSurf,Phys.T0)
    SD[RhoVS] = p_vs / (Phys.Rv * TSurf)
    SD[TS] = TSurf
  end  
  @inline function SurfaceFluxValues(SD,z,U,p,nS,z0M,z0H,LandClass)
    FT = eltype(U)
    SD[uStar] = uStarCoefficientGPU(U[uPos],U[vPos],U[wPos],nS)
    SD[CM] = FT(Param.CM)
    SD[CT] = FT(Param.CE)
    SD[CH] = FT(Param.CH)
  end
  return SurfaceValues, SurfaceFluxValues
end  

Base.@kwdef struct HeldSuarezDrySurface <: SurfaceValues end

function (::HeldSuarezDrySurface)(Phys,Param,uPos,vPos,wPos)
  @inline function SurfaceValues(x,U,p)
    FT = eltype(x)
    (Lon,Lat,R)= Grids.cart2sphere(x[1],x[2],x[3])
    TSurf = Param.DeltaTS * exp(-FT(0.5) * Lat^2 / Param.DeltaLat^2) + Param.TSMin
    return TSurf, FT(0)
  end
  @inline function SurfaceFluxValues(z,U,p,nS,TS,z0M,z0H,LandClass)
    FT = eltype(U)
    uStar = uStarCoefficientGPU(U[uPos],U[vPos],U[wPos],nS)
    CM = FT(Param.CM)
    CT = FT(Param.CE)
    CH = FT(0)
    RiBSurf = FT(0)
    return uStar, CM, CT, CH
  end
  return SurfaceValues, SurfaceFluxValues
end

Base.@kwdef struct FriersonSurface <: SurfaceValues end

function (::FriersonSurface)(Phys,Param,RhoPos,uPos,vPos,wPos,ThPos,TS,RhoVS,CM,CT,CH,RibSurf)
  @inline function SurfaceValues!(SD,xS,U,p)
    FT = eltype(SD)
    Lat = xS[2] 
    SD[TS] = Param.DeltaTS * exp(-FT(0.5) * Lat^2 / Param.DeltaLat^2) + Param.TSMin
    SD[RhoVS] = FT(0)
  end
  @inline function SurfaceFluxValues!(SD,z,U,p,nS,z0M,z0H,LandClass)
    FT = eltype(U)
    norm_uh = U[uPos]^2 + U[vPos]^2
    Th = U[ThPos] / U[RhoPos]
    SD[RiBSurf] = Phys.Grav * z * (Th - SD[TS]) / SD[TS] / norm_uh
    uStar = uStarCoefficientGPU(U[uPos],U[vPos],U[wPos],nS)
    if RiBSurf < FT(0)
      SD[CM] = Phys.Karman^2 / log(z / Param.z_0)^2  
    elseif RiBSurf < Param.Ri_C
      SD[CM] = Phys.Karman^2 / log(z / Param.z_0)^2 * (FT(1) - SD[RiBSurf] / Param.Ri_C)
    else
      SD[CM] = FT(0)
    end  
    SD[CT] = SD[CM]
    SD[CH] = SD[CM]
  end
  return SurfaceValues!, SurfaceFluxValues!
end


@inline function uStarCoefficientGPU(v1,v2,w,nS)
# Computation norm_v_a
# |v_a| = |v - n(n*v)| = sqrt(v*v -(n*v)^2)
  wS = -(nS[1]* v1 + nS[2] * v2) / nS[3]
  wC = eltype(v1)(0.5) * (wS + w)
  nU = nS[1] * v1 + nS[2] * v2 + nS[3] * wC
  uStar = sqrt((v1 - nS[1] * nU) * (v1 - nS[1] * nU) +
    (v2 - nS[2] * nU) * (v2 - nS[2] * nU) +
    (wC - nS[3] * nU) * (wC - nS[3] * nU))
  uStar = max(uStar, 0.1)
  return uStar
end

Base.@kwdef struct MOSurface <: SurfaceValues end

function (::MOSurface)(uf,Phys,RhoPos,uPos,vPos,wPos,ThPos)
  @inline function SurfaceValues(x,U,p)
    FT = eltype(x)
    return FT(300), FT(0)
  end
  @inline function SurfaceFluxValues(z,U,p,nSS,TS,z0M,z0H,LandClass)
    FT = eltype(U)
    uStar = uStarCoefficientGPU(U[uPos],U[vPos],U[wPos],nSS)
    theta = U[ThPos] / U[RhoPos]
    CM, CT = Surfaces.MOSTIteration(uf,z0M,z0H,z,uStar,theta,TS,LandClass,Phys)
    return uStar, CM, CT, CT
  end  
  return SurfaceValues, SurfaceFluxValues
end  



@kernel inbounds = true function SurfaceFluxDataKernel!(SurfaceFluxValues!,SurfaceData,@Const(U),@Const(p),@Const(dz),
  @Const(nSS),@Const(z0M),@Const(z0H),@Const(LandClass))

  IC, = @index(Global, NTuple)

  NumG = @uniform @ndrange()[1]

  if IC <= NumG
    SurfaceFluxValues!(view(SurfaceData,:,IC),dz[1,IC],
      view(U,1,IC,:),p[1,IC], view(nSS,:,IC),
      z0M[IC],z0H[IC],LandClass[IC])
  end
end

function SurfaceFluxData!(U,p,T,PotT,dz,nSS,SurfaceData,LandUseData,Model,NumberThreadGPU)
  z0M = LandUseData.z0M
  z0H = LandUseData.z0H
  LandClass = LandUseData.LandClass
  backend = get_backend(U)
  NumG = size(U,2)
  groupS = (max(div(NumG,NumberThreadGPU),1))
  ndrangeS = (NumG)
  KSurfaceFluxDataKernel! = SurfaceFluxDataKernel!(backend,groupS)
  KSurfaceFluxDataKernel!(Model.SurfaceFluxValues,SurfaceData,U,p,dz,nSS,
    z0M,z0H,LandClass,ndrange=ndrangeS)
  KernelAbstractions.synchronize(backend)
end

@kernel inbounds = true function SurfaceDataKernel!(SurfaceValues!,SurfaceData,@Const(U),@Const(p),@Const(xS),@Const(Glob))

  IC, = @index(Global, NTuple)

  NumG = @uniform @ndrange()[1]

  if IC <= NumG
    SurfaceValues!(view(SurfaceData,:,IC),view(xS,:,IC),view(U,1,IC,:),p[1,IC])
  end
end

function SurfaceData!(U,p,xS,Glob,SurfaceData,Model,NumberThreadGPU)
  backend = get_backend(U)
  NumG = size(U,2)
  groupS = (max(div(NumG,NumberThreadGPU),1))
  ndrangeS = (NumG)
  KSurfaceDataKernel! = SurfaceDataKernel!(backend,groupS)
  KSurfaceDataKernel!(Model.SurfaceValues,SurfaceData,U,p,xS,Glob,ndrange=ndrangeS)
  KernelAbstractions.synchronize(backend)
end

function SurfaceFlux(Phys,Param,ThPos,RhoPos,RhoVPos,uStar,CT,CH,TSurf,RhoVSurf)
  @inline function SurfaceFlux!(FU,U,p,dz,SD)
    FT = eltype(U)
    Rho = U[RhoPos]
    RhoTh = U[ThPos]
    RhoV = U[RhoVPos]
    RhoD = Rho - RhoV
    Rm = Phys.Rd * RhoD + Phys.Rv * RhoV
    Cpml = Phys.Cpd * RhoD + Phys.Cpv * RhoV
    T = p / Rm
    LatFlux = -SD[CT] * SD[uStar] * (RhoV - SD[RhoVSurf])
    SensFlux = -SD[CH] * SD[uStar] * (T - SD[TSurf])
    FRho = LatFlux
    FRhoV = LatFlux
    PrePi=(p / Phys.p0)^(Rm / Cpml)
    FRhoTh = RhoTh * (SensFlux / T + ((Phys.Rv / Rm) - FT(1) / Rho -
      log(PrePi)*(Phys.Rv / Rm - Phys.Cpv / Cpml)) *  LatFlux)
    FU[RhoPos] += FRho / dz
    FU[ThPos] += FRhoTh  / dz
    FU[RhoVPos] += FRhoV  / dz
  end  
  return SurfaceFlux!
end     
function SurfaceFlux(Phys,Param,ThPos,RhoPos)
  @inline function SurfaceFlux!(FU,U,p,dz,uStar,CT,CH,TSurf,RhoVSurf)
    FT = eltype(U)
    Rho = U[RhoPos]
    RhoTh = U[ThPos]
    RhoD = Rho
    Rm = Phys.Rd * RhoD 
    Cpml = Phys.Cpd * RhoD 
    T = p / Rm
    SensFlux = -CH * uStar * (T - TSurf)
    PrePi=(p / Phys.p0)^(Rm / Cpml)
    FRhoTh = RhoTh * SensFlux / T 
    FU[ThPos] += FRhoTh  / dz
  end  
  return SurfaceFlux!
end     
