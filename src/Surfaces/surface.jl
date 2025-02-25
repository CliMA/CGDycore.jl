abstract type SurfaceValues end
abstract type SurfaceFluxValues end

Base.@kwdef struct DefaultSurface <: SurfaceValues end

function (::DefaultSurface)(Phys,Param,uPos,vPos,wPos)
  @inline function SurfaceValues(SD,xS,U,p)
    SD[TSurfPos] = Param.TSurf
  end  
  return SurfaceValues
end  

Base.@kwdef struct HeldSuarezMoistSurface <: SurfaceValues end

function (::HeldSuarezMoistSurface)(Phys,Param,uPos,vPos,wPos)
  @inline function SurfaceValues(SD,xS,U,p)
    FT = eltype(xS)
    Lon = xS[1]
    Lat = xS[2]
    TSurf = Param.DeltaTS * exp(-FT(0.5) * Lat^2 / Param.DeltaLat^2) + Param.TSMin
    p_vs = Thermodynamics.fpws(TSurf,Phys)
    SD[RhoVSurfPos] = p_vs / (Phys.Rv * TSurf)
    SD[TSurfPos] = TSurf
  end  
  return SurfaceValues
end  
Base.@kwdef struct HeldSuarezMoistSurfaceFlux <: SurfaceFluxValues end

function (::HeldSuarezMoistSurfaceFlux)(Phys,Param,uPos,vPos,wPos)
  @inline function SurfaceFluxValues(SD,z,U,p,nS,z0M,z0H,LandClass)
    FT = eltype(U)
    SD[uStarPos] = uStarCoefficientGPU(U[uPos],U[vPos],U[wPos],nS)
    SD[CMPos] = FT(Param.CM)
    SD[CTPos] = FT(Param.CE)
    SD[CHPos] = FT(Param.CH)
  end
  return SurfaceFluxValues
end  

Base.@kwdef struct HeldSuarezDrySurface <: SurfaceValues end

function (::HeldSuarezDrySurface)(Phys,Param,uPos,vPos,wPos)
  @inline function SurfaceValues(x,U,p)
    FT = eltype(x)
    (Lon,Lat,R)= Grids.cart2sphere(x[1],x[2],x[3])
    TSurf = Param.DeltaTS * exp(-FT(0.5) * Lat^2 / Param.DeltaLat^2) + Param.TSMin
    return TSurf, FT(0)
  end
end

Base.@kwdef struct HeldSuarezDrySurfaceFlux <: SurfaceFluxValues end

function (::HeldSuarezDrySurfaceFlux)(Phys,Param,uPos,vPos,wPos)
  @inline function SurfaceFluxValues(z,U,p,nS,TS,z0M,z0H,LandClass)
    FT = eltype(U)
    uStar = uStarCoefficientGPU(U[uPos],U[vPos],U[wPos],nS)
    CM = FT(Param.CM)
    CT = FT(Param.CE)
    CH = FT(0)
    RiBSurf = FT(0)
    return uStar, CM, CT, CH
  end
  return SurfaceFluxValues
end

Base.@kwdef struct FriersonSurface <: SurfaceValues end

function (::FriersonSurface)(Phys,Param,RhoPos,uPos,vPos,wPos,ThPos)
  @inline function SurfaceValues!(SD,xS,U,p)
    FT = eltype(SD)
    Lat = xS[2] 
    SD[TSurfPos] = Param.DeltaTS * exp(-FT(0.5) * Lat^2 / Param.DeltaLat^2) + Param.TSMin
    SD[RhoVSurfPos] = FT(0)
  end
end

Base.@kwdef struct FriersonSurfaceFlux <: SurfaceFluxValues end

function (::FriersonSurfaceFlux)(Phys,Param,RhoPos,uPos,vPos,wPos,ThPos)
  @inline function SurfaceFluxValues!(SD,z,U,p,nS,z0M,z0H,LandClass)
    FT = eltype(U)
    norm_uh = U[uPos]^2 + U[vPos]^2
    Th = U[ThPos] / U[RhoPos]
    SD[RiBSurfPos] = Phys.Grav * z * (Th - SD[TS]) / SD[TS] / norm_uh
    SD[uStarPos] = uStarCoefficientGPU(U[uPos],U[vPos],U[wPos],nS)
    if RiBSurf < FT(0)
      SD[CMPos] = Phys.Karman^2 / log(z / Param.z_0)^2  
    elseif RiBSurf < Param.Ri_C
      SD[CMPos] = Phys.Karman^2 / log(z / Param.z_0)^2 * (FT(1) - SD[RiBSurf] / Param.Ri_C)
    else
      SD[CMPos] = FT(0)
    end  
    SD[CTPos] = SD[CMPos]
    SD[CHPos] = SD[CMPos]
  end
  return SurfaceFluxValues!
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

Base.@kwdef struct MOSurfaceFlux <: SurfaceFluxValues end

function (::MOSurfaceFlux)(uf,Phys,RhoPos,uPos,vPos,wPos,ThPos)
  @inline function SurfaceFluxValues(SD,z,U,p,nSS,z0M,z0H,LandClass)
    FT = eltype(U)
    TS = SD[TSurfPos]
    uStar = uStarCoefficientGPU(U[uPos],U[vPos],U[wPos],nSS)
    theta = U[ThPos] / U[RhoPos]
    CM, CT = MOSTIteration(uf,z0M,z0H,z,uStar,theta,TS,LandClass,Phys)
    SD[uStarPos] = uStar
    SD[CMPos] = CM
    SD[CTPos] = CT
    SD[CHPos] = CT
  end  
  return SurfaceFluxValues
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

function SurfaceFlux(Phys,Param,ThPos,RhoPos,RhoVPos)
  @inline function SurfaceFlux!(FU,U,p,dz,SD)
    FT = eltype(U)
    Rho = U[RhoPos]
    RhoTh = U[ThPos]
    RhoV = U[RhoVPos]
    RhoD = Rho - RhoV
    Rm = Phys.Rd * RhoD + Phys.Rv * RhoV
    Cpml = Phys.Cpd * RhoD + Phys.Cpv * RhoV
    T = p / Rm
    LatFlux = -SD[CHPos] * SD[uStarPos] * (RhoV - SD[RhoVSurfPos])
    SensFlux = -SD[CTPos] * SD[uStarPos] * (T - SD[TSurfPos])
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
  @inline function SurfaceFlux!(FU,U,p,dz,SD)
    FT = eltype(U)
    Rho = U[RhoPos]
    RhoTh = U[ThPos]
    RhoD = Rho
    Rm = Phys.Rd * RhoD 
    Cpml = Phys.Cpd * RhoD 
    T = p / Rm
    SensFlux = -SD[CTPos] * SD[uStarPos] * (T - SD[TSurfPos])
    PrePi=(p / Phys.p0)^(Rm / Cpml)
    FRhoTh = RhoTh * SensFlux / T 
    FU[ThPos] += FRhoTh  / dz
  end  
  return SurfaceFlux!
end     
