abstract type SurfaceValues end

Base.@kwdef struct HeldSuarezMoistSurface <: SurfaceValues end


function (::HeldSuarezMoistSurface)(Phys,Param,uPos,vPos,wPos)
  @inline function SurfaceValues(x,U,p)
    FT = eltype(x)
    (Lon,Lat,R)= Grids.cart2sphere(x[1],x[2],x[3])
    TSurf = Param.DeltaTS * exp(-FT(0.5) * Lat^2 / Param.DeltaLat^2) + Param.TSMin
    p_vs = Thermodynamics.fpvs(TSurf,Phys.T0)
    RhoVSurf = p_vs / (Phys.Rv * TSurf)
    return TSurf, RhoVSurf
  end  
  @inline function SurfaceFluxValues(z,U,p,dXdxI,nS,TS,z0M,z0H,LandClass)
    FT = eltype(U)
    uStar = uStarCoefficientGPU(U[uPos],U[vPos],U[wPos],dXdxI,nS)
    CT = FT(Param.CE)
    CH = FT(Param.CH)
    return uStar, CT, CH
  end
  return SurfaceValues, SurfaceFluxValues
end  


@inline function uStarCoefficientGPU(v1,v2,w,dXdxI,nS)
# Computation norm_v_a
# |v_a| = |v - n(n*v)| = sqrt(v*v -(n*v)^2)
  @inbounds wS = -(dXdxI[1]* v1 + dXdxI[2] * v2) / dXdxI[3]
  wC = eltype(v1)(0.5) * (wS + w)
  @inbounds nU = nS[1] * v1 + nS[2] * v2 + nS[3] * wC
  @inbounds sqrt((v1 - nS[1] * nU) * (v1 - nS[1] * nU) +
    (v2 - nS[2] * nU) * (v2 - nS[2] * nU) +
    (wC - nS[3] * nU) * (wC - nS[3] * nU))
end

Base.@kwdef struct MOSurface <: SurfaceValues end

function (::MOSurface)(uf,Phys,RhoPos,uPos,vPos,wPos,ThPos)
  @inline function SurfaceValues(x,U,p)
    FT = eltype(x)
    return FT(300), FT(0)
  end
  @inline function SurfaceFluxValues(z,U,p,dXdxI,nS,TS,z0M,z0H,LandClass)
    FT = eltype(U)
    uStar = uStarCoefficientGPU(U[uPos],U[vPos],U[wPos],dXdxI,nS)
    theta = U[ThPos] / U[RhoPos]
    CT, CH = Surfaces.MOSTIteration(uf,z0M,z0H,z,uStar,theta,TS,LandClass,Phys)
    return uStar, CT, CH
  end  
  return SurfaceValues, SurfaceFluxValues
end  



@kernel function SurfaceFluxDataKernel!(SurfaceFluxValues,uStar,CT,CH,@Const(U),@Const(p),@Const(X),
  @Const(dXdxI),@Const(nS),@Const(Glob),@Const(TS),@Const(z0M),@Const(z0H),@Const(LandClass))

  ID,IF = @index(Global, NTuple)

  NumF = @uniform @ndrange()[2]

  if IF <= NumF
    @inbounds ind = Glob[ID,IF]
    @inbounds xS = SVector{3}(X[ID,1,1,1,IF], X[ID,1,2,1,IF], X[ID,1,3,1,IF])
    @inbounds dz = sqrt(X[ID,2,1,1,IF]^2 + X[ID,2,2,1,IF]^2 + X[ID,2,3,1,IF]^2) -
      sqrt(X[ID,1,1,1,IF]^2 + X[ID,1,2,1,IF]^2 + X[ID,1,3,1,IF]^2)
    uStar[ID,IF], CT[ID,IF], CH[ID,IF] = SurfaceFluxValues(dz,
      view(U,1,ind,:),p[1,ind], view(dXdxI,3,:,1,ID,1,IF),view(nS,ID,:,IF),
      TS[ID,IF],z0M[ID,IF],z0H[ID,IF],LandClass[ID,IF])
  end
end

function SurfaceFluxData!(U,p,X,dXdxI,nS,Glob,SurfaceData,LandUseData,Model,NumberThreadGPU)
  TS = SurfaceData.TS
  uStar = SurfaceData.uStar
  CT = SurfaceData.CT
  CH = SurfaceData.CH
  z0M = LandUseData.z0M
  z0H = LandUseData.z0H
  LandClass = LandUseData.LandClass
  backend = get_backend(TS)
  DoF, NumF = size(TS)
  NFG = min(div(NumberThreadGPU,DoF),NumF)
  groupS = (DoF, NFG)
  ndrangeS = (DoF, NumF)
  KSurfaceFluxDataKernel! = SurfaceFluxDataKernel!(backend,groupS)
  KSurfaceFluxDataKernel!(Model.SurfaceFluxValues,uStar,CT,CH,U,p,X,dXdxI,nS,Glob,TS,
    z0M,z0H,LandClass,ndrange=ndrangeS)
  KernelAbstractions.synchronize(backend)
end

@kernel function SurfaceDataKernel!(SurfaceValues,TS,RhoVS,@Const(U),@Const(p),@Const(X),@Const(Glob))

  ID,IF = @index(Global, NTuple)

  NumF = @uniform @ndrange()[2]

  if IF <= NumF
    @inbounds ind = Glob[ID,IF]
    @inbounds xS = SVector{3}(X[ID,1,1,1,IF], X[ID,1,2,1,IF], X[ID,1,3,1,IF])
    @inbounds TS, RhoVS =  SurfaceValues(xS,view(U,1,ind,:),p[1,ind])
  end
end

function SurfaceData!(U,p,X,Glob,SurfaceData,Model,NumberThreadGPU)
  TS = SurfaceData.TS
  RhoVS = SurfaceData.RhoVS
  backend = get_backend(TS)
  DoF, NumF = size(TS)
  NFG = min(div(NumberThreadGPU,DoF),NumF)
  groupS = (DoF, NFG)
  ndrangeS = (DoF, NumF)
  KSurfaceDataKernel! = SurfaceDataKernel!(backend,groupS)
  KSurfaceDataKernel!(Model.SurfaceValues,TS,RhoVS,U,p,X,Glob,ndrange=ndrangeS)
  KernelAbstractions.synchronize(backend)
end
