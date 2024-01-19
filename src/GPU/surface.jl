abstract type SurfaceValues end

Base.@kwdef struct HeldSuarezMoistSurface <: SurfaceValues end


function (::HeldSuarezMoistSurface)(Phys,Param,uPos,vPos,wPos)
  function SurfaceValues(x,U,p)
    FT = eltype(x)
    (Lon,Lat,R)= Grids.cart2sphere(x[1],x[2],x[3])
    TSurf = Param.DeltaTS * exp(-FT(0.5) * Lat^2 / Param.DeltaLat^2) + Param.TSMin
    p_vs = Thermodynamics.fpvs(TSurf,Phys.T0)
    RhoVSurf = p_vs / (Phys.Rv * TSurf)
    return TSurf, RhoVSurf
  end  
  function SurfaceData(U,p,dXdxI,nS)
    FT = eltype(U)
    uStar = uStarCoefficientGPU(U[uPos],U[vPos],U[wPos],dXdxI,nS)
    CT = FT(Param.CE)
    CH = FT(Param.CH)
    return uStar, CT, CH
  end
  return SurfaceValues, SurfaceData
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

function (::MOSurface)(Phys,Param,uPos,vPos,wPos)
  function SurfaceData(U,p,dXdxI,nS)
    FT = eltype(U)
    uStar = uStarCoefficientGPU(U[uPos],U[vPos],U[wPos],dXdxI,nS)
    CT = FT(Param.CE)
    CH = FT(Param.CH)
    return uStar, CT, CH
  end
  return SurfaceData
end  
