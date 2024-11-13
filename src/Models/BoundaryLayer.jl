"""
Notation:

CM Drag coefficient momentum
CT Drag coefficient temperature
CH Drag coefficient water vapor
uStar
RiBSurf surface bulk Richardson number
hBL Boundary layer height
"""
abstract type BoundaryLayerHeight end

Base.@kwdef struct BLHRichardson{FT} <: BoundaryLayerHeight 
  Ri_C::FT
end

function (::BLHRichardson)(Phys,uPos,vPos,ThPos)
  @inline function BoundaryLayerHeight(hBL,U,z)
    Theta_1 = U[1,ThPos]
    hBL = z[end]
    for iz = 2 : length(z)
      norm_uh = U[iz,uPos]^2 + U[iz,vPos]^2
      Ri = Phys.Grav * z[iz] * (U[iz,ThPos] - Theta_1) / Theta_1 / norm_uh  
      if Ri > Ri_C
        hBL = z[i]
        break
      end
    end  
    return hBL
  end
  return BoundaryLayerHeight
end

abstract type EddyKoefficient end

Base.@kwdef struct SimpleKoefficient <: EddyKoefficient end

function (profile::SimpleKoefficient)(Param,Phys,uStar)
  @inline function Eddy(U,Surf,p,dz,LenScale)
    K = Param.CE * Surf[uStar] * dz / 2
    if p < Param.p_pbl
      dpR = (Param.p_pbl - p) / Param.p_strato
      K = K * exp(-dpR * dpR)
    end
    return K * U[1]
  end
  return Eddy
end

Base.@kwdef struct TkeKoefficient <: EddyKoefficient end

function (profile::TkeKoefficient)(Param,Phys,TkePos,RhoPos)
  @inline function Eddy(U,uStar,p,dz,LenScale)
    sqrTke = sqrt(max(abs(U[TkePos] / U[RhoPos]),eltype(U)(1.e-8)))
    K = Phys.Cd * U[RhoPos] * sqrTke / LenScale
    return K 
  end
  return Eddy
end  

Base.@kwdef struct FriersonKoefficient <: EddyKoefficient end

function (profile::FriersonKoefficient)(Param,Phys,RhoPos,uPos,vPos,ThPos,TS,hBL,uStar,CM,RiBSurf)
  @inline function Eddy(U,Surf,p,dz)
    FT = eltype(U)
    norm_uh = U[uPos]^2 + U[vPos]^2
    Th = U[ThPos] / U[RhoPos]
    Ri = Phys.Grav * z[iz] * (Th - Surf[TS]) / Surf[TS] / norm_uh
    zbl = Param.f_b * Surf[hBL]
    if z < zbl
      if Surf[RiBSurf] < FT(0)
        K = Phys.Karm * Surf[uStar] * sqrt(Surf[CM]) * z  
      else
        K = Phys.Karm * Surf[uStar] * sqrt(Surf[CM]) * z /
          (FT(1) + Ri / Param.Ri_C * log(z / Param.zM) / (FT(1) - Ri / Param.Ri_C))
      end    
    elseif z < hBL  
      if Surf[RiBSurf] < FT(0)
        Kb = Phys.Karm * Surf[uStar] * sqrt(Surf[CM]) * zbl  
      else
        Kb = Phys.Karm * Surf[uStar] * sqrt(Surf[CM]) * zbl /
          (FT(1) + Ri / Param.Ri_C * log(zbl / zM) / (FT(1) - Ri / Param.Ri_C))
      end    
      K = Kb * z / zbl * (FT(1) - (z - zbl) / (h - zbl))^2
    end    
    return K * U[RhoPos]
  end
  return Eddy
end  
