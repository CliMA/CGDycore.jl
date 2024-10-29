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

function (profile::SimpleKoefficient)(Param,Phys)
  @inline function Eddy(U,uStar,p,dz,LenScale)
    K = Param.CE * uStar * dz / 2
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

Base.@kwdef struct FriersonKoefficient{FT} <: EddyKoefficient 
  f_b = FT(0.1)
  z0 = FT(3.21e-5)
  Ri_C = FT(1)
end

function (profile::FriersonKoefficient)(Param,Phys,TkePos,RhoPos)
  @inline function Eddy(U,CM,uStar,RiBSurf,hBL,ThS,p,z,dz,LenScale)
    FT = eltype(U)
    norm_uh = U[uPos]^2 + U[vPos]^2
    Th = U[ThPos] / U[RhoPos]
    Ri = Phys.Grav * z[iz] * (Th - ThS) / ThS / norm_uh
    zbl = f_b * hBL
    if z < zbl
      if RiBSurf < FT(0)
        K = Phys.Karm * uStar * sqrt(CM) * z  
      else
        K = Phys.Karm * uStar * sqrt(CM) * z /
          (FT(1) + Ri / Ri_C * log(z / zM) / (FT(1) - Ri / Ri_C))
      end    
    elseif z < hBL  
      if RiBSurf < FT(0)
        Kb = Phys.Karm * uStar * sqrt(CM) * zbl  
      else
        Kb = Phys.Karm * uStar * sqrt(CM) * zbl /
          (FT(1) + Ri / Ri_C * log(zbl / zM) / (FT(1) - Ri / Ri_C))
      end    
      K = Kb * z / zbl * (FT(1) - (z - zbl) / (h - zbl))^2
    end    
    return K 
  end
  return Eddy
end  
