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
