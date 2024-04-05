abstract type TurbulenceModel end

Base.@kwdef struct TKEModel <: TurbulenceModel end

function (profile::TKEModel)(Param,Phys,RhoPos,uPos,vPos,ThPos,TkePos)
  @inline function Model(UB,UT,dzB,dzT)
    FT = eltype(UB)
    dzF = FT(0.5) * (dzB + dzT)
    dudzF = (UT[uPos] - UB[uPos]) / dzF
    dvdzF = (UT[vPos] - UB[vPos]) / dzF
    ThT = UT[ThPos] / UT[RhoPos]
    ThB = UB[ThPos] / UB[RhoPos]
    N2 = Phys.Grav * FT(0.5) * (ThT- ThB) / dzF / (ThT + ThB)
    S = sqrt(dudzF * dudzF + dvdzF * dvdzF)
    LenScale = FT(1000.0)
    TkeF = FT(0.5) * (UT[TkePos] / UT[RhoPos] + UB[TkePos] / UB[RhoPos])
    TkeFAbs = max(TkeF, FT(1.e-8))
    RhoF = FT(0.5) * (UT[RhoPos] + UB[RhoPos])
    sqrTkeFAbs = sqrt(TkeFAbs)
    DiffKoeff = Phys.Cd * sqrTkeFAbs / LenScale
    S = DiffKoeff * S + FT(1.e-8)
    P = S * (FT(1) - min(Phys.PrTke * DiffKoeff * N2 / S, FT(0.75)))
    if N2 > FT(0)
      Len = min(LenScale,FT(0.76) * sqrt(TkeF / N2))
    else  
      Len = LenScale
    end
    Cd = FT(1.9) * Phys.Cd +(FT(0.93) - FT(1.9) * Phys.Cd) * Len / LenScale
    Diss = Phys.Cd * TkeF * sqrTkeFAbs / LenScale
    F = P - Diss
    KV = DiffKoeff * RhoF
    return F, KV
  end
  return Model
end  
