abstract type TurbulenceModel end

Base.@kwdef struct TKEModel <: TurbulenceModel end

#=
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
    LenScale = dzF
    RhoF = FT(0.5) * (UT[RhoPos] + UB[RhoPos])
    TkeF = UB[TkePos] / RhoF
    TkeFAbs = max(TkeF, FT(1.e-8))
    sqrTkeFAbs = sqrt(TkeFAbs)
    DiffKoeff = Phys.Cd * sqrTkeFAbs / LenScale
    S = DiffKoeff * S * S + FT(1.e-8)
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
=#


function (profile::TKEModel)(Param,Phys,RhoPos,uPos,vPos,ThPos,TkePos)
  @inline function Model(UB,UT,dzB,dzT)

#   Constants
    SigT = 0.74e0
    fTke = 0.1

    FT = eltype(UB)
    dzF = FT(0.5) * (dzB + dzT)

    dudzF = (UT[uPos] - UB[uPos]) / dzF
    dvdzF = (UT[vPos] - UB[vPos]) / dzF
    S = sqrt(dudzF * dudzF + dvdzF * dvdzF)

    ThT = UT[ThPos] / UT[RhoPos]
    ThB = UB[ThPos] / UB[RhoPos]
    N2 = Phys.Grav * FT(2.0) * (ThT- ThB) / dzF / (ThT + ThB)

    LenScale = dzF
#   Diffusion Koefficient
    RhoF = FT(0.5) * (UT[RhoPos] + UB[RhoPos])
    TkeF = FT(0.5) * (UB[TkePos] / UB[RhoPos] + UT[TkePos] / UT[RhoPos])
    TkeFAbs = max(TkeF, FT(1.e-8))
    sqrTkeFAbs = sqrt(TkeFAbs)
    DiffKoeff = RhoF * max(Phys.Cd * sqrTkeFAbs * LenScale, 1.e-2)

#   Richardson-number and production terms
    Rich = N2 / (S * S + 1.e-3)
    Rich = max(Rich,-4.0/3.0)
    Ptke = max(1.0 - DiffKoeff*Rich/(SigT*DiffKoeff),fTke)

#   Local Length Scale
    if N2 > 0.0
      Len = min(LenScale, 0.76 * sqrTkeFAbs / N2)
    else
      Len = LenScale
    end

    CDis = FT(1.9) * Phys.Cd +(FT(0.93) - FT(1.9) * Phys.Cd) * Len / LenScale
    Diss = Phys.Cd * TkeF * sqrTkeFAbs / LenScale
    F = DiffKoeff * S * S * Ptke - CDis * TkeF * sqrTkeFAbs * RhoF / Len

    return F, DiffKoeff
  end
  return Model
end

