abstract type TurbulenceModel end

"""
    TKEModel{FT}

One-and-a-half order closure with a prognostic turbulent kinetic energy e (stored as rho*e).

  K    = Cm * l * sqrt(e)                      eddy diffusivity (momentum, heat and e)
  P    = K * (S^2 - N^2)                        shear and buoyancy production
  eps  = Ceps * e^(3/2) / l                     dissipation
  l    = kappa z / (1 + kappa z / LMax)         Blackadar, limited by Crit*sqrt(e/N^2) for N^2 > 0

With Ceps = Cm^3 the neutral surface layer gives K = kappa u* z and e = u*^2 / Cm^2
(Mellor-Yamada level 2.5: Cm = 0.554, Ceps = 0.17, Crit = 0.75).
"""
Base.@kwdef struct TKEModel{FT} <: TurbulenceModel
  Cm::FT = 0.554
  Ceps::FT = 0.554^3
  LMax::FT = 100.0
  LMin::FT = 0.1
  Crit::FT = 0.75
  EMin::FT = 1.0e-4
  Karman::FT = 0.4
end

@inline function MixingLength(Turb::TKEModel,z,e,N2)
  FT = eltype(e)
  kz = Turb.Karman * max(z, FT(0))
  l = kz / (FT(1) + kz / Turb.LMax)
  if N2 > FT(0)
    l = min(l, Turb.Crit * sqrt(e / N2))
  end
  return max(l, Turb.LMin)
end

@inline EddyDiffusivity(Turb::TKEModel,l,e) = Turb.Cm * l * sqrt(e)

@inline Dissipation(Turb::TKEModel,l,e) = Turb.Ceps * e * sqrt(e) / l

# derivative of the sinks (dissipation, buoyancy destruction) with respect to e;
# it goes on the diagonal of the implicit column solve
@inline function DissipationDerivative(Turb::TKEModel,l,e,N2)
  FT = eltype(e)
  D = FT(1.5) * Turb.Ceps * sqrt(e) / l
  if N2 > FT(0)
    D += FT(0.5) * Turb.Cm * l * N2 / sqrt(e)
  end
  return D
end
