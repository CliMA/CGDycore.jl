################################################################################
# Two-point means for split-form fluxes.
#
# For nearly equal arguments the closed formulas cancel badly, in particular
# in Float32. With f = (x - y) / (x + y) all means are evaluated by a series
# in f^2 when f^2 < eps(FT)^(1/4); the truncation error of the series is then
# below eps(FT) (threshold 1.2e-4 for Float64, 1.9e-2 for Float32).
################################################################################

@inline series_threshold(::Type{FT}) where {FT} = sqrt(sqrt(eps(FT)))

# f^2 with f = (x - y) / (x + y), without cancellation in x - y
@inline mean_f2(x, y) = (x * (x - 2 * y) + y * y) / (x * (x + 2 * y) + y * y)

# Logarithmic mean (x - y) / (log(x) - log(y)), Ismail and Roe (2009)
@inline function ln_mean(x, y)
  FT = typeof(x)
  f2 = mean_f2(x, y)
  if f2 < series_threshold(FT)
    return (x + y) / @evalpoly(f2, FT(2), FT(2 / 3), FT(2 / 5), FT(2 / 7))
  else
    return (y - x) / log(y / x)
  end
end

# Inverse logarithmic mean (log(x) - log(y)) / (x - y)
@inline function inv_ln_mean(x, y)
  FT = typeof(x)
  f2 = mean_f2(x, y)
  if f2 < series_threshold(FT)
    return @evalpoly(f2, FT(2), FT(2 / 3), FT(2 / 5), FT(2 / 7)) / (x + y)
  else
    return log(y / x) / (y - x)
  end
end

# Stolarsky mean (gamma - 1) / gamma * (y^gamma - x^gamma) / (y^(gamma-1) - x^(gamma-1)),
# Winters et al. (2020)
@inline function stolarsky_mean(x, y, gamma)
  FT = typeof(x)
  f2 = mean_f2(x, y)
  if f2 < series_threshold(FT)
    c1 = FT(1 / 3) * (gamma - FT(2))
    c2 = FT(-1 / 15) * (gamma + FT(1)) * (gamma - FT(3)) * c1
    c3 = FT(-1 / 21) * (FT(2) * gamma * (gamma - FT(2)) - FT(9)) * c2
    return FT(0.5) * (x + y) * @evalpoly(f2, FT(1), c1, c2, c3)
  else
    return (gamma - FT(1)) / gamma * (y^gamma - x^gamma) /
      (y^(gamma - FT(1)) - x^(gamma - FT(1)))
  end
end

# Density mean used in the density flux and in the gravity term rho * [[phi]]
abstract type DensityMean end
struct ArithMean <: DensityMean end
struct LogMean <: DensityMean end       # well-balanced for T = const
struct GammaMean <: DensityMean end     # well-balanced for theta = const

@inline density_mean(::ArithMean, x, y, gamma) = typeof(x)(0.5) * (x + y)
@inline density_mean(::LogMean, x, y, gamma) = ln_mean(x, y)
@inline density_mean(::GammaMean, x, y, gamma) = stolarsky_mean(x, y, gamma)
