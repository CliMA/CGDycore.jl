abstract type AbstractTransportType end
struct MomentumTransport <: AbstractTransportType end
struct HeatTransport <: AbstractTransportType end

abstract type UniversalFunction end

Base.@kwdef struct Businger{FT} <: UniversalFunction 
  am::FT = 4.7
  ah::FT = 4.7
  Pr::FT = 0.74
end

# Nishizawa2018 Eq. A7
f_momentum(uf::Businger, zeta) = sqrt(sqrt(1 - 15 * zeta))

# Nishizawa2018 Eq. A8
f_heat(uf::Businger, zeta) = sqrt(1 - 9 * zeta)

function psi(uf::Businger, zeta, ::MomentumTransport)
    FT = eltype(zeta)
    if zeta < 0
        # Businger1971 Eq. A3 (zeta < 0)
        f_m = f_momentum(uf, zeta)
        log_term = log((1 + f_m)^2 * (1 + f_m^2) / 8)
        return log_term - 2 * atan(f_m) + FT(pi) / 2
    else
        # Businger1971 Eq. A3 (zeta >= 0)
        return -uf.am * zeta
    end
end

function psi(uf::Businger, zeta, tt::HeatTransport)
    if zeta < 0
        # Businger1971 Eq. A4 (zeta < 0)
        f_h = f_heat(uf, zeta)
        return 2 * log((1 + f_h) / 2)
    else
        # Businger1971 Eq. A4 (zeta >= 0)
        return -uf.ah * zeta / uf.Pr
    end
end

