#=
DCMIP-2025 mountain-generated mesoscale test (Andrews et al., EGUsphere-2026-2293)
Initial conditions on a height (z) coordinate for
  Case A: gap flow         (mountain chain with a gap, h0 = 1500 m)
  Case B: vortex shedding  (isolated Gaussian mountain, h0 = 2000 m)

Returned prognostic state at (λ, φ, z):
    ρ       density                     [kg m^-3]
    θ       potential temperature       [K]
    u, v, w local zonal / meridional / vertical wind [m s^-1]
  (+ optional geocentric Cartesian components Ux, Uy, Uz, see `wind_cartesian`)

Equations used (paper numbering):
    u  = u0 cos φ, v = 0                                   (1a,1b)
    T  = T0                                                 (1c)
    ps = psp exp( -a N² u0/(2 g² κ) (u0/a + 2Ω)(sin²φ - 1)
                  - N²/(g² κ) Φs )                          (1e)
    p  = ps exp( -g (z - zs) / (Rd T0) )                    (2)
    ρ  = p / (Rd T0)                                        (equation of state)
    θ  = T0 (p0/p)^κ                                        (p0 = 1e5 Pa)
=#

# ----------------------------------------------------------------------------
# Constants (Table 1)
# ----------------------------------------------------------------------------
Base.@kwdef struct Params{FT<:AbstractFloat}
    X::FT   = 20            # small-Earth factor
    u0::FT  = 10            # max zonal wind at equator [m/s]
    T0::FT  = 288           # isothermal temperature [K]
    g::FT   = 9.80616       # gravity [m/s²]
    a::FT   = 6.371229e6/20 # small-Earth radius [m]
    cp::FT  = 1004.64
    Rd::FT  = 287.04
    psp::FT = 1e5           # surface pressure at the pole where Φs≈0 [Pa]
    p0::FT  = 1e5           # reference pressure for θ [Pa]
    Ω::FT   = 20*7.2921e-5  # X·Ω_reg; set 0 for the non-rotating runs
end

κ(P)  = P.Rd/P.cp
N2(P) = P.g^2/(P.cp*P.T0)          # Brunt–Väisälä frequency², N ≈ 0.0182 1/s (Eq. 3)

"Parameters with or without rotation (and an optional different X)."
function make_params(; rotation::Bool=true, X::Real=20, FT=Float64)
    return Params{FT}(X=X, a=6.371229e6/X, Ω = rotation ? X*7.2921e-5 : 0.0)
end

# ----------------------------------------------------------------------------
# Orography
# ----------------------------------------------------------------------------
@inline wrap_dλ(Δ) = mod(Δ + π, 2π) - π     # wrap longitude difference to [-π, π)

# ---- Case A: mountain chain with gap (Eqs. 18-20) ---------------------------
Base.@kwdef struct GapOrography{FT<:AbstractFloat}
    h0::FT   = 1500
    λc::FT   = π
    φc::FT   = 0
    e1::Int  = 10
    e2::Int  = 10
    e3::Int  = 10
    d1::FT
    d2::FT
    d3::FT
end

function GapOrography(P::Params{FT}; h0=1500, xlon=800e3/P.X, xlat=6000e3/P.X,
                      xgap=1000e3/P.X, e=10) where {FT}
    # length factors, Eq. (20), with h0/zs = 10
    di(x) = x/(2P.a) * log(FT(10))^(-1/FT(e))
    return GapOrography{FT}(h0=FT(h0), e1=e, e2=e, e3=e,
                            d1=di(xlon), d2=di(xlat), d3=di(xgap))
end

function surface_height(o::GapOrography, λ, φ)
    dλ = wrap_dλ(λ - o.λc)
    f  = o.h0*exp(-(dλ/o.d1)^o.e1 - ((φ - o.φc)/o.d2)^o.e2)
    gp = 1 - exp(-((φ - o.φc)/o.d3)^o.e3)
    return f*gp
end

# ∂zs/∂λ (Eq. A4; the gap factor g(φ) is already contained in zs)
function dzs_dλ(o::GapOrography, λ, φ)
    dλ = wrap_dλ(λ - o.λc)
    return -o.e1/o.d1 * (dλ/o.d1)^(o.e1-1) * surface_height(o, λ, φ)
end

# ---- Case B: isolated Gaussian mountain (Eqs. 23-24) -----------------------
Base.@kwdef struct VortexOrography{FT<:AbstractFloat}
    h0::FT = 2000
    λc::FT = π
    φc::FT = π/9
    d::FT  = 250e3/20          # Gaussian half-width, X⁻¹·250 km = 12.5 km
    a::FT  = 6.371229e6/20
end

VortexOrography(P::Params{FT}; h0=2000, d=250e3/P.X) where {FT} =
    VortexOrography{FT}(h0=FT(h0), d=FT(d), a=P.a)

# great-circle distance (haversine form, equivalent to Eq. 24 but accurate near r = 0)
function great_circle(o::VortexOrography, λ, φ)
    s = sin((φ - o.φc)/2)^2 + cos(φ)*cos(o.φc)*sin((λ - o.λc)/2)^2
    return 2o.a*asin(min(one(s), sqrt(s)))
end

function surface_height(o::VortexOrography, λ, φ)
    r = great_circle(o, λ, φ)
    return o.h0*exp(-(r/o.d)^2)
end

# ∂zs/∂λ (Eqs. A6-A7)
function dzs_dλ(o::VortexOrography, λ, φ)
    r = great_circle(o, λ, φ)
    r == 0 && return zero(r)
    drdλ = o.a*sin(λ - o.λc)*cos(o.φc)*cos(φ) / sin(r/o.a)
    return -2r/o.d^2 * drdλ * surface_height(o, λ, φ)
end

# ----------------------------------------------------------------------------
# Thermodynamic fields
# ----------------------------------------------------------------------------
"Balanced surface pressure, Eq. (1e)."
function surface_pressure(P::Params, orog, λ, φ)
    Φs = P.g*surface_height(orog, λ, φ)
    return P.psp*exp( -(P.a*N2(P)*P.u0)/(2P.g^2*κ(P)) * (P.u0/P.a + 2P.Ω)*(sin(φ)^2 - 1)
                      - N2(P)/(P.g^2*κ(P))*Φs )
end

"Pressure on a height level z (above sea level), Eq. (2)."
function pressure(P::Params, orog, λ, φ, z)
    zs = surface_height(orog, λ, φ)
    return surface_pressure(P, orog, λ, φ)*exp(-P.g*(z - zs)/(P.Rd*P.T0))
end

# ----------------------------------------------------------------------------
# Optional nonzero initial vertical velocity for terrain-following coordinates
# (Appendix A).  Blending functions A(z): z = z_base + A(z_base)·zs  (Eq. B5)
# ----------------------------------------------------------------------------
blend_mpas(z; zT=20007.5)   = cos(π/2*z/zT)^6      # Eq. (B6)
blend_gungho(z; zT=20007.5) = (zT - z)/zT            # Eq. (B7)

# ----------------------------------------------------------------------------
# Full initial state at one point
# ----------------------------------------------------------------------------
"""
    initial_state(P, orog, λ, φ, z; w_mode=:zero, blend=blend_mpas)

Return `(ρ, θ, u, v, w)` at longitude λ, latitude φ [rad] and height z [m].

* `w_mode = :zero`   → w = 0 (the choice used for MPAS/GungHo in the paper).
* `w_mode = :follow` → w = u/(a cos φ) · ∂zs/∂λ · A(z)  (Eq. A2/A3), i.e. flow follows the
                       sloped coordinate surfaces of a hybrid-z grid.  `z` should then be the
                       *base* (undeformed) height used in A(z).
"""
function initial_state(P::Params, orog, λ, φ, z; w_mode::Symbol=:zero, blend=blend_mpas)
    p   = pressure(P, orog, λ, φ, z)
    ρ   = p/(P.Rd*P.T0)
    θ   = P.T0*(P.p0/p)^κ(P)
    u   = P.u0*cos(φ)
    v   = zero(u)
    w   = zero(u)
    if w_mode === :follow
        c = cos(φ)
        w = abs(c) < 1e-12 ? zero(u) : u/(P.a*c)*dzs_dλ(orog, λ, φ)*blend(z)
    end
    return ρ, θ, u, v, w
end

"""
    wind_cartesian(u, v, w, λ, φ) -> (Ux, Uy, Uz)

Convert local (east, north, up) wind to geocentric Cartesian components.
"""
function wind_cartesian(u, v, w, λ, φ)
    sλ, cλ, sφ, cφ = sin(λ), cos(λ), sin(φ), cos(φ)
    Ux = -sλ*u - sφ*cλ*v + cφ*cλ*w
    Uy =  cλ*u - sφ*sλ*v + cφ*sλ*w
    Uz =           cφ*v   + sφ*w
    return Ux, Uy, Uz
end

# ----------------------------------------------------------------------------
# Base vertical grid (Eqs. 9-10): 57 layers, top ≈ 20007.5 m
# ----------------------------------------------------------------------------
function base_vertical_grid(; zL=1000.0, zU=6000.0, zT=20000.0,
                              dzmin=100.0, dzmax=500.0, γ=1.01679)
    zi = [0.0]; dz = dzmin
    while zi[end] < zT - 1e-6
        zb = zi[end]                              # lower interface z_{k-1/2}
        if zb < zL - 1e-9
            dz = dzmin
        elseif zb <= zU + 1e-9
            dz = min(dz^γ, dzmax)
        else
            dz = dzmax
        end
        push!(zi, zb + dz)
    end
    zm = 0.5 .* (zi[1:end-1] .+ zi[2:end])        # mid-levels, Eq. (B8)
    return zi, zm
end

# ----------------------------------------------------------------------------
# Fill 3-D arrays on a regular lon-lat grid and pure z levels
# ----------------------------------------------------------------------------
"""
    fill_initial_conditions(P, orog, λs, φs, zs_levels; kwargs...)

Returns a NamedTuple of arrays `(ρ, θ, u, v, w, Ux, Uy, Uz, zsurf)` sized (nλ, nφ, nz).
Levels below the surface (z < zs) are filled with the analytic extension and are flagged
in the `inside` Bool array so that you can mask them.
"""
function fill_initial_conditions(P::Params, orog, λs, φs, zlev; kwargs...)
    nλ, nφ, nz = length(λs), length(φs), length(zlev)
    FT = typeof(P.a)
    ρ  = zeros(FT, nλ, nφ, nz); θ = similar(ρ)
    u  = similar(ρ); v = similar(ρ); w = similar(ρ)
    Ux = similar(ρ); Uy = similar(ρ); Uz = similar(ρ)
    inside = falses(nλ, nφ, nz)
    zsurf  = [surface_height(orog, λ, φ) for λ in λs, φ in φs]
    for k in 1:nz, j in 1:nφ, i in 1:nλ
        ρ[i,j,k], θ[i,j,k], u[i,j,k], v[i,j,k], w[i,j,k] =
            initial_state(P, orog, λs[i], φs[j], zlev[k]; kwargs...)
        Ux[i,j,k], Uy[i,j,k], Uz[i,j,k] =
            wind_cartesian(u[i,j,k], v[i,j,k], w[i,j,k], λs[i], φs[j])
        inside[i,j,k] = zlev[k] < zsurf[i,j]
    end
    return (; ρ, θ, u, v, w, Ux, Uy, Uz, zsurf, inside)
end

# ----------------------------------------------------------------------------
# Example
# ----------------------------------------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    P  = make_params(rotation=true)            # make_params(rotation=false) for Ω = 0
    zi, zm = base_vertical_grid()
    println("levels = ", length(zm), ",  model top = ", zi[end], " m")

    # 0.5° lon-lat grid (720 x 361)
    λs = range(0, 2π, length=721)[1:end-1]
    φs = range(-π/2, π/2, length=361)

    for (name, orog) in (("gap flow", GapOrography(P)),
                         ("vortex shedding", VortexOrography(P)))
        ic = fill_initial_conditions(P, orog, λs, φs, zm)
        println(name, ": max zs = ", maximum(ic.zsurf), " m, ",
                "ρ range = ", extrema(ic.ρ), ", θ range = ", extrema(ic.θ))
    end
end
