module AbstractSphericalHEVIDG

using LinearAlgebra
using StaticArrays

export AbstractGeometryMap, CubedSphereRadialMap
export MetricTerms, CartesianState
export extract_cartesian_flux, transform_to_contravariant!
export backtransform_increment, solve_hevi_column!

# -------------------------------------------------------------------
# 1. STATE & ABSTRACT GEOMETRY INTERFACE
# -------------------------------------------------------------------

"""
    CartesianState{T}
State vector q = [ρ, ρu, ρv, ρw, ρθ]^T in Earth-centered 3D Cartesian coordinates.
"""
struct CartesianState{T} <: FieldVector{5, T}
    ρ::T
    ρu::T
    ρv::T
    ρw::T
    ρθ::T
end

"""
    MetricTerms{T}
Holds local metric terms at a DG quadrature point, abstracting the mesh.
- J_inv: ∂ξ_i / ∂x_j (3x3 Matrix mapping Cartesian to Computational derivatives)
- detJ: |J| (Volume element scale factor)
"""
struct MetricTerms{T}
    J_inv::SMatrix{3, 3, T, 9} 
    detJ::T                     
end

abstract type AbstractGeometryMap end

"""
    CubedSphereRadialMap
Concrete implementation of a Cubed-Sphere grid extruded radially 
with terrain-following topography h(x,y,z).
"""
struct CubedSphereRadialMap <: AbstractGeometryMap
    R_earth::Float64
    z_top::Float64
end

# -------------------------------------------------------------------
# 2. KERNEL 1: PURE CARTESIAN FLUX EXTRACTOR
# -------------------------------------------------------------------

"""
    extract_cartesian_flux(q::CartesianState, p::Real)
Computes pure 3D Cartesian flux tensor F_cart ∈ R^{5 × 3}.
"""
@inline function extract_cartesian_flux(q::CartesianState{T}, p::T) where {T}
    inv_ρ = one(T) / q.ρ
    u = SVector{3, T}(q.ρu * inv_ρ, q.ρv * inv_ρ, q.ρw * inv_ρ)
    
    F_mass  = q.ρ * u
    F_mom_x = SVector{3, T}(q.ρu * u[1] + p, q.ρu * u[2],     q.ρu * u[3])
    F_mom_y = SVector{3, T}(q.ρv * u[1],     q.ρv * u[2] + p, q.ρv * u[3])
    F_mom_z = SVector{3, T}(q.ρw * u[1],     q.ρw * u[2],     q.ρw * u[3] + p)
    F_θ     = q.ρθ * u

    return (F_mass, F_mom_x, F_mom_y, F_mom_z, F_θ)
end

# -------------------------------------------------------------------
# 3. METRIC SUBROUTINE: TRANSFORM CARTESIAN TO CONTRAVARIANT FLUXES
# -------------------------------------------------------------------

"""
    transform_to_contravariant(F_cart, metric::MetricTerms)

Maps standard Cartesian fluxes F_cart to computational contravariant 
fluxes F_ξ = [F_ξ1, F_ξ2, F_ξ3] using the geometric inverse Jacobian J_inv.
"""
@inline function transform_to_contravariant(
    F_cart::NTuple{5, SVector{3, T}}, 
    metric::MetricTerms{T}
) where {T}
    # Extract gradient vectors ∇_x(ξ_i) from inverse Jacobian rows
    ∇ξ1 = metric.J_inv[1, :]
    ∇ξ2 = metric.J_inv[2, :]
    ∇ξ3 = metric.J_inv[3, :]

    # Project Cartesian flux tensor along metric directional gradients
    F_ξ1 = SVector{5, T}(dot(F_cart[1], ∇ξ1), dot(F_cart[2], ∇ξ1), dot(F_cart[3], ∇ξ1), dot(F_cart[4], ∇ξ1), dot(F_cart[5], ∇ξ1))
    F_ξ2 = SVector{5, T}(dot(F_cart[1], ∇ξ2), dot(F_cart[2], ∇ξ2), dot(F_cart[3], ∇ξ2), dot(F_cart[4], ∇ξ2), dot(F_cart[5], ∇ξ2))
    F_ξ3 = SVector{5, T}(dot(F_cart[1], ∇ξ3), dot(F_cart[2], ∇ξ3), dot(F_cart[3], ∇ξ3), dot(F_cart[4], ∇ξ3), dot(F_cart[5], ∇ξ3))

    # Scale by Volume Jacobian |J| for conservation statement
    return (F_ξ1 .* metric.detJ, F_ξ2 .* metric.detJ, F_ξ3 .* metric.detJ)
end

# -------------------------------------------------------------------
# 4. KERNEL 2: INCREMENT BACKTRANSFORMATION
# -------------------------------------------------------------------

"""
    backtransform_increment(dq_comp, metric::MetricTerms)

Converts stage update increment dq_comp computed in the reference element 
back into the geocentric physical Cartesian state increment dq_cart.
"""
@inline function backtransform_increment(
    dq_comp::SVector{5, T}, 
    metric::MetricTerms{T}
) where {T}
    inv_detJ = one(T) / metric.detJ
    
    # Primitive mapping back to geocentric Cartesian Frame
    return CartesianState{T}(
        dq_comp[1] * inv_detJ,
        dq_comp[2] * inv_detJ,
        dq_comp[3] * inv_detJ,
        dq_comp[4] * inv_detJ,
        dq_comp[5] * inv_detJ
    )
end

# -------------------------------------------------------------------
# 5. HEVI COLUMN DRIVER
# -------------------------------------------------------------------

"""
    solve_hevi_column!(q_col, metrics_col, dt)

Solves a vertical column without explicit knowledge of whether the grid 
is spherical, cubed-sphere, or flat.
"""
function solve_hevi_column!(
    q_col::Vector{CartesianState{T}}, 
    metrics_col::Vector{MetricTerms{T}},
    dt::T
) where {T}
    N_vert = length(q_col)
    
    dq_comp_col = Vector{SVector{5, T}}(undef, N_vert)
    
    for k in 1:N_vert
        # 1. Equation of state (in Cartesian Variables)
        p_k = 100000.0 * (q_col[k].ρθ * 287.0 / 100000.0)^1.4 
        
        # 2. Kernel 1: Extract Cartesian Flux
        F_cart = extract_cartesian_flux(q_col[k], p_k)
        
        # 3. Subroutine: Map using hidden metric layer
        F_ξ1, F_ξ2, F_ξ3 = transform_to_contravariant(F_cart, metrics_col[k])
        
        # 4. HEVI Column Solve Step (Vertical Implicit Operator acts on F_ξ3)
        # Placeholder for 1D block-tridiagonal solve along radial index k
        dq_comp_col[k] = -dt * F_ξ3 
    end

    # 5. Kernel 2: Backtransform & Update Cartesian State
    for k in 1:N_vert
        dq_cart = backtransform_increment(dq_comp_col[k], metrics_col[k])
        q_col[k] = q_col[k] + dq_cart
    end

    return nothing
end

end # module AbstractSphericalHEVIDG
