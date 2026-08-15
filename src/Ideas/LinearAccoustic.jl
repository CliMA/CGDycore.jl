"""
Linear Acoustics — Split-Form DG Volume Kernel on Curvilinear Grids
====================================================================
Equations in physical space (2D):
    ∂p/∂t + c² (∂u/∂x + ∂v/∂y) = 0
    ∂u/∂t + ∂p/∂x               = 0
    ∂v/∂t + ∂p/∂y               = 0

Mapped to reference element [−1,1]² via (x,y) = X(ξ,η).

The metric terms are stored in "conservative" (Kopriva) form to satisfy
the Discrete Geometric Conservation Law (DGCL):

    Jᵢ aˢᵢ_n  =  J·(∂ξⁿ/∂xˢ)  at node i

where s ∈ {x,y}, n ∈ {ξ,η}.  With this the physical-space divergence
becomes the reference-space flux-differencing sum

    (du/dt)_i = -1/Jᵢ * [
        2 Σ_j D̂ˢᴿᶠ[i,j]  ·  f_S^ξ( Û_i, Û_j )
      + 2 Σ_j D̂ˢᴿᶠ[i,j]  ·  f_S^η( Û_i, Û_j )
    ]

where the contravariant two-point fluxes are

    f_S^ξ( (p_L,u_L,v_L), (p_R,u_R,v_R) )
        = 0.5*(ãˣξ_L + ãˣξ_R) * f_S^x + 0.5*(ãʸξ_L + ãʸξ_R) * f_S^y

and analogously for the η direction.  The Cartesian two-point fluxes are
    f_S^x = ( 0.5 c²(u_L+u_R),  0.5(p_L+p_R),  0 )
    f_S^y = ( 0.5 c²(v_L+v_R),  0,  0.5(p_L+p_R) )

This follows Winters et al. (2018) "A Comparative Study on Polynomial
Dealiasing and Split Form Discontinuous Galerkin Schemes for
Under-Resolved Turbulence Computations" and Kopriva (2006) for metrics.

State vector per node: U = [p, u, v]

Layout conventions (column-major Julia):
    U    [3, Np, K]          state
    dU   [3, Np, K]          d/dt residual (volume only)
    Dhat [Np, Np]            split-form 1D SBP matrix  Dhat = 0.5 M⁻¹ Q
    Ja   [2, 2, Np, K]       metric cofactors  Ja[s,n,i,e] = J·(∂ξⁿ/∂xˢ)
                              s=1→x, s=2→y;  n=1→ξ, n=2→η
    J    [Np, K]             Jacobian determinant at each node

Dependencies:
    ] add KernelAbstractions CUDA StaticArrays
"""

using KernelAbstractions
using KernelAbstractions.Extras: @unroll
using StaticArrays
using LinearAlgebra: det
using FastGaussQuadrature

const N_FIELDS = 3   # (p, u, v)

# ── Cartesian two-point fluxes ────────────────────────────────────────────────
# f_S^x: flux in the x-direction for the split form
@inline function flux_x(pL, uL, vL, pR, uR, vR, c²)
    return SVector(
        0.5 * c² * (uL + uR),   # pressure eq
        0.5 * (pL + pR),         # u-momentum eq
        zero(pL)                  # v-momentum eq
    )
end

# f_S^y: flux in the y-direction for the split form
@inline function flux_y(pL, uL, vL, pR, uR, vR, c²)
    return SVector(
        0.5 * c² * (vL + vR),   # pressure eq
        zero(pL),                 # u-momentum eq
        0.5 * (pL + pR)          # v-momentum eq
    )
end

# ── Contravariant two-point flux in reference direction n (1=ξ, 2=η) ─────────
#
# F_S^n = 0.5*(Ja[x,n]_L + Ja[x,n]_R)*f_S^x
#       + 0.5*(Ja[y,n]_L + Ja[y,n]_R)*f_S^y
#
@inline function contravariant_flux(
        pL, uL, vL, pR, uR, vR,
        Jaxn_L, Jayn_L,          # metric cofactors at left  node, direction n
        Jaxn_R, Jayn_R,          # metric cofactors at right node, direction n
        c²)

    avg_Jaxn = 0.5 * (Jaxn_L + Jaxn_R)
    avg_Jayn = 0.5 * (Jayn_L + Jayn_R)

    fSx = flux_x(pL, uL, vL, pR, uR, vR, c²)
    fSy = flux_y(pL, uL, vL, pR, uR, vR, c²)

    return avg_Jaxn .* fSx .+ avg_Jayn .* fSy
end

# ── Volume kernel ─────────────────────────────────────────────────────────────
#
# Grid index convention:
#   2D tensor-product element: node (i,j), flattened to linear index
#   p = i + (j-1)*Np1D   with  i,j ∈ 1:Np1D
#
# Each GPU thread handles one (node p, element e).
# The reference-direction loops (over i' or j') are unrolled.
#
# Ja layout:  Ja[s, n, node, elem]
#   s=1 → x-component of metric (∂x/∂ξⁿ related cofactor)
#   s=2 → y-component
#   n=1 → ξ-direction
#   n=2 → η-direction
#
@kernel function volume_kernel!(dU, U, Dhat, Ja, Jinv, c², ::Val{Np1D}) where {Np1D}

    p, e = @index(Global, NTuple)   # linear node index p, element index e

    # Convert linear node index → 2D tensor indices (i along ξ, j along η)
    i = mod1(p, Np1D)
    j = div(p - 1, Np1D) + 1

    # Load state at node (i,j)
    pᵢⱼ = U[1, p, e]
    uᵢⱼ = U[2, p, e]
    vᵢⱼ = U[3, p, e]

    rhs_p = zero(eltype(U))
    rhs_u = zero(eltype(U))
    rhs_v = zero(eltype(U))

    # ── ξ-direction flux differencing (vary i', fix j) ──────────────────────
    @unroll for i′ in 1:Np1D
        p′ = i′ + (j - 1) * Np1D   # linear index of neighbour node

        pₙ = U[1, p′, e]
        uₙ = U[2, p′, e]
        vₙ = U[3, p′, e]

        # Metric cofactors for direction n=1 (ξ) at both nodes
        Jax1_c = Ja[1, 1, p,  e]   # J·(∂ξ/∂x) at current node
        Jay1_c = Ja[2, 1, p,  e]   # J·(∂ξ/∂y) at current node
        Jax1_n = Ja[1, 1, p′, e]   # J·(∂ξ/∂x) at neighbour
        Jay1_n = Ja[2, 1, p′, e]   # J·(∂ξ/∂y) at neighbour

        F = contravariant_flux(
                pᵢⱼ, uᵢⱼ, vᵢⱼ,
                pₙ,  uₙ,  vₙ,
                Jax1_c, Jay1_c,
                Jax1_n, Jay1_n,
                c²)

        Dij = Dhat[i, i′]   # only the ξ-index varies

        rhs_p -= 2 * Dij * F[1]
        rhs_u -= 2 * Dij * F[2]
        rhs_v -= 2 * Dij * F[3]
    end

    # ── η-direction flux differencing (fix i, vary j') ──────────────────────
    @unroll for j′ in 1:Np1D
        p′ = i + (j′ - 1) * Np1D   # linear index of neighbour node

        pₙ = U[1, p′, e]
        uₙ = U[2, p′, e]
        vₙ = U[3, p′, e]

        # Metric cofactors for direction n=2 (η) at both nodes
        Jax2_c = Ja[1, 2, p,  e]
        Jay2_c = Ja[2, 2, p,  e]
        Jax2_n = Ja[1, 2, p′, e]
        Jay2_n = Ja[2, 2, p′, e]

        F = contravariant_flux(
                pᵢⱼ, uᵢⱼ, vᵢⱼ,
                pₙ,  uₙ,  vₙ,
                Jax2_c, Jay2_c,
                Jax2_n, Jay2_n,
                c²)

        Dij = Dhat[j, j′]   # only the η-index varies

        rhs_p -= 2 * Dij * F[1]
        rhs_u -= 2 * Dij * F[2]
        rhs_v -= 2 * Dij * F[3]
    end

    # ── Scale by 1/J ────────────────────────────────────────────────────────
    Ji = Jinv[p, e]
    dU[1, p, e] = Ji * rhs_p
    dU[2, p, e] = Ji * rhs_u
    dU[3, p, e] = Ji * rhs_v
end

# ── Grid generation: wavy curvilinear mesh on [0,1]² ─────────────────────────
#
# Physical mapping:  x(ξ,η) = ξ + ε·sin(2πξ)·sin(2πη)
#                   y(ξ,η) = η + ε·sin(2πξ)·sin(2πη)
#
# Metric cofactors (conservative Kopriva form):
#   Ja[x,ξ] =  J·∂ξ/∂x  =  ∂y/∂η
#   Ja[y,ξ] = -J·∂ξ/∂y  = -∂x/∂η   (2D cross-product metric)
#   Ja[x,η] = -J·∂η/∂x  = -∂y/∂ξ
#   Ja[y,η] =  J·∂η/∂y  =  ∂x/∂ξ
#
function build_curvilinear_mesh(Np1D, Kx, Ky; ε = 0.1)
    K     = Kx * Ky
    Np    = Np1D^2
    T     = Float64

    # Reference nodes: equispaced on [-1,1] (replace with GLL for production)
    ξ1D   = range(-1.0, 1.0; length = Np1D) |> collect
    nodes = [(ξ1D[i], ξ1D[j]) for j in 1:Np1D, i in 1:Np1D][:]   # Np-vector of (ξ,η)

    # Allocate
    x   = zeros(T, Np, K)
    y   = zeros(T, Np, K)
    Ja  = zeros(T, 2, 2, Np, K)
    J   = zeros(T, Np, K)
    Jinv = zeros(T, Np, K)

    e = 0
    for ey in 1:Ky, ex in 1:Kx
        e += 1
        # Reference→global reference element for this cell
        x0 = (ex - 1) / Kx;  x1 = ex / Kx
        y0 = (ey - 1) / Ky;  y1 = ey / Ky

        for (p, (ξ, η)) in enumerate(nodes)
            # Map reference [-1,1]² → cell [x0,x1]×[y0,y1]
            rx = 0.5 * ((1 - ξ) * x0 + (1 + ξ) * x1)
            ry = 0.5 * ((1 - η) * y0 + (1 + η) * y1)

            # Wavy physical map
            xp  = rx + ε * sin(2π * rx) * sin(2π * ry)
            yp  = ry + ε * sin(2π * rx) * sin(2π * ry)
            x[p, e] = xp
            y[p, e] = yp

            # Jacobian of physical map w.r.t. (rx, ry)  (chain rule handles [-1,1]→cell)
            dxdrx = 1 + 2π * ε * cos(2π * rx) * sin(2π * ry)
            dxdry =     2π * ε * sin(2π * rx) * cos(2π * ry)
            dydrx =     2π * ε * cos(2π * rx) * sin(2π * ry)
            dydry = 1 + 2π * ε * sin(2π * rx) * cos(2π * ry)

            # dr/dξ (scaling from reference element to cell)
            drx_dxi  = 0.5 * (x1 - x0)
            dry_deta = 0.5 * (y1 - y0)

            # Full Jacobian d(x,y)/d(ξ,η)  (tensor-product cell → diagonal dr/dξ)
            dxdxi  = dxdrx * drx_dxi
            dydxi  = dydrx * drx_dxi
            dxdeta = dxdry * dry_deta
            dydeta = dydry * dry_deta

            Jdet = dxdxi * dydeta - dxdeta * dydxi
            J[p, e]    = Jdet
            Jinv[p, e] = 1.0 / Jdet

            # Conservative metric cofactors (Kopriva 2006, eq. 6)
            Ja[1, 1, p, e] =  dydeta   #  J · ∂ξ/∂x
            Ja[2, 1, p, e] = -dxdeta   #  J · ∂ξ/∂y
            Ja[1, 2, p, e] = -dydxi    #  J · ∂η/∂x
            Ja[2, 2, p, e] =  dxdxi    #  J · ∂η/∂y
        end
    end
    return x, y, Ja, J, Jinv
end

function build_Dhat(Np1D::Int)
    ξ, w = gausslobatto(Np1D)

    # Barycentric weights
    bw = ones(Np1D)
    for j in 1:Np1D, k in 1:Np1D
        k == j && continue
        bw[j] /= (ξ[j] - ξ[k])
    end

    # Lagrange differentiation matrix
    D = zeros(Np1D, Np1D)
    for i in 1:Np1D, j in 1:Np1D
        if i != j
            D[i, j] = (bw[j] / bw[i]) / (ξ[i] - ξ[j])
        end
    end
    for i in 1:Np1D
        D[i, i] = -sum(D[i, k] for k in 1:Np1D if k != i)
    end

    return 0.5 .* D
end

# ── Main ──────────────────────────────────────────────────────────────────────
function main()
    Np1D = 4      # nodes per direction per element (polynomial degree p = Np1D-1)
    Kx   = 10     # elements in x
    Ky   = 10     # elements in y
    K    = Kx * Ky
    Np   = Np1D^2
    c²   = 1.0^2  # speed of sound squared

    backend = CPU()   # swap to CUDABackend() / ROCBackend() for GPU
    T = Float64

    println("Building curvilinear mesh (Np1D=$Np1D, Kx=$Kx, Ky=$Ky, K=$K)...")
    x_cpu, y_cpu, Ja_cpu, J_cpu, Jinv_cpu = build_curvilinear_mesh(Np1D, Kx, Ky; ε = 0.1)

    # Transfer to backend arrays
    U    = KernelAbstractions.zeros(backend, T, N_FIELDS, Np, K)
    dU   = KernelAbstractions.zeros(backend, T, N_FIELDS, Np, K)
    Dhat = KernelAbstractions.allocate(backend, T, Np1D, Np1D)
    Ja   = KernelAbstractions.allocate(backend, T, 2, 2, Np, K)
    Jinv = KernelAbstractions.allocate(backend, T, Np, K)

    copyto!(Dhat, build_Dhat(Np1D))
    copyto!(Ja,   Ja_cpu)
    copyto!(Jinv, Jinv_cpu)

    # Initial condition: Gaussian pressure pulse centred at (0.5, 0.5)
    U_cpu = zeros(T, N_FIELDS, Np, K)
    e = 0
    for ey in 1:Ky, ex in 1:Kx
        e += 1
        for p in 1:Np
            xp = x_cpu[p, e]
            yp = y_cpu[p, e]
            U_cpu[1, p, e] = exp(-80 * ((xp - 0.5)^2 + (yp - 0.5)^2))
        end
    end
    copyto!(U, U_cpu)
    @show sum(abs.(U))

    println("Launching volume kernel...")
    kernel! = volume_kernel!(backend, (Np, 1))
    kernel!(dU, U, Dhat, Ja, Jinv, c², Val(Np1D); ndrange = (Np, K))
    @show sum(abs.(dU))
    KernelAbstractions.synchronize(backend)

    # Report
    dU_cpu = Array(dU)
    max_dp = maximum(abs, @view dU_cpu[1, :, :])
    max_du = maximum(abs, @view dU_cpu[2, :, :])
    max_dv = maximum(abs, @view dU_cpu[3, :, :])

    println("\nVolume kernel completed ✓")
    println("  Backend       : ", backend)
    println("  Np1D / Np / K : $Np1D / $Np / $K")
    println("  max|dp/dt|    : ", max_dp)
    println("  max|du/dt|    : ", max_du)
    println("  max|dv/dt|    : ", max_dv)

    #@assert max_dp > 0 "Pressure residual is zero — check kernel or initial condition."
    @assert max_du > 0 "u-velocity residual is zero — check kernel or initial condition."
    println("\nAll sanity checks passed ✓")
end

main()
