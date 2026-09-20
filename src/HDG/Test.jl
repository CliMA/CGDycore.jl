using LinearAlgebra
using SparseArrays
using FastGaussQuadrature
using Plots

# ==============================================================================
# 1. PHYSICAL CONSTANTS & DOMAIN SETUP
# ==============================================================================
const g     = 9.81              # Gravitational acceleration [m/s²]
const R_d   = 287.05            # Gas constant for dry air [J/(kg K)]
const c_p   = 1004.6            # Specific heat at constant pressure [J/(kg K)]
const gamma = 1.4               # Heat capacity ratio
const p0    = 100000.0          # Reference pressure [Pa]

# Domain
const H     = 10000.0           # Height domain [m]
const Ne    = 40                # Number of elements
const p     = 4                 # Polynomial degree (p=4 -> 5 GLL nodes/elem)
const Np    = p + 1

# ==============================================================================
# 2. QUADRATURE & REFERENCE ELEMENT OPERATORS
# ==============================================================================
# GLL Nodes and Weights on [-1, 1]
nodes, weights = gausslobatto(Np)

# Lagrange basis function values and derivatives at GLL nodes
function lagrange_basis(xi, nodes)
    n = length(nodes)
    phi = zeros(n)
    dphi = zeros(n)
    for i in 1:n
        num = 1.0; den = 1.0
        sum_diff = 0.0
        for j in 1:n
            if j != i
                num *= (xi - nodes[j])
                den *= (nodes[i] - nodes[j])
                
                # Derivative product rule sum
                prod_k = 1.0
                for k in 1:n
                    if k != i && k != j
                        prod_k *= (xi - nodes[k])
                    end
                end
                sum_diff += prod_k
            end
        end
        phi[i] = num / den
        dphi[i] = sum_diff / den
    end
    return phi, dphi
end

# Differentiation Matrix D_ref (dphi_j/dxi at node i)
D_ref = zeros(Np, Np)
for i in 1:Np
    _, dphi = lagrange_basis(nodes[i], nodes)
    D_ref[i, :] = dphi
end

# Reference Diagonal Mass Matrix (due to GLL quadrature)
M_ref = Diagonal(weights)

# ==============================================================================
# 3. MESH & INITIAL STATE (Hydrostatic Atmosphere + Perturbation)
# ==============================================================================
x_faces = range(0.0, stop=H, length=Ne+1)
dx = H / Ne
J = dx / 2.0 # Coordinate transformation Jacobian (dx/dxi)

# Physical grid coordinates [Np, Ne]
X = zeros(Np, Ne)
for e in 1:Ne
    for i in 1:Np
        X[i, e] = x_faces[e] + (nodes[i] + 1.0) * J
    end
end

# Initial Hydrostatic State with T0 = 280 K
T0 = 280.0
rho   = zeros(Np, Ne)
m     = zeros(Np, Ne)
Theta = zeros(Np, Ne)
theta_bg_grid = zeros(Np, Ne)

for e in 1:Ne
    for i in 1:Np
        z = X[i, e]
        # Hydrostatic Pressure and Density profiles
        p_hydro = p0 * (1.0 - (g * z) / (c_p * T0))^(c_p / R_d)
        exner = (p_hydro / p0)^(R_d / c_p)
        theta_bg = T0 / exner
        theta_bg_grid[i, e] = theta_bg
        
        # Localized potential temperature perturbation (Gaussian bump)
        z_c = 3000.0
        a = 1000.0
        dtheta = 0.0 * exp(-((z - z_c) / a)^2)
        
        theta_total = theta_bg + dtheta
        rho[i, e]   = p_hydro / (R_d * T0)
        m[i, e]     = 0.0
        Theta[i, e] = rho[i, e] * theta_total
    end
end

# Group conservative variables U = [rho, m, Theta]
U = cat(rho, m, Theta, dims=3)

# ==============================================================================
# 4. THERMODYNAMICS & LMARS RIEMANN FLUX
# ==============================================================================
@inline function pressure(rho, Theta)
    # Positivity guard for intermediate numerical stage evaluations
    return p0 * (R_d * Theta / p0)^gamma
end

@inline function speed_of_sound(rho, Theta)
    rho_pos = max(1e-10, rho)
    p_val   = pressure(rho_pos, Theta)
    return 360.0
end

# 1D LMARS Numerical Flux at interface
function lmars_flux(U_L, U_R, n_L)
    rho_L, m_L, Theta_L = U_L[1], U_L[2], U_L[3]
    rho_R, m_R, Theta_R = U_R[1], U_R[2], U_R[3]
    
    u_L = m_L / rho_L
    u_R = m_R / rho_R
    p_L = pressure(rho_L, Theta_L)
    p_R = pressure(rho_R, Theta_R)
    
    a_L = speed_of_sound(rho_L, Theta_L)
    a_R = speed_of_sound(rho_R, Theta_R)
    a_avg = 0.5 * (a_L + a_R)
    rho_avg = 0.5 * (rho_L + rho_R)
    
    # Acoustic interface values
    p_hat = 0.5 * (p_L + p_R) - 0.5 * a_avg * rho_avg * (u_R - u_L) #* n_L
    u_hat = 0.5 * (u_L + u_R) - 0.5 / (a_avg * rho_avg) * (p_R - p_L) #* n_L
    
    # Upwind flux evaluation based on interface velocity sign
#   if u_hat * n_L >= 0.0
    if u_hat  >= 0.0
        return [rho_L * u_hat, m_L * u_hat + p_hat, Theta_L * u_hat]
    else
        return [rho_R * u_hat, m_R * u_hat + p_hat, Theta_R * u_hat]
    end
end

# ==============================================================================
# 5. EXPLICIT RIGHT-HAND SIDE (RHS) EVALUATION
# ==============================================================================
function compute_rhs(U)
    dU = zeros(size(U))
    
    # Loop over all elements
    for e in 1:Ne
        # Local state extract
        rho_e   = U[:, e, 1]
        m_e     = U[:, e, 2]
        Theta_e = U[:, e, 3]
        
        # Physical Fluxes inside volume
        F_vol = zeros(Np, 3)
        for i in 1:Np
            u_i = m_e[i] / rho_e[i]
            p_i = pressure(rho_e[i], Theta_e[i])
            
            F_vol[i, 1] = m_e[i]
            F_vol[i, 2] = m_e[i] * u_i + p_i
            F_vol[i, 3] = Theta_e[i] * u_i
        end
        
        # Volume Derivative Term: D_ref * F_vol / J
        dF_vol = (D_ref * F_vol) ./ J
        
        # Source Term (Gravity acting on momentum equation)
        S_vol = zeros(Np, 3)
        for i in 1:Np
            S_vol[i, 2] = -rho_e[i] * g
        end
        
        # Boundary States & Surface Fluxes
        # Left boundary face (i = 1, normal = -1)
        if e == 1
            # Free-slip boundary condition: reflect momentum
            U_ext_L = [rho_e[1], -m_e[1], Theta_e[1]]
            F_hat_L = lmars_flux(U_ext_L, U[1, e, :], -1.0)
        else
            F_hat_L = lmars_flux(U[Np, e-1, :], U[1, e, :], -1.0)
        end
        
        # Right boundary face (i = Np, normal = +1)
        if e == Ne
            # Free-slip boundary condition: reflect momentum
            U_ext_R = [rho_e[Np], -m_e[Np], Theta_e[Np]]
            F_hat_R = lmars_flux(U[Np, e, :], U_ext_R, 1.0)
        else
            F_hat_R = lmars_flux(U[Np, e, :], U[1, e+1, :], 1.0)
        end
        
        # Weak DG Surface penalty assembly
        @show e
        @show F_hat_L[1]
        @show F_hat_R[1]
        @show F_hat_L[2]
        @show F_hat_R[2]
        @show F_hat_L[3]
        @show F_hat_R[3]
        for var in 1:3
            # Volume divergence contribution
            vol_term = -dF_vol[:, var] + S_vol[:, var]
            @show vol_term 
            @show -dF_vol[:, var]
            @show S_vol[:, var]
            
            # Surface lifting vectors
            surf_L = zeros(Np); surf_L[1]  = (F_vol[1, var] - F_hat_L[var]) / (weights[1] * J)
            surf_R = zeros(Np); surf_R[Np] = (F_vol[Np, var] - F_hat_R[var]) / (weights[Np] * J)
            
            dU[:, e, var] = vol_term - surf_L - surf_R
        end
        @show dU[:, e, 1]
        @show dU[:, e, 2]
        @show dU[:, e, 3]
    end
    stop
    
    return dU
end

# ==============================================================================
# 6. TIME INTEGRATION (Explicit RK4) & MAIN LOOP
# ==============================================================================
c_sound_max = 340.0
min_gll_spacing = nodes[2] - nodes[1]  # ~ 0.172 for p=4
dx_min = J * min_gll_spacing


dt = 0.01 * dx_min / c_sound_max  # Courant stability limit for explicit GLL acoustic waves
@show dx_min,dt
t_final = 0.4 #10.0
nt = ceil(Int, t_final / dt)
nt = 180

println("Simulating $(Ne) elements (p=$(p)) up to t = $(t_final)s with dt = $(dt)s ($(nt) steps)...")

# Explicit RK4 Loop
for step in 1:nt
    @show step
    k1 = compute_rhs(U)
    k2 = compute_rhs(U .+ 0.5 .* dt .* k1)
    k3 = compute_rhs(U .+ 0.5 .* dt .* k2)
    k4 = compute_rhs(U .+ dt .* k3)
    
    global U += (dt / 6.0) .* (k1 .+ 2.0 .* k2 .+ 2.0 .* k3 .+ k4)
end

# Compute final potential temperature perturbation (theta_final - theta_background)
theta_final = U[:, :, 3] ./ U[:, :, 1]
dtheta_final = theta_final .- theta_bg_grid

# Plot initial vs final potential temperature perturbation
plot(X[:], dtheta_final[:], label="Final Perturbation (t=$(t_final)s)", 
     xlabel="Height [m]", ylabel="Δθ [K]", lw=2, legend=:topright)
