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

# Reference Diagonal Mass Matrix
M_ref = Diagonal(weights)

# ==============================================================================
# 3. MESH & EXACT ISOTHERMAL HYDROSTATIC INITIAL STATE
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

const T0 = 280.0
rho   = zeros(Np, Ne)
m     = zeros(Np, Ne)
Theta = zeros(Np, Ne)
theta_bg_grid = zeros(Np, Ne)

for e in 1:Ne
    for i in 1:Np
        z = X[i, e]
        
        # Exact Isothermal Atmosphere Hydrostatic Profile
        p_hydro   = p0 * exp(-g * z / (R_d * T0))
        rho_hydro = p_hydro / (R_d * T0)
        
        # Background potential temperature: theta = T0 * (p0 / p_hydro)^(R_d / c_p)
        theta_bg = T0 * exp(g * z / (c_p * T0))
        theta_bg_grid[i, e] = theta_bg
        
        # Thermal perturbation (Gaussian bump)
        z_c = 3000.0
        a_width = 1000.0
        dtheta = 1.0 * exp(-((z - z_c) / a_width)^2)
        
        theta_total = theta_bg + dtheta
        rho[i, e]   = rho_hydro
        m[i, e]     = 0.0
        Theta[i, e] = rho[i, e] * theta_total
    end
end

# Conservative state vector U = [rho, m, Theta]
U = cat(rho, m, Theta, dims=3)

# ==============================================================================
# 4. UNTHRESHOLDED THERMODYNAMICS & LMARS FLUX
# ==============================================================================
@inline function pressure(rho, Theta)
    return p0 * (R_d * Theta / p0)^gamma
end

@inline function speed_of_sound(rho, Theta)
#   return sqrt(gamma * pressure(rho, Theta) / rho)
    return 340.0
end

"""
Outward LMARS numerical flux along normal direction n (+1 or -1)
"""
function lmars_outward_flux(U_int, U_ext, n)
    rho_i, m_i, Theta_i = U_int[1], U_int[2], U_int[3]
    rho_e, m_e, Theta_e = U_ext[1], U_ext[2], U_ext[3]
    
    u_i_n = (m_i / rho_i) * n
    u_e_n = (m_e / rho_e) * n
    
    p_i = pressure(rho_i, Theta_i)
    p_e = pressure(rho_e, Theta_e)
    
    a_i = speed_of_sound(rho_i, Theta_i)
    a_e = speed_of_sound(rho_e, Theta_e)
    
    a_avg   = 0.5 * (a_i + a_e)
    rho_avg = 0.5 * (rho_i + rho_e)
    
    # Interface acoustic state
    p_hat   = 0.5 * (p_i + p_e) - 0.5 * a_avg * rho_avg * (u_e_n - u_i_n)
    u_hat_n = 0.5 * (u_i_n + u_e_n) - 0.5 / (a_avg * rho_avg) * (p_e - p_i)
    
    # Upwind flux along outward normal n
    if u_hat_n >= 0.0
        return [rho_i * u_hat_n, 
                rho_i * u_hat_n * (m_i / rho_i) + p_hat * n, 
                Theta_i * u_hat_n]
    else
        return [rho_e * u_hat_n, 
                rho_e * u_hat_n * (m_e / rho_e) + p_hat * n, 
                Theta_e * u_hat_n]
    end
end

# ==============================================================================
# 5. EXPLICIT RIGHT-HAND SIDE (RHS) EVALUATION
# ==============================================================================
function compute_rhs(U)
    dU = zeros(size(U))
    
    for e in 1:Ne
        rho_e   = U[:, e, 1]
        m_e     = U[:, e, 2]
        Theta_e = U[:, e, 3]
        
        # Physical Fluxes inside element
        F_vol = zeros(Np, 3)
        for i in 1:Np
            u_i = m_e[i] / rho_e[i]
            p_i = pressure(rho_e[i], Theta_e[i])
            
            F_vol[i, 1] = m_e[i]
            F_vol[i, 2] = m_e[i] * u_i + p_i
            F_vol[i, 3] = Theta_e[i] * u_i
        end
        
        # Volume derivative: D_ref * F / J
        dF_vol = (D_ref * F_vol) ./ J
        
        # Pointwise Gravity Source Term
        S_vol = zeros(Np, 3)
        for i in 1:Np
            S_vol[i, 2] = -rho_e[i] * g
        end
        
        # Boundary States & Outward Surface Fluxes
        # Left boundary face (i = 1, outward normal n = -1.0)
        if e == 1
            U_ext_L = [rho_e[1], -m_e[1], Theta_e[1]]
            F_hat_L = lmars_outward_flux(U[1, e, :], U_ext_L, -1.0)
        else
            F_hat_L = lmars_outward_flux(U[1, e, :], U[Np, e-1, :], -1.0)
        end
        
        # Right boundary face (i = Np, outward normal n = +1.0)
        if e == Ne
            U_ext_R = [rho_e[Np], -m_e[Np], Theta_e[Np]]
            F_hat_R = lmars_outward_flux(U[Np, e, :], U_ext_R, 1.0)
        else
            F_hat_R = lmars_outward_flux(U[Np, e, :], U[1, e+1, :], 1.0)
        end
        
        # Weak DG Assembly with addition (+) of surface jump penalty terms
        for var in 1:3
            vol_term = -dF_vol[:, var] + S_vol[:, var]
            
            # Physical flux projected along outward normal
            F_phys_L_n = -F_vol[1, var]
            F_phys_R_n =  F_vol[Np, var]
            
            surf_L = zeros(Np); surf_L[1]  = (F_phys_L_n - F_hat_L[var]) / (weights[1] * J)
            surf_R = zeros(Np); surf_R[Np] = (F_phys_R_n - F_hat_R[var]) / (weights[Np] * J)
            
            # Surface jump terms added to volume integral
            dU[:, e, var] = vol_term + surf_L + surf_R
        end
    end
    
    return dU
end

# ==============================================================================
# 6. TIME INTEGRATION (Explicit RK4) & MAIN LOOP
# ==============================================================================
c_sound_max = 340.0
min_gll_spacing = nodes[2] - nodes[1]
dx_min = J * min_gll_spacing

dt = 0.4 * dx_min / c_sound_max
t_final = 10000.0
nt = ceil(Int, t_final / dt)

println("Simulating $(Ne) elements (p=$(p)) up to t = $(t_final)s with dt = $(dt)s ($(nt) steps)...")

# Explicit RK4 Loop
for step in 1:nt
    k1 = compute_rhs(U)
    k2 = compute_rhs(U .+ 0.5 .* dt .* k1)
    k3 = compute_rhs(U .+ 0.5 .* dt .* k2)
    k4 = compute_rhs(U .+ dt .* k3)
    
    global U += (dt / 6.0) .* (k1 .+ 2.0 .* k2 .+ 2.0 .* k3 .+ k4)
end

# Final potential temperature perturbation
theta_final = U[:, :, 3] ./ U[:, :, 1]
dtheta_final = theta_final .- theta_bg_grid

# Plot result
plot(X[:], dtheta_final[:], label="Final Perturbation (t=$(t_final)s)", 
     xlabel="Height [m]", ylabel="Δθ [K]", lw=2, legend=:topright)
