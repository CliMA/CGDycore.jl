using LinearAlgebra
using SparseArrays
using FastGaussQuadrature
using Plots

# ==============================================================================
# 1. PHYSICAL CONSTANTS & DOMAIN SETUP
# ==============================================================================
const g         = 9.81              # Gravitational acceleration [m/s²]
const R_d       = 287.05            # Gas constant [J/(kg K)]
const c_p       = 1004.6            # Specific heat [J/(kg K)]
const gamma_gas = 1.4               # Heat capacity ratio
const p0        = 100000.0          # Reference pressure [Pa]

const H         = 10000.0           # Height domain [m]
const Ne        = 40                # Number of elements
const p         = 4                 # Polynomial degree
const Np        = p + 1             # Local volume nodes per element
const Nfaces    = Ne + 1           # Total unique global trace interfaces

nodes, weights = gausslobatto(Np)

function lagrange_basis(xi, nodes)
    n = length(nodes)
    phi, dphi = zeros(n), zeros(n)
    for i in 1:n
        num, den, sum_diff = 1.0, 1.0, 0.0
        for j in 1:n
            if j != i
                num *= (xi - nodes[j]); den *= (nodes[i] - nodes[j])
                prod_k = 1.0
                for k in 1:n
                    if k != i && k != j; prod_k *= (xi - nodes[k]); end
                end
                sum_diff += prod_k
            end
        end
        phi[i], dphi[i] = num / den, sum_diff / den
    end
    return phi, dphi
end

D_ref = zeros(Np, Np)
for i in 1:Np
    _, dphi = lagrange_basis(nodes[i], nodes)
    D_ref[i, :] = dphi
end

# ==============================================================================
# 2. DATA STRUCTURES & THERMODYNAMICS
# ==============================================================================
struct HDGState
    U::Array{Float64, 3}        # Volume conservative variables [Np, Ne, 3]
    Lambda::Matrix{Float64}     # Interface trace variables [3, Nfaces]
end

struct HDGBlockJacobian
    A_inv::Vector{Matrix{Float64}}        # Local inverse volume block A_e^-1 (3Np x 3Np)
    B_vol::Vector{Matrix{Float64}}        # Trace-to-Volume coupling B_e (3Np x 6)
    C_trace::Vector{Matrix{Float64}}      # Volume-to-Trace projection C_e (6 x 3Np)
    K_global_fact::Factorization          # Pre-factorized global trace Schur complement matrix
end

@inline pressure(rho, Theta) = p0 * (R_d * Theta / p0)^gamma_gas
@inline speed_of_sound(rho, Theta) = sqrt(gamma_gas * pressure(rho, Theta) / rho)

# ==============================================================================
# 3. PURE EXPLICIT LMARS FLUX & RHS (NO LAMBDA DEPENDENCY)
# ==============================================================================
"""
Standard explicit LMARS numerical flux evaluated between left element state U_L 
and right element state U_R across an interface face.
"""
function lmars_explicit_flux(U_L::Vector{Float64}, U_R::Vector{Float64})
    rho_L, m_L, Theta_L = U_L[1], U_L[2], U_L[3]
    rho_R, m_R, Theta_R = U_R[1], U_R[2], U_R[3]
    
    u_L = m_L / rho_L; p_L = pressure(rho_L, Theta_L)
    u_R = m_R / rho_R; p_R = pressure(rho_R, Theta_R)
    
    c_s = 0.5 * (speed_of_sound(rho_L, Theta_L) + speed_of_sound(rho_R, Theta_R))
    
    F_L = [m_L, m_L * u_L + p_L, Theta_L * u_L]
    F_R = [m_R, m_R * u_R + p_R, Theta_R * u_R]
    
    # Standard LMARS central average + acoustic wave dissipation jump
    return 0.5 .* (F_L .+ F_R) .- 0.5 * c_s .* (U_R .- U_L)
end

"""
Computes the full explicit strong-form residual R(U) using pure element-to-element 
LMARS numerical fluxes.
"""
function compute_rhs_explicit!(RHS_total::Array{Float64, 3}, U::Array{Float64, 3}, J::Float64)
    # Compute interface fluxes across all internal and boundary faces (Ne + 1 faces)
    F_hat_faces = zeros(3, Ne + 1)
    
    for f in 1:(Ne + 1)
        if f == 1
            # Left boundary (reflective wall)
            U_L = [U[1, 1, 1], -U[1, 1, 2], U[1, 1, 3]]
            U_R = U[1, 1, :]
        elseif f == Ne + 1
            # Right boundary (reflective wall)
            U_L = U[Np, Ne, :]
            U_R = [U[Np, Ne, 1], -U[Np, Ne, 2], U[Np, Ne, 3]]
        else
            # Interior face between element (f-1) and element f
            U_L = U[Np, f-1, :]
            U_R = U[1, f, :]
        end
        
        F_hat_faces[:, f] = lmars_explicit_flux(U_L, U_R)
    end
    
    # Assemble strong DG residual for each element
    for e in 1:Ne
        rho_e, m_e, Theta_e = U[:, e, 1], U[:, e, 2], U[:, e, 3]
        
        # Volume flux derivative
        F_vol = zeros(Np, 3)
        for i in 1:Np
            u_i = m_e[i] / rho_e[i]
            p_i = pressure(rho_e[i], Theta_e[i])
            F_vol[i, 1] = m_e[i]
            F_vol[i, 2] = m_e[i] * u_i + p_i
            F_vol[i, 3] = Theta_e[i] * u_i
        end
        dF_vol = (D_ref * F_vol) ./ J
        
        # Boundary flux jumps
        F_hat_L = F_hat_faces[:, e]      # Flux at left face
        F_hat_R = F_hat_faces[:, e + 1]  # Flux at right face
        
        for var in 1:3
            vol_term = -dF_vol[:, var]
            if var == 2; vol_term .-= rho_e .* g; end
            
            # Strong-form surface terms: (F_phys * n - F*_hat) / (w * J)
            F_phys_L_n = -F_vol[1, var]   # n = -1
            F_phys_R_n =  F_vol[Np, var]  # n = +1
            
            surf_L = zeros(Np); surf_L[1]  = (F_phys_L_n - (-F_hat_L[var])) / (weights[1] * J)
            surf_R = zeros(Np); surf_R[Np] = (F_phys_R_n - F_hat_R[var])     / (weights[Np] * J)
            
            RHS_total[:, e, var] = vol_term + surf_L + surf_R
        end
    end
    return RHS_total
end

# ==============================================================================
# 4. EXTENDED BLOCK JACOBIAN ASSEMBLY & FACTORIZATION FOR LINEAR ALGEBRA
# ==============================================================================
function build_block_jacobian(state::HDGState, dt::Float64, gamma_dt::Float64, J::Float64)
    n_vol = Np * 3
    n_trace_elem = 6
    N_global_trace = 3 * Nfaces
    
    A_inv_list   = Vector{Matrix{Float64}}(undef, Ne)
    B_vol_list   = Vector{Matrix{Float64}}(undef, Ne)
    C_trace_list = Vector{Matrix{Float64}}(undef, Ne)
    
    I_global, J_global, V_global = Int[], Int[], Float64[]
    
    for e in 1:Ne
        f_L, f_R = e, e + 1
        rho_e, Theta_e = state.U[:, e, 1], state.U[:, e, 3]
        
        # Volume Linearization A = I - gamma_dt * (dR_vol / dU)
        A_vol = Matrix(1.0I, n_vol, n_vol)
        for i in 1:Np
            c_s2 = gamma_gas * pressure(rho_e[i], Theta_e[i]) / rho_e[i]
            dp_dTheta = (gamma_gas - 1.0) * (pressure(rho_e[i], Theta_e[i]) / Theta_e[i])
            
            for j in 1:Np
                d_ij = D_ref[i, j] / J
                row_r, col_m = (i-1)*3 + 1, (j-1)*3 + 2
                row_m, col_r, col_t = (i-1)*3 + 2, (j-1)*3 + 1, (j-1)*3 + 3
                row_t = (i-1)*3 + 3
                
                A_vol[row_r, col_m] -= gamma_dt * (-d_ij)
                A_vol[row_m, col_r] -= gamma_dt * (-d_ij * c_s2)
                A_vol[row_m, col_t] -= gamma_dt * (-d_ij * dp_dTheta)
                A_vol[row_t, col_m] -= gamma_dt * (-d_ij * (Theta_e[j] / rho_e[j]))
            end
            A_vol[(i-1)*3 + 2, (i-1)*3 + 1] -= gamma_dt * (-g)
        end
        
        # Extended Interface Coupling Blocks B and C (Blow-up with trace variables)
        B_vol = zeros(n_vol, n_trace_elem)
        C_trace = zeros(n_trace_elem, n_vol)
        
        tau_L = 0.5 * speed_of_sound(rho_e[1], Theta_e[1])
        tau_R = 0.5 * speed_of_sound(rho_e[Np], Theta_e[Np])
        
        inv_wJ_L = 1.0 / (weights[1] * J)
        inv_wJ_R = 1.0 / (weights[Np] * J)
        
        for var in 1:3
            B_vol[(1-1)*3 + var, var]      = gamma_dt * tau_L * inv_wJ_L
            B_vol[(Np-1)*3 + var, 3 + var] = gamma_dt * tau_R * inv_wJ_R
            
            C_trace[var, (1-1)*3 + var]     = tau_L
            C_trace[3+var, (Np-1)*3 + var] = tau_R
        end
        
        A_inv = inv(A_vol)
        K_elem = C_trace * A_inv * B_vol
        
        A_inv_list[e]   = A_inv
        B_vol_list[e]   = B_vol
        C_trace_list[e] = C_trace
        
        trace_indices = [ (f_L-1)*3 + var for var in 1:3 ]
        append!(trace_indices, [ (f_R-1)*3 + var for var in 1:3 ])
        
        for ti in 1:n_trace_elem, tj in 1:n_trace_elem
            push!(I_global, trace_indices[ti])
            push!(J_global, trace_indices[tj])
            push!(V_global, K_elem[ti, tj])
        end
    end
    
    for gi in 1:N_global_trace
        push!(I_global, gi); push!(J_global, gi); push!(V_global, 1.0)
    end
    
    K_global = sparse(I_global, J_global, V_global, N_global_trace, N_global_trace)
    K_global_fact = lu(K_global)
    
    return HDGBlockJacobian(A_inv_list, B_vol_list, C_trace_list, K_global_fact)
end

# ==============================================================================
# 5. GENERAL ROSENBROCK STAGE SOLVE (STATIC CONDENSATION)
# ==============================================================================
function solve_rosenbrock_stage!(k_stage::HDGState, Jac::HDGBlockJacobian, F_stage::Array{Float64,3})
    N_global_trace = 3 * Nfaces
    R_global = zeros(N_global_trace)
    
    # 1. Project explicit stage RHS onto global interface trace system: R_global = - C * A^-1 * F_stage
    for e in 1:Ne
        f_L, f_R = e, e + 1
        f_elem_vec = reshape(F_stage[:, e, :]', Np * 3)
        
        r_trace_elem = Jac.C_trace[e] * (Jac.A_inv[e] * f_elem_vec)
        
        trace_indices = [ (f_L-1)*3 + var for var in 1:3 ]
        append!(trace_indices, [ (f_R-1)*3 + var for var in 1:3 ])
        
        for ti in 1:6
            R_global[trace_indices[ti]] += r_trace_elem[ti]
        end
    end
    
    # 2. Solve global interface stage update using pre-factorized Schur Complement
    dL_flat = Jac.K_global_fact \ R_global
    k_stage.Lambda .= reshape(dL_flat, 3, Nfaces)
    
    # 3. Local Volume Back-substitution: k_U = A^-1 * (F_stage - B * k_Lambda)
    for e in 1:Ne
        f_L, f_R = e, e + 1
        dL_elem = [k_stage.Lambda[:, f_L]; k_stage.Lambda[:, f_R]]
        f_elem_vec = reshape(F_stage[:, e, :]', Np * 3)
        
        dU_elem_vec = Jac.A_inv[e] * (f_elem_vec .- Jac.B_vol[e] * dL_elem)
        k_stage.U[:, e, :] .= reshape(dU_elem_vec, 3, Np)'
    end
    return k_stage
end

# ==============================================================================
# 6. ROS2 INTEGRATOR DRIVER STEP
# ==============================================================================
function ros2_step!(state::HDGState, dt::Float64, J::Float64, 
                    RHS_workspace::Array{Float64,3}, F_stage_workspace::Array{Float64,3},
                    k1::HDGState, k2::HDGState, stage_state::HDGState)
    
    gamma_ros = 1.0 + 1.0 / sqrt(2.0)
    gamma_dt  = gamma_ros * dt
    
    # --- BUILD & FACTORIZE EXTENDED BLOCK JACOBIAN ONCE PER STEP ---
    Jac = build_block_jacobian(state, dt, gamma_dt, J)
    
    # --- STAGE 1 ---
    # Pure explicit RHS computation
    compute_rhs_explicit!(RHS_workspace, state.U, J)
    F_stage_workspace .= dt .* RHS_workspace
    solve_rosenbrock_stage!(k1, Jac, F_stage_workspace)
    
    # --- STAGE 2 ---
    # Intermediate State: U_stage = U^n + gamma * k1_U
    stage_state.U .= state.U .+ gamma_ros .* k1.U
    stage_state.Lambda .= state.Lambda .+ gamma_ros .* k1.Lambda
    
    compute_rhs_explicit!(RHS_workspace, stage_state.U, J)
    
    # f2 = dt * R_explicit(U_stage) - 2 * k1_U
    for e in 1:Ne, var in 1:3, i in 1:Np
        F_stage_workspace[i, e, var] = dt * RHS_workspace[i, e, var] - 2.0 * k1.U[i, e, var]
    end
    
    solve_rosenbrock_stage!(k2, Jac, F_stage_workspace)
    
    # --- COMBINE STAGES: U^{n+1} = U^n + 1.5 * k1 + 0.5 * k2 ---
    state.U .+= 1.5 .* k1.U .+ 0.5 .* k2.U
    state.Lambda .+= 1.5 .* k1.Lambda .+ 0.5 .* k2.Lambda
    
    return state
end

# ==============================================================================
# 7. INITIAL CONDITIONS & SIMULATION LOOP
# ==============================================================================
x_faces = range(0.0, stop=H, length=Ne+1)
dx = H / Ne
J_geom = dx / 2.0

X = zeros(Np, Ne)
for e in 1:Ne, i in 1:Np
    X[i, e] = x_faces[e] + (nodes[i] + 1.0) * J_geom
end

const T0 = 280.0
rho = zeros(Np, Ne); m = zeros(Np, Ne); Theta = zeros(Np, Ne)
theta_bg_grid = zeros(Np, Ne)

for e in 1:Ne, i in 1:Np
    z = X[i, e]
    p_hydro = p0 * exp(-g * z / (R_d * T0))
    rho_hydro = p_hydro / (R_d * T0)
    theta_bg = T0 * exp(g * z / (c_p * T0))
    theta_bg_grid[i, e] = theta_bg
    
    # Gaussian Potential Temperature Anomaly
    dtheta = 1.0 * exp(-((z - 3000.0) / 1000.0)^2)
    rho[i, e]   = rho_hydro
    m[i, e]     = 0.0
    Theta[i, e] = rho_hydro * (theta_bg + dtheta)
end

U = cat(rho, m, Theta, dims=3)
Lambda = zeros(3, Nfaces)
for f in 1:Nfaces
    Lambda[:, f] = (f == 1) ? U[1, 1, :] : (f == Nfaces ? U[Np, Ne, :] : 0.5 .* (U[Np, f-1, :] .+ U[1, f, :]))
end

state = HDGState(U, Lambda)

# Workspace Allocations
RHS_ws     = zeros(Np, Ne, 3)
F_stage_ws = zeros(Np, Ne, 3)
k1_state   = HDGState(zeros(Np, Ne, 3), zeros(3, Nfaces))
k2_state   = HDGState(zeros(Np, Ne, 3), zeros(3, Nfaces))
stg_state  = HDGState(zeros(Np, Ne, 3), zeros(3, Nfaces))

dt = 2.0
t_final = 10.0
nt = ceil(Int, t_final / dt)

println("Running ROS2 HDG Solver with Pure Explicit RHS & Extended Interface Linear Algebra...")

dt = 0.04 * 50 / 340
nt = 2000
for step in 1:nt
    @show step
    ros2_step!(state, dt, J_geom, RHS_ws, F_stage_ws, k1_state, k2_state, stg_state)
end

# Visualizing Solution
theta_final = state.U[:, :, 3] ./ state.U[:, :, 1]
dtheta_final = theta_final .- theta_bg_grid
dtheta_final = state.U[:, :, 2] ./ state.U[:, :, 1]

plot(X[:], dtheta_final[:], 
     label="ROS2 Explicit-RHS / HDG-Solves", 
     xlabel="Height [m]", 
     ylabel="Δθ [K]", 
     lw=2, 
     title="Potential Temperature Anomaly after Time Integration")
