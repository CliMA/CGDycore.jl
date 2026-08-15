# ==============================================================================
# 1D Acoustic-Advection Solver via Symplectic Sub-Stepping
# ==============================================================================

# Helper function for periodic indexing
@inline p_idx(i, N) = mod1(i, N)

"""
Computes the 3rd-order upwind advection derivative for a given field.
Assumes fixed advection velocity a > 0.
"""
function compute_advection_component(field, N, dx, a)
    N_field = zero(field)
    denom = 6.0 * dx
    for i in 1:N
        # 3rd-order upwind stencil
        df_dx = (field[p_idx(i+2, N)] - 6.0*field[p_idx(i+1, N)] + 
                 3.0*field[i] + 2.0*field[p_idx(i-1, N)]) / denom
        N_field[i] = -a * df_dx
    end
    return N_field
end

"""
Applies 2nd-order central spatial differences for the acoustic operators
directly using the staggered grid logic.
"""
function compute_dp_dx!(dp, p, N, dx)
    for i in 1:N
        # u[i] sits between p[i] and p[i-1]
        dp[i] = (p[i] - p[p_idx(i-1, N)]) / dx
    end
end

function compute_du_dx!(du, u, N, dx)
    for i in 1:N
        # p[i] sits between u[i+1] and u[i]
        du[i] = (u[p_idx(i+1, N)] - u[i]) / dx
    end
end

# ==============================================================================
# Simulation Parameters & Initialization
# ==============================================================================
const N = 200
const L_domain = 2.0
const dx = L_domain / N

const rho0 = 1.0
const c = 340.0        # Acoustic speed
const a_vel = 10.0      # Advection velocity

# Time stepping parameters
const dt = 0.002       # Global macro time step
const t_end = 0.1
const M_substeps = 10  # Internal steps to ensure symplectic stability for acoustics
const dtau = dt / M_substeps

# Initial conditions
x_p = [(i - 0.5) * dx for i in 1:N]
p = exp.(-50.0 * (x_p .- 1.0).^2)
u = zeros(N)

# Allocations for derivative arrays
dp = zero(p)
du = zero(u)

t = 0.0
step = 0
println("Starting Symplectic Forward-Backward Sub-stepped ETD Simulation...")

# ==============================================================================
# Main Time Loop
# ==============================================================================
while t < t_end
    global p, u, t, step
    
    # --- Evaluate frozen advection tendencies at time level n ---
    Np_n = compute_advection_component(p, N, dx, a_vel)
    Nu_n = compute_advection_component(u, N, dx, a_vel)
    
    # --------------------------------------------------------------------------
    # SUB-STEP 1: Predictor Stage (Find p_a, u_a)
    # --------------------------------------------------------------------------
    p_v = copy(p)
    u_v = copy(u)
    
    for m in 1:M_substeps
        # 1. Update Pressure (Forward step using current velocity)
        compute_du_dx!(du, u_v, N, dx)
        @. p_v += dtau * (-rho0 * c^2 * du + Np_n)
        
        # 2. Update Velocity (Backward step using newly updated pressure)
        compute_dp_dx!(dp, p_v, N, dx)
        @. u_v += dtau * (-(1.0 / rho0) * dp + Nu_n)
    end
    
    p_a = p_v
    u_a = u_v
    
    # --- Evaluate advection tendencies at the predictor state ---
    Np_a = compute_advection_component(p_a, N, dx, a_vel)
    Nu_a = compute_advection_component(u_a, N, dx, a_vel)
    
    # Difference in advective updates
    dNp = Np_a .- Np_n
    dNu = Nu_a .- Nu_n
    
    # --------------------------------------------------------------------------
    # SUB-STEP 2: Corrector Stage (Find p_next, u_next)
    # --------------------------------------------------------------------------
    p_z = copy(p_a)
    u_z = copy(u_a)
    
    for m in 1:M_substeps
        tau = (m - 0.5) * dtau  # Midpoint evaluation for the linear time ramp
        ramp = tau / dt
        
        # 1. Update Pressure
        compute_du_dx!(du, u_z, N, dx)
        @. p_z += dtau * (-rho0 * c^2 * du + ramp * dNp)
        
        # 2. Update Velocity
        compute_dp_dx!(dp, p_z, N, dx)
        @. u_z += dtau * (-(1.0 / rho0) * dp + ramp * dNu)
    end
    
    # Update global states
    p .= p_z
    u .= u_z
    
    t += dt
    step += 1
    
    if step % 10 == 0
        println("Step: $step | Time: $(round(t, digits=4)) | Max Pressure: $(round(maximum(p), digits=4))")
    end
end

println("Simulation finished successfully using fully explicit Symplectic Splitting!")

