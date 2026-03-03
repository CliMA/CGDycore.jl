using LinearAlgebra
using PolynomialBases
using Plots

# --- Parameter ---
p = 4              # Polynomgrad
N_elements = 5     # Anzahl der Elemente
a = 1.0            # Advektionsgeschwindigkeit
flux_type = :upwind # :upwind oder :central

# --- 1. LGL Knoten und Operatoren auf dem Referenzelement [-1, 1] ---
basis = LobattoLegendre(p)
x_lgl = basis.nodes
D = basis.D        # Ableitungsmatrix
M = diagm(basis.weights) # Massenmatrix (diagonal bei LGL)

# SBP Eigenschaft: M*D + D'*M = B, wobei B = diag(-1, 0, ..., 0, 1)
B = diagm([-1.0; zeros(p-1); 1.0])

# --- 2. Aufbau des globalen Operators ---
N_nodes = p + 1
Total_nodes = N_elements * N_nodes
L = zeros(Total_nodes, Total_nodes)

# Metrik: dx/dxi = h/2
h = 2.0 / N_elements
inv_J = 2.0 / h

for e in 1:N_elements
    idx = ((e-1)*N_nodes + 1):(e*N_nodes)
    
    # Volumenterm in Split-Form: 0.5 * (D*u + M^-1 * D' * M * u)
    # Für konstantes 'a' entspricht dies der Standard-Ableitung
    V = a * inv_J * D 
    L[idx, idx] .-= V
    
    # Oberflächen-Kopplung (Flüsse)
    # Wir nutzen Upwind-Flüsse an den Elementgrenzen
    # Linker Rand (L) und Rechter Rand (R) des Elements
    L_node = idx[1]
    R_node = idx[end]
    
    # Zentraler Anteil + Dissipation (Upwind)
    beta = (flux_type == :upwind) ? 1.0 : 0.0
    num_flux_coeff = a * 0.5
    dissipation = abs(a) * 0.5 * beta
    
    # Inter-Element Kopplung (periodisch)
    prev_R = (e == 1) ? (N_elements * N_nodes) : ((e-1)*N_nodes)
    next_L = (e == N_elements) ? 1 : (e*N_nodes + 1)

    # Strafe/Flux am linken Rand des Elements
    L[L_node, L_node] -= inv_J * (1/basis.weights[1]) * (num_flux_coeff + dissipation)
    L[L_node, prev_R] += inv_J * (1/basis.weights[1]) * (num_flux_coeff + dissipation)

    # Strafe/Flux am rechten Rand des Elements
    L[R_node, R_node] += inv_J * (1/basis.weights[end]) * (num_flux_coeff - dissipation)
    L[R_node, next_L] -= inv_J * (1/basis.weights[end]) * (num_flux_coeff - dissipation)
end

# --- 3. Eigenwerte berechnen ---
ev = eigvals(L)

# --- 4. Plotten ---
scatter(real.(ev), imag.(ev), 
    title="DGSEM LGL Spektrum (p=$p, Elements=$N_elements)",
    xlabel="Re(λ) - Dissipation", ylabel="Im(λ) - Dispersion",
    label="Eigenwerte", marker=:circle, markersize=4, aspect_ratio=:equal)
vline!([0], color=:black, linestyle=:dash, label="")
