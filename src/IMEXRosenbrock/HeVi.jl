# --- Struktur für einen GLL-Knoten in der Kugelschale ---
struct GLLNode
    x::Float64; y::Float64; z::Float64  # Kartesische Position (Erdzentrum)
    r::Float64                          # Radius
    n_r::Vector{Float64}                # Radiale Einheitsvektor [nx, ny, nz]
end

# Hilfsfunktion: Erstelle Geometrie für einen Knoten
function create_node(x, y, z)
    r = sqrt(x^2 + y^2 + z^2)
    return GLLNode(x, y, z, r, [x/r, y/r, z/r])
end

# --- Die HEVI-Rekonstruktions-Logik ---
"""
    reconstruct_cartesian!(state_cart, q_red_new, nodes)

Überträgt die Ergebnisse des impliziten 3x3 Solvers (rho, m_r, E) 
zurück in den 5-komponentigen kartesischen Zustandsvektor.
"""
function reconstruct_cartesian!(state_cart, q_red_new, nodes)
    n_nodes = length(nodes)
    
    for i in 1:n_nodes
        node = nodes[i]
        
        # 1. Neue Werte aus dem 3x3 Solver extrahieren
        rho_new = q_red_new[(i-1)*3 + 1]
        mr_new  = q_red_new[(i-1)*3 + 2]
        E_new   = q_red_new[(i-1)*3 + 3]
        
        # 2. Aktuellen kartesischen Impuls extrahieren (aus dem expliziten Vorhersageschritt)
        # state_cart ist hier [rho, mx, my, mz, E]
        mx_exp = state_cart[i, 2]
        my_exp = state_cart[i, 3]
        mz_exp = state_cart[i, 4]
        
        # 3. Berechne den aktuellen radialen Impulsanteil (vor der Korrektur)
        mr_exp = mx_exp * node.n_r[1] + my_exp * node.n_r[2] + mz_exp * node.n_r[3]
        
        # 4. Differenz-Korrektur (Delta des radialen Impulses)
        dm_r = mr_new - mr_exp
        
        # 5. Kartesische Korrektur: Addiere das Delta nur in radialer Richtung
        # Tangentialer Wind bleibt exakt so, wie er aus dem expliziten Schritt kam!
        state_cart[i, 1] = rho_new
        state_cart[i, 2] = mx_exp + dm_r * node.n_r[1]
        state_cart[i, 3] = my_exp + dm_r * node.n_r[2]
        state_cart[i, 4] = mz_exp + dm_r * node.n_r[3]
        state_cart[i, 5] = E_new
    end
end

# --- Beispielhafte Anwendung im Zeitschritt ---

# 1. Erstelle Test-Knoten für eine vertikale Säule
nodes = [create_node(6371e3, 0.0, z) for z in range(0.0, 10000.0, length=4)]

# 2. Kartesischer Zustand nach explizitem Horizontal-Schritt (Vorhersage)
# Nehmen wir an: [rho, mx, my, mz, E]
state_exp = ones(4, 5) 

# 3. Berechne mr_exp für den RHS des Solvers
mr_exp_vec = [dot(state_exp[i, 2:4], nodes[i].n_r) for i in 1:4]
rhs_red = vcat([[state_exp[i,1], mr_exp_vec[i], state_exp[i,5]] for i in 1:4]...)

# 4. Löse das implizite System S * q_new = rhs
# q_red_new = S_inv \ rhs_red (S_inv ist der vorfaktorisierte Block-Solver)
q_red_new = rhs_red * 1.01 # Dummy-Lösung

# 5. Zurück in kartesische Welt
reconstruct_cartesian!(state_exp, q_red_new, nodes)

println("Kartesischer Impuls nach radialer Korrektur:")
display(state_exp[:, 2:4])

