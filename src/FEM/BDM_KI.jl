using TypedPolynomials
using LinearAlgebra

"""
Generiert die BDM_1 duale Basis auf [-1, 1]^2 mit TypedPolynomials.
"""
function generate_bdm1_dual_basis()
    # 1. Algebraische Variablen deklarieren
    @polyvar x y
    
    # 2. Monombasis für BDM_1 definieren ((P_1)^2 + 2 zusätzliche curl-Erweiterungen)
    # Jede Funktion ist ein 2-elementiger Vektor aus echten Polynom-Objekten
    basis = [
        [1 + 0*x, 0 + 0*x],  # phi_1 (Vollständige P_1-Basis)
        [0 + 0*x, 1 + 0*x],  # phi_2
        [y,      0 + 0*x],  # phi_3
        [0 + 0*x, y     ],  # phi_4
        [x,      0 + 0*x],  # phi_5
        [0 + 0*x, x     ],  # phi_6
        [x^2,    -2*x*y],  # phi_7 (Zusatzterm 1)
        [2*x*y,  -y^2  ]   # phi_8 (Zusatzterm 2)
    ]
    
    # 3. Analytische 1D-Integration für die Kanten (da Integrationsgrenzen fix [-1, 1])
    # Hilfsfunktion, um ein univariates Polynom p(t) über [-1, 1] exakt zu integrieren
    function integrate_1d(p)
        # Bestimme das Integral durch Auswertung der Stammfunktion an den Grenzen
        P = differentiate(p, x) # Platzhalter-Integration via Koeffizienten-Shift:
        # Da TypedPolynomials Arithmetik beherrscht, werten wir die Monome direkt aus:
        total = 0.0
        for term in terms(p)
            c = coefficient(term)
            # Finde die Exponenten von x und y (einer ist fixiert, der andere die Variable)
            # Da wir die Variablen auf den Kanten durch Zahlen ersetzen, bleibt ein 1D-Polynom übrig
            deg = degree(term)
            total += c * (1^(deg + 1) - (-1)^(deg + 1)) / (deg + 1)
        end
        return total
    end

    # 4. Vandermonde-DoF-Matrix aufbauen
    # dofs[i](f) berechnet DoF_i für ein Vektorfeld f
    V = zeros(8, 8)
    
    for j in 1:8
        f = basis[j]
        
        # Kante 1 (unten): y = -1, n = [0, -1]. Fluss = -f[2]
        f_e1 = subs(f, y => -1.0)
        V[1, j] = integrate_1d(-f_e1[2])       # Moment 0 (1)
        V[2, j] = integrate_1d(-f_e1[2] * x)   # Moment 1 (x)
        
        # Kante 2 (rechts): x = 1, n =. Fluss = f[1]
        f_e2 = subs(f, x => 1.0)
        V[3, j] = integrate_1d(subs(f_e2[1], y => x))     # Variable auf y umbiegen für integration
        V[4, j] = integrate_1d(subs(f_e2[1] * y, y => x))
        
        # Kante 3 (oben): y = 1, n =. Fluss = f[2]
        f_e3 = subs(f, y => 1.0)
        V[5, j] = integrate_1d(f_e3[2])
        V[6, j] = integrate_1d(f_e3[2] * x)
        
        # Kante 4 (links): x = -1, n = [-1, 0]. Fluss = -f[1]
        f_e4 = subs(f, x => -1.0)
        V[7, j] = integrate_1d(subs(-f_e4[1], y => x))
        V[8, j] = integrate_1d(subs(-f_e4[1] * y, y => x))
    end
    
    # 5. Invertieren, um die Koeffizienten der dualen Basis zu erhalten
    V_inv = inv(V)
    
    # 6. Formfunktionen als Linearkombination zusammensetzen
    dual_basis = []
    for i in 1:8
        phi_x = sum(V_inv[j, i] * basis[j][1] for j in 1:8)
        phi_y = sum(V_inv[j, i] * basis[j][2] for j in 1:8)
        push!(dual_basis, [phi_x, phi_y])
    end
    
    return dual_basis
end

# --- Ausführen und Anzeigen ---
dual_basis = generate_bdm1_dual_basis()

kanten = ["Kante 1 (unten)", "Kante 2 (rechts)", "Kante 3 (oben)", "Kante 4 (links)"]
for i in 1:8
    k_idx = div(i-1, 2) + 1
    m_type = mod(i-1, 2) == 0 ? "Konstantes Moment" : "Lineares Moment"
    println("=== $(kanten[k_idx]) - $m_type ===")
    println("Phi_x = ", dual_basis[i][1])
    println("Phi_y = ", dual_basis[i][2], "\n")
end

