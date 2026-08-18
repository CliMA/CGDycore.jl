@kernel inbounds = true function BlockThomasKernel!(
    @Const(DL),         # Untere Diagonale (gepackt 4x2):   (num_systems, 4, 2, nz-1)
    @Const(DD),         # Hauptdiagonale (4x4):            (num_systems, 4, 4, nz)
    @Const(DU),         # Obere Diagonale (gepackt 4x2):    (num_systems, 4, 2, nz-1)
    D,                  # Eingang & Ausgang (RHS / Lösung): (num_systems, 4, nz)
    ::Val{nz}, 
) where {nz} 
    iD = @index(Local, Linear)
    ID = @index(Global, Linear)
    NDloc = @uniform @groupsize()[1]
    ND = @uniform @ndrange()[1]

    
    # Lokaler On-Chip-Speicher via @localmem
    C_local = @localmem eltype(DD) (NDloc, 4, 4, nz)
    D_local = @localmem eltype(DD) (NDloc, 4, nz)

    
    if ID <= ND
        
        # --- 1. VORWÄRTSELIMINATION (i = 1) ---
        B_1 = SMatrix{4, 4, eltype(DD)}(
            DD[ID, 1, 1, 1], DD[ID, 2, 1, 1], DD[ID, 3, 1, 1], DD[ID, 4, 1, 1],
            DD[ID, 1, 2, 1], DD[ID, 2, 2, 1], DD[ID, 3, 2, 1], DD[ID, 4, 2, 1],
            DD[ID, 1, 3, 1], DD[ID, 2, 3, 1], DD[ID, 3, 3, 1], DD[ID, 4, 3, 1],
            DD[ID, 1, 4, 1], DD[ID, 2, 4, 1], DD[ID, 3, 4, 1], DD[ID, 4, 4, 1]
        )
        B_inv = inv(B_1)
        
        # DU hat nur N-1 Blöcke. Block 1 existiert immer, da N >= 2 vorausgesetzt wird.
        # Spalten 1-2 besetzt, Spalten 3-4 sind Null (Padding)
        C_1 = SMatrix{4, 4, eltype(DD)}(
            DU[ID, 1, 1, 1], DU[ID, 2, 1, 1], DU[ID, 3, 1, 1], DU[ID, 4, 1, 1],
            DU[ID, 1, 2, 1], DU[ID, 2, 2, 1], DU[ID, 3, 2, 1], DU[ID, 4, 2, 1],
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0
        )
        
        C_prime = B_inv * C_1
        D_1 = SVector{4, eltype(DD)}(D[ID, 1, 1], D[ID, 2, 1], D[ID, 3, 1], D[ID, 4, 1])
        D_prime = B_inv * D_1
        
        for j in 1:4, i in 1:4 C_local[iD, i, j, 1] = C_prime[i, j] end
        for i in 1:4 D_local[iD, i, 1] = D_prime[i] end
        
        # --- Schleife für i = 2 bis N ---
        for iz in 2:nz
            B_i = SMatrix{4, 4, eltype(DD)}(
                DD[ID, 1, 1, iz], DD[ID, 2, 1, iz], DD[ID, 3, 1, iz], DD[ID, 4, 1, iz],
                DD[ID, 1, 2, iz], DD[ID, 2, 2, iz], DD[ID, 3, 2, iz], DD[ID, 4, 2, iz],
                DD[ID, 1, 3, iz], DD[ID, 2, 3, iz], DD[ID, 3, 3, iz], DD[ID, 4, 3, iz],
                DD[ID, 1, 4, iz], DD[ID, 2, 4, iz], DD[ID, 3, 4, iz], DD[ID, 4, 4, iz]
            )
            D_i = SVector{4, eltype(DD)}(D[ID, 1, iz], D[ID, 2, iz], D[ID, 3, iz], D[ID, 4, iz])
            
            # DL_idx greift auf iz - 1 zu (läuft von 1 bis N-1)
            izm1 = iz - 1
            
            # A_i aufbauen aus DL (Spalten 1-2 sind 0, Spalten 3-4 gepackt)
            A_i = SMatrix{4, 4, eltype(DD)}(
                0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0,
                DL[ID, 1, 1, izm1], DL[ID, 2, 1, izm1], DL[ID, 3, 1, izm1], DL[ID, 4, 1, izm1],
                DL[ID, 1, 2, izm1], DL[ID, 2, 2, izm1], DL[ID, 3, 2, izm1], DL[ID, 4, 2, izm1]
            )
            
            B_mod_inv = inv(B_i - A_i * C_prime)
            D_prime = B_mod_inv * (D_i - A_i * D_prime)
            for i in 1:4 D_local[iD, i, iz] = D_prime[i] end
            
            # C_i / DU_i existiert nur bis einschließlich iz == N-1
            if iz < nz
                C_i = SMatrix{4, 4, eltype(DD)}(
                    DU[ID, 1, 1, iz], DU[ID, 2, 1, iz], DU[ID, 3, 1, iz], DU[ID, 4, 1, iz],
                    DU[ID, 1, 2, iz], DU[ID, 2, 2, iz], DU[ID, 3, 2, iz], DU[ID, 4, 2, iz],
                    0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0
                )
                C_prime = B_mod_inv * C_i
                for j in 1:4, i in 1:4 C_local[iD, i, j, iz] = C_prime[i, j] end
            end
        end
        
        # --- 2. RÜCKWÄRTSSUBSTITUTION (In-Place in D) ---
        X_curr = SVector{4, eltype(DD)}(D_local[iD, 1, nz], D_local[iD, 2, nz], D_local[iD, 3, nz], D_local[iD, 4, nz])
        for i in 1:4 D[ID, i, nz] = X_curr[i] end  
        
        for iz in (nz-1):-1:1
            C_prime = SMatrix{4, 4, eltype(DD)}(
                C_local[iD, 1, 1, iz], C_local[iD, 2, 1, iz], C_local[iD, 3, 1, iz], C_local[iD, 4, 1, iz],
                C_local[iD, 1, 2, iz], C_local[iD, 2, 2, iz], C_local[iD, 3, 2, iz], C_local[iD, 4, 2, iz],
                C_local[iD, 1, 3, iz], C_local[iD, 2, 3, iz], C_local[iD, 3, 3, iz], C_local[iD, 4, 3, iz],
                C_local[iD, 1, 4, iz], C_local[iD, 2, 4, iz], C_local[iD, 3, 4, iz], C_local[iD, 4, 4, iz]
            )
            D_prime = SVector{4, eltype(DD)}(D_local[iD, 1, iz], D_local[iD, 2, iz], D_local[iD, 3, iz], D_local[iD, 4, iz])
            X_curr = D_prime - C_prime * X_curr
            for i in 1:4 D[ID, i, iz] = X_curr[i] end  
        end
    end
end
