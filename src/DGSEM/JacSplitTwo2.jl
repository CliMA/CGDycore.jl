# Solve SA^T * CC_row[row, :]^T = C_eff^T
@inline function solve_sa_transpose!(CC_row, row, SA, iz, ID, ::Val{M}) where {M}
    # 1. Forward substitution with U^T (non-unit lower triangular matrix)
    @unroll for k = 1 : M - 2
        acc = CC_row[row, k]

        # Subtract off-diagonal terms: U^T[k, i] = U[i, k] = SA[i, k]
        # CC_row[row, i] is y_i, which was already scaled by 1 / U[i, i] in step i
        @unroll for i = 1 : k - 1
            acc -= SA[i, k, iz, ID] * CC_row[row, i]
        end

        # Divide the net accumulation by diagonal element U[k, k] = SA[k, k]
        CC_row[row, k] = acc / SA[k, k, iz, ID]
    end

    # 2. Backward substitution with L^T (unit upper triangular matrix)
    @unroll for k = M - 2 : -1 : 1
        acc = CC_row[row, k]

        # Subtract off-diagonal terms: L^T[k, i] = L[i, k] = SA[i, k] (for i > k)
        @unroll for i = k + 1 : M - 2
            acc -= SA[i, k, iz, ID] * CC_row[row, i]
        end

        # No division needed here since L^T has 1s on its diagonal
        CC_row[row, k] = acc
    end
end

@kernel inbounds = true function SchurBoundaryKernel!(@Const(A13), @Const(A23), @Const(A31), @Const(A32),
  @Const(B1m_34), @Const(B1_1), @Const(B1_23), @Const(B1_4), @Const(B2_23), @Const(B3_14),
  @Const(B1p_12), @Const(C23_1), @Const(C23_2), @Const(C14_3), @Const(SA), DU, DD, DL, fac, ::Val{M}) where {M}

    iz, ID = @index(Global, NTuple)
    nz = @uniform @ndrange()[1]
    ND = @uniform @ndrange()[2]

    CC1 = @private eltype(DD) (4, M,)
    CC2 = @private eltype(DD) (4, M - 2,)
    CC3 = @private eltype(DD) (4, M - 2,)

    if ID <= ND
        FT = eltype(DD)
        invfac = FT(1) / fac

        # Allocate private memory for CC components
        # Block structure: CC = [CC1 | CC2 | CC3]
        # CC1: 4 x M,  CC2: 4 x (M-2),  CC3: 4 x (M-2)

        # ===================================================================
        # STEP A: COMPUTE CC = C * inv(A)
        # ===================================================================

        # 1. Compute CC3 by solving SA^T * CC3^T = C_eff^T
        @unroll for row = 1:4
            if row == 1
                @unroll for j = 1 : M - 2
                    CC3[1, j] = C14_3[1, j, iz, ID]
                end
            elseif row == 4
                @unroll for j = 1 : M - 2
                    CC3[4, j] = C14_3[2, j, iz, ID]
                end
            else
                r1_idx = (row == 2) ? 1 : 2
                @unroll for j = 1 : M - 2
                    val = FT(0)
                    @unroll for k = 1:M
                        val -= invfac * C23_1[r1_idx, k, iz, ID] * A13[k, j, iz, ID]
                    end
                    @unroll for k = 1:M - 2
                        val -= invfac * C23_2[r1_idx, k, iz, ID] * A23[k, j, iz, ID]
                    end
                    CC3[row, j] = val
                end
            end

            # In-place solve: CC3[row, :] = C_eff[row, :] * inv(SA)
            solve_sa_transpose!(CC3, row, SA, iz, ID, Val(M))
        end

        # 2. Compute CC1 = (C1 - CC3 * A31) / fac
        @unroll for row = 1:4
            r1_idx = (row == 2) ? 1 : 2
            @unroll for j = 1:M
                val = FT(0)
                if row == 2 || row == 3
                    val += C23_1[r1_idx, j, iz, ID]
                end
                @unroll for k = 1:M - 2
                    val -= CC3[row, k] * A31[k, j, iz, ID]
                end
                CC1[row, j] = val * invfac
            end
        end

        # 3. Compute CC2 = (C2 - CC3 * A32) / fac
        @unroll for row = 1:4
            r1_idx = (row == 2) ? 1 : 2
            @unroll for j = 1:M - 2
                val = FT(0)
                if row == 2 || row == 3
                    val += C23_2[r1_idx, j, iz, ID]
                end
                @unroll for k = 1:M - 2
                    val -= CC3[row, k] * A32[k, j, iz, ID]
                end
                CC2[row, j] = val * invfac
            end
        end

        # ===================================================================
        # STEP B & C: COMPUTE CC * B & INSERT INTO DD, DL, DU
        # ===================================================================

        # -------------------------------------------------------------------
        # 1. Diagonal Block DD: Columns 1 & 4
        # -------------------------------------------------------------------
        @unroll for row = 1:4
            # Col 1: B = [B1_1; 0; B3_14[:, 1]]
            prod_col1 = CC1[row, 1] * B1_1[iz, ID]
            @unroll for k = 1:M - 2
                prod_col1 += CC3[row, k] * B3_14[k, 1, iz, ID]
            end

            # Col 4: B = [B1_4; 0; B3_14[:, 2]]
            prod_col4 = CC1[row, M] * B1_4[iz, ID]
            @unroll for k = 1:M - 2
                prod_col4 += CC3[row, k] * B3_14[k, 2, iz, ID]
            end

            DD[ID, row, 1, iz] -= prod_col1
            DD[ID, row, 4, iz] -= prod_col4
        end

        # -------------------------------------------------------------------
        # 2. Diagonal Block DD: Columns 2 & 3
        # -------------------------------------------------------------------
        @unroll for row = 1:4
            # Col 2: B = [B1_23[:, 1]; B2_23[:, 1]; 0]
            prod_col2 = FT(0)
            @unroll for k = 1:M
                prod_col2 += CC1[row, k] * B1_23[k, 1, iz, ID]
            end
            @unroll for k = 1:M - 2
                prod_col2 += CC2[row, k] * B2_23[k, 1, iz, ID]
            end

            # Col 3: B = [B1_23[:, 2]; B2_23[:, 2]; 0]
            prod_col3 = FT(0)
            @unroll for k = 1:M
                prod_col3 += CC1[row, k] * B1_23[k, 2, iz, ID]
            end
            @unroll for k = 1:M - 2
                prod_col3 += CC2[row, k] * B2_23[k, 2, iz, ID]
            end

            DD[ID, row, 2, iz] -= prod_col2
            DD[ID, row, 3, iz] -= prod_col3
        end

        # -------------------------------------------------------------------
        # 3. Lower Block DL (iz > 1)
        # -------------------------------------------------------------------
        if iz > 1
            @unroll for row = 1:4
                prod_col1 = CC1[row, 1] * B1m_34[1, iz, ID]
                prod_col2 = CC1[row, 1] * B1m_34[2, iz, ID]

                DL[ID, row, 1, iz - 1] -= prod_col1
                DL[ID, row, 2, iz - 1] -= prod_col2
            end
        end

        # -------------------------------------------------------------------
        # 4. Upper Block DU (iz < nz)
        # -------------------------------------------------------------------
        if iz < nz
            @unroll for row = 1:4
                prod_col1 = CC1[row, M] * B1p_12[1, iz, ID]
                prod_col2 = CC1[row, M] * B1p_12[2, iz, ID]

                DU[ID, row, 1, iz] -= prod_col1
                DU[ID, row, 2, iz] -= prod_col2
            end
        end
    end
end
