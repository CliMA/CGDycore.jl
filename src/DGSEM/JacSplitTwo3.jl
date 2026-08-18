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

    if ID <= ND
        FT = eltype(DD)
        invfac = FT(1) / fac

        # ONLY allocate CC3 in thread registers (4 x (M-2))
        # This drastically reduces register pressure on GPUs!
        CC3 = @private FT (4, M - 2)

        # ===================================================================
        # STEP A: COMPUTE CC3 = C_eff * inv(SA)
        # ===================================================================
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

            # In-place solve SA^T * CC3[row, :]^T = C_eff^T
            solve_sa_transpose!(CC3, row, SA, iz, ID, Val(M))
        end

        # ===================================================================
        # STEP B: FUSED CC * B & SCHUR COMPLEMENT UPDATE
        # ===================================================================

        # Cache scalar boundary vectors to avoid repeated global memory reads
        b1_1_val = B1_1[iz, ID]
        b1_4_val = B1_4[iz, ID]

        # -------------------------------------------------------------------
        # 1. DD Columns 1 & 4 (Uses CC1[row, 1], CC1[row, M], and CC3)
        # -------------------------------------------------------------------
        @unroll for row = 1:4
            r1_idx = (row == 2) ? 1 : 2
            
            # Compute CC1[row, 1] on the fly
            cc1_1 = (row == 2 || row == 3) ? C23_1[r1_idx, 1, iz, ID] : FT(0)
            @unroll for k = 1:M - 2
                cc1_1 -= CC3[row, k] * A31[k, 1, iz, ID]
            end
            cc1_1 *= invfac

            # Compute CC1[row, M] on the fly
            cc1_M = (row == 2 || row == 3) ? C23_1[r1_idx, M, iz, ID] : FT(0)
            @unroll for k = 1:M - 2
                cc1_M -= CC3[row, k] * A31[k, M, iz, ID]
            end
            cc1_M *= invfac

            # Accumulate dot products into scalar registers
            prod_col1 = cc1_1 * b1_1_val
            prod_col4 = cc1_M * b1_4_val
            @unroll for k = 1:M - 2
                prod_col1 += CC3[row, k] * B3_14[k, 1, iz, ID]
                prod_col4 += CC3[row, k] * B3_14[k, 2, iz, ID]
            end

            DD[ID, row, 1, iz] -= prod_col1
            DD[ID, row, 4, iz] -= prod_col4

            # ---------------------------------------------------------------
            # Off-diagonal DL & DU updates using cc1_1 and cc1_M
            # ---------------------------------------------------------------
            if iz > 1
                DL[ID, row, 1, iz - 1] -= cc1_1 * B1m_34[1, iz, ID]
                DL[ID, row, 2, iz - 1] -= cc1_1 * B1m_34[2, iz, ID]
            end

            if iz < nz
                DU[ID, row, 1, iz] -= cc1_M * B1p_12[1, iz, ID]
                DU[ID, row, 2, iz] -= cc1_M * B1p_12[2, iz, ID]
            end
        end

        # -------------------------------------------------------------------
        # 2. DD Columns 2 & 3 (Fused evaluation of CC1 and CC2)
        # -------------------------------------------------------------------
        @unroll for row = 1:4
            r1_idx = (row == 2) ? 1 : 2
            prod_col2 = FT(0)
            prod_col3 = FT(0)

            # CC1 component contribution
            @unroll for k = 1:M
                cc1_k = (row == 2 || row == 3) ? C23_1[r1_idx, k, iz, ID] : FT(0)
                @unroll for j = 1:M - 2
                    cc1_k -= CC3[row, j] * A31[j, k, iz, ID]
                end
                cc1_k *= invfac

                prod_col2 += cc1_k * B1_23[k, 1, iz, ID]
                prod_col3 += cc1_k * B1_23[k, 2, iz, ID]
            end

            # CC2 component contribution
            @unroll for k = 1:M - 2
                cc2_k = (row == 2 || row == 3) ? C23_2[r1_idx, k, iz, ID] : FT(0)
                @unroll for j = 1:M - 2
                    cc2_k -= CC3[row, j] * A32[j, k, iz, ID]
                end
                cc2_k *= invfac

                prod_col2 += cc2_k * B2_23[k, 1, iz, ID]
                prod_col3 += cc2_k * B2_23[k, 2, iz, ID]
            end

            DD[ID, row, 2, iz] -= prod_col2
            DD[ID, row, 3, iz] -= prod_col3
        end
    end
end
