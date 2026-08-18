@inline function solve_lu_transpose_inplace!(W3, row, SA, iz, ID, ::Val{M}) where {M}
    # 1. Forward substitution with U^T (in-place update on W3)
    @unroll for k = 1 : M - 2
        @unroll for i = 1 : k - 1
            W3[row, k] -= SA[i, k, iz, ID] * W3[row, i]
        end
        W3[row, k] /= SA[k, k, iz, ID]
    end

    # 2. Backward substitution with L^T
    @unroll for k = M - 2 : -1 : 1
        @unroll for i = k + 1 : M - 2
            W3[row, k] -= SA[k, i, iz, ID] * W3[row, i]
        end
    end
end

@kernel inbounds = true function SchurBoundaryKernel!(@Const(A13), @Const(A23), @Const(A31), @Const(A32),
  @Const(B1m_34), @Const(B1_1), @Const(B1_23), @Const(B1_4), @Const(B2_23), @Const(B3_14),
  @Const(B1p_12), @Const(C23_1), @Const(C23_2), @Const(C14_3), @Const(SA), DU, DD, DL, fac, ::Val{M}) where {M}

    iz, ID = @index(Global, NTuple)
    nz = @uniform @ndrange()[1]
    ND = @uniform @ndrange()[2]
    W3 = @private eltype(DD) (4, M - 2,)

    if ID <= ND
        FT = eltype(DD)
        invfac = FT(1) / fac


        # ===================================================================
        # Step 1: Precompute Effective Left Operator C_eff & Perform 4 Transpose Solves
        # C_eff = C14_3 - (1/fac) * (C23_1 * A13 + C23_2 * A23)
        # ===================================================================
        # ===================================================================
        @unroll for row = 1:4
            if row == 1
                @unroll for j = 1 : M - 2
                    W3[row, j] = C14_3[1, j, iz, ID]
                end
            elseif row == 4
                @unroll for j = 1 : M - 2
                    W3[row, j] = C14_3[2, j, iz, ID]
                end
            else # row == 2 or row == 3
                r1_idx = (row == 2) ? 1 : 2
                @unroll for j = 1 : M - 2
                    val = FT(0)

                    # - (1/fac) * C23_1 * A13 contribution
                    @unroll for k = 1:M
                        val -= invfac * C23_1[r1_idx, k, iz, ID] * A13[k, j, iz, ID]
                    end

                    # - (1/fac) * C23_2 * A23 contribution
                    @unroll for k = 1:M - 2
                        val -= invfac * C23_2[r1_idx, k, iz, ID] * A23[k, j, iz, ID]
                    end

                    W3[row, j] = val
                end
            end
            # In-place LU Transpose Solve directly on W3[row, :]
            solve_lu_transpose_inplace!(W3, row, SA, iz, ID, Val(M))
        end  

        # ===================================================================
        # Step 2: Compute Direct Diagonal Shifts & Contract W3 * B_eff
        # ===================================================================

        # -------------------------------------------------------------------
        # Block Columns 1 & 4 (DD)
        # -------------------------------------------------------------------
        b1_val = B1_1[iz, ID] * invfac
        bm_val = B1_4[iz, ID] * invfac

        # Direct diagonal shift from (1/fac) * C23_1 * B1
        DD[ID, 2, 1, iz] += -invfac * C23_1[1, 1, iz, ID] * B1_1[iz, ID]
        DD[ID, 3, 1, iz] += -invfac * C23_1[2, 1, iz, ID] * B1_1[iz, ID]
        DD[ID, 2, 4, iz] += -invfac * C23_1[1, M, iz, ID] * B1_4[iz, ID]
        DD[ID, 3, 4, iz] += -invfac * C23_1[2, M, iz, ID] * B1_4[iz, ID]

        # Schur complement contraction: W3 * B_eff
        @unroll for row = 1:4
            s1 = FT(0); s2 = FT(0)
            @unroll for j = 1 : M - 2
                b_eff_col1 = B3_14[j, 1, iz, ID] - invfac * A31[j, 1, iz, ID] * B1_1[iz, ID]
                b_eff_col4 = B3_14[j, 2, iz, ID] - invfac * A31[j, M, iz, ID] * B1_4[iz, ID]

                s1 -= W3[row, j] * b_eff_col1
                s2 -= W3[row, j] * b_eff_col4
            end
            DD[ID, row, 1, iz] += s1
            DD[ID, row, 4, iz] += s2
        end

        # -------------------------------------------------------------------
        # Block Columns 2 & 3 (DD)
        # -------------------------------------------------------------------
        # Direct diagonal shift from (1/fac) * (C23_1 * B1_23 + C23_2 * B2_23)
        @unroll for j = 1 : M
            DD[ID, 2, 2, iz] += -invfac * C23_1[1, j, iz, ID] * B1_23[j, 1, iz, ID]
            DD[ID, 3, 2, iz] += -invfac * C23_1[2, j, iz, ID] * B1_23[j, 1, iz, ID]
            DD[ID, 2, 3, iz] += -invfac * C23_1[1, j, iz, ID] * B1_23[j, 2, iz, ID]
            DD[ID, 3, 3, iz] += -invfac * C23_1[2, j, iz, ID] * B1_23[j, 2, iz, ID]
        end
        @unroll for j = 1 : M - 2
            DD[ID, 2, 2, iz] += -invfac * C23_2[1, j, iz, ID] * B2_23[j, 1, iz, ID]
            DD[ID, 3, 2, iz] += -invfac * C23_2[2, j, iz, ID] * B2_23[j, 1, iz, ID]
            DD[ID, 2, 3, iz] += -invfac * C23_2[1, j, iz, ID] * B2_23[j, 2, iz, ID]
            DD[ID, 3, 3, iz] += -invfac * C23_2[2, j, iz, ID] * B2_23[j, 2, iz, ID]
        end

        # Schur complement contraction: W3 * B_eff
        @unroll for row = 1:4
            s1 = FT(0); s2 = FT(0)
            @unroll for j = 1 : M - 2
                b_eff_col2 = FT(0); b_eff_col3 = FT(0)
                @unroll for k = 1:M
                    b_eff_col2 -= A31[j, k, iz, ID] * B1_23[k, 1, iz, ID]
                    b_eff_col3 -= A31[j, k, iz, ID] * B1_23[k, 2, iz, ID]
                end
                @unroll for k = 1:M - 2
                    b_eff_col2 -= A32[j, k, iz, ID] * B2_23[k, 1, iz, ID]
                    b_eff_col3 -= A32[j, k, iz, ID] * B2_23[k, 2, iz, ID]
                end

                s1 -= W3[row, j] * (b_eff_col2 * invfac)
                s2 -= W3[row, j] * (b_eff_col3 * invfac)
            end
            DD[ID, row, 2, iz] += s1
            DD[ID, row, 3, iz] += s2
        end

        # -------------------------------------------------------------------
        # Lower Off-Diagonal Coupling (DL, iz > 1)
        # -------------------------------------------------------------------
        if iz > 1
            # Direct diagonal shift
            DL[ID, 2, 1, iz - 1] += -invfac * C23_1[1, 1, iz, ID] * B1m_34[1, iz, ID]
            DL[ID, 3, 1, iz - 1] += -invfac * C23_1[2, 1, iz, ID] * B1m_34[1, iz, ID]
            DL[ID, 2, 2, iz - 1] += -invfac * C23_1[1, 1, iz, ID] * B1m_34[2, iz, ID]
            DL[ID, 3, 2, iz - 1] += -invfac * C23_1[2, 1, iz, ID] * B1m_34[2, iz, ID]

            # Schur complement contraction
            b1m1 = B1m_34[1, iz, ID] * invfac
            b1m2 = B1m_34[2, iz, ID] * invfac
            @unroll for row = 1:4
                s1 = FT(0); s2 = FT(0)
                @unroll for j = 1 : M - 2
                    b_eff1 = -A31[j, 1, iz, ID] * b1m1
                    b_eff2 = -A31[j, 1, iz, ID] * b1m2
                    s1 -= W3[row, j] * b_eff1
                    s2 -= W3[row, j] * b_eff2
                end
                DL[ID, row, 1, iz - 1] += s1
                DL[ID, row, 2, iz - 1] += s2
            end
        end

        # -------------------------------------------------------------------
        # Upper Off-Diagonal Coupling (DU, iz < nz)
        # -------------------------------------------------------------------
        if iz < nz
            # Direct diagonal shift
            DU[ID, 2, 1, iz] += -invfac * C23_1[1, M, iz, ID] * B1p_12[1, iz, ID]
            DU[ID, 3, 1, iz] += -invfac * C23_1[2, M, iz, ID] * B1p_12[1, iz, ID]
            DU[ID, 2, 2, iz] += -invfac * C23_1[1, M, iz, ID] * B1p_12[2, iz, ID]
            DU[ID, 3, 2, iz] += -invfac * C23_1[2, M, iz, ID] * B1p_12[2, iz, ID]

            # Schur complement contraction
            b1p1 = B1p_12[1, iz, ID] * invfac
            b1p2 = B1p_12[2, iz, ID] * invfac
            @unroll for row = 1:4
                s1 = FT(0); s2 = FT(0)
                @unroll for j = 1 : M - 2
                    b_eff1 = -A31[j, M, iz, ID] * b1p1
                    b_eff2 = -A31[j, M, iz, ID] * b1p2
                    s1 -= W3[row, j] * b_eff1
                    s2 -= W3[row, j] * b_eff2
                end
                DU[ID, row, 1, iz] += s1
                DU[ID, row, 2, iz] += s2
            end
        end
    end
end
