@inline function solve_sa_transpose!(CC3_row, SA, iz, ID, ::Val{M}) where {M}
    # 1. Forward substitution with U^T (non-unit lower triangular)
    @unroll for k = 1 : M - 2
        acc = CC3_row[k]
        @unroll for i = 1 : k - 1
            acc -= SA[i, k, iz, ID] * CC3_row[i]
        end
        CC3_row[k] = acc / SA[k, k, iz, ID]
    end

    # 2. Backward substitution with L^T (unit upper triangular)
    @unroll for k = M - 2 : -1 : 1
        acc = CC3_row[k]
        @unroll for i = k + 1 : M - 2
            acc -= SA[i, k, iz, ID] * CC3_row[i]
        end
        CC3_row[k] = acc
    end
end

@kernel inbounds = true function SchurBoundaryKernel!(@Const(A13), @Const(A23), @Const(A31), @Const(A32),
  @Const(B1m_34), @Const(B1_1), @Const(B1_23), @Const(B1_4), @Const(B2_23), @Const(B3_14),
  @Const(B1p_12), @Const(C23_1), @Const(C23_2), @Const(C14_3), @Const(SA), DU, DD, DL, fac, ::Val{M}) where {M}

    iz, ID = @index(Global, NTuple)
    nz = @uniform @ndrange()[1]
    ND = @uniform @ndrange()[2]
    # 1D buffer storing ONLY the active row (length M - 2)
    CC3_row = @private eltype(DD) (M - 2)

    if ID <= ND
        FT = eltype(DD)
        invfac = FT(1) / fac


        # Cache invariant boundary scalar reads
        b1_1_val = B1_1[iz, ID]
        b1_4_val = B1_4[iz, ID]

        # Cache iz boundary conditions outside the loop
        is_iz_gt_1 = iz > 1
        is_iz_lt_nz = iz < nz

        # ===================================================================
        # SINGLE UNIFIED ROW LOOP (1 to 4)
        # ===================================================================
        @unroll for row = 1:4
            r1_idx = (row == 2) ? 1 : 2

            # ---------------------------------------------------------------
            # 1. Compute CC3_row = C_eff[row, :] * inv(SA)
            # ---------------------------------------------------------------
            if row == 1
                @unroll for j = 1 : M - 2
                    CC3_row[j] = C14_3[1, j, iz, ID]
                end
            elseif row == 4
                @unroll for j = 1 : M - 2
                    CC3_row[j] = C14_3[2, j, iz, ID]
                end
            else
                @unroll for j = 1 : M - 2
                    val = FT(0)
                    @unroll for k = 1:M
                        val -= invfac * C23_1[r1_idx, k, iz, ID] * A13[k, j, iz, ID]
                    end
                    @unroll for k = 1:M - 2
                        val -= invfac * C23_2[r1_idx, k, iz, ID] * A23[k, j, iz, ID]
                    end
                    CC3_row[j] = val
                end
            end

            # In-place solve SA^T * CC3_row^T = C_eff^T
            solve_sa_transpose!(CC3_row, SA, iz, ID, Val(M))

            # ---------------------------------------------------------------
            # 2. DD Columns 1 & 4 + Off-diagonal DL & DU Updates
            # ---------------------------------------------------------------
            # Compute CC1[row, 1] on the fly
            cc1_1 = (row == 2 || row == 3) ? C23_1[r1_idx, 1, iz, ID] : FT(0)
            @unroll for k = 1:M - 2
                cc1_1 -= CC3_row[k] * A31[k, 1, iz, ID]
            end
            cc1_1 *= invfac

            # Compute CC1[row, M] on the fly
            cc1_M = (row == 2 || row == 3) ? C23_1[r1_idx, M, iz, ID] : FT(0)
            @unroll for k = 1:M - 2
                cc1_M -= CC3_row[k] * A31[k, M, iz, ID]
            end
            cc1_M *= invfac

            # Accumulate dot products into scalar registers for Col 1 & 4
            prod_col1 = cc1_1 * b1_1_val
            prod_col4 = cc1_M * b1_4_val
            @unroll for k = 1:M - 2
                prod_col1 += CC3_row[k] * B3_14[k, 1, iz, ID]
                prod_col4 += CC3_row[k] * B3_14[k, 2, iz, ID]
            end

            DD[ID, row, 1, iz] -= prod_col1
            DD[ID, row, 4, iz] -= prod_col4

            if is_iz_gt_1
                DL[ID, row, 1, iz - 1] -= cc1_1 * B1m_34[1, iz, ID]
                DL[ID, row, 2, iz - 1] -= cc1_1 * B1m_34[2, iz, ID]
            end

            if is_iz_lt_nz
                DU[ID, row, 1, iz] -= cc1_M * B1p_12[1, iz, ID]
                DU[ID, row, 2, iz] -= cc1_M * B1p_12[2, iz, ID]
            end

            # ---------------------------------------------------------------
            # 3. DD Columns 2 & 3 (Fused CC1 and CC2 Evaluation)
            # ---------------------------------------------------------------
            prod_col2 = FT(0)
            prod_col3 = FT(0)

            # CC1 component contribution
            @unroll for k = 1:M
                cc1_k = (row == 2 || row == 3) ? C23_1[r1_idx, k, iz, ID] : FT(0)
                @unroll for j = 1:M - 2
                    cc1_k -= CC3_row[j] * A31[j, k, iz, ID]
                end
                cc1_k *= invfac

                prod_col2 += cc1_k * B1_23[k, 1, iz, ID]
                prod_col3 += cc1_k * B1_23[k, 2, iz, ID]
            end

            # CC2 component contribution
            @unroll for k = 1:M - 2
                cc2_k = (row == 2 || row == 3) ? C23_2[r1_idx, k, iz, ID] : FT(0)
                @unroll for j = 1:M - 2
                    cc2_k -= CC3_row[j] * A32[j, k, iz, ID]
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
