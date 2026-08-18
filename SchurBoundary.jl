@kernel inbounds = true function SchurBoundaryKernel!(@Const(A13),@Const(A23),@Const(A31),@Const(A32),
  @Const(B1m_34),@Const(B1_1),@Const(B1_23),@Const(B1_4), @Const(B2_23),@Const(B3_14),
  @Const(B1p_12),@Const(C23_1),@Const(C23_2),@Const(C14_3),@Const(SA),DU,DD,DL,fac, ::Val{M}) where {M}

  iz,ID = @index(Global, NTuple)

  nz = @uniform @ndrange()[1]
  ND = @uniform @ndrange()[2]

  @uniform FT = eltype(DD)

  r1 = @private FT (M,2,)
  r2 = @private FT (M-2,2,)
  r3 = @private FT (M-2,2,)
  r11 = @private FT (2,)
  r1M = @private FT (2,)
  s = @private FT (4,2,)
  invfac = FT(1) / fac

  if ID <= ND
      
#   Column 1 and 4
    @unroll for i = 1 : M
      r1[i,1] = FT(0)  
      r1[i,2] = FT(0)  
    end  
    r1[1,1] = B1_1[iz,ID]
    r1[M,2] = B1_4[iz,ID]
    @unroll for i = 1 : M - 2
      r3[i,1] = B3_14[i,1,iz,ID] - invfac * A31[i,1,iz,ID] * r1[1,1]
      r3[i,2] = B3_14[i,2,iz,ID] - invfac * A31[i,M,iz,ID] * r1[M,2]
    end  
  

    @unroll for k = 1 : M - 3
      @unroll for i = k + 1 : M - 2
       r3[i,1] -= SA[i,k,iz,ID] * r3[k,1]
       r3[i,2] -= SA[i,k,iz,ID] * r3[k,2]
      end
    end
#   Backward loop
    @unroll for k = M - 2 : -1 : 1
      r3[k,1] /= SA[k,k,iz,ID]
      r3[k,2] /= SA[k,k,iz,ID]
      @unroll for i = 1 : k - 1
        r3[i,1] -= SA[i,k,iz,ID] * r3[k,1]
        r3[i,2] -= SA[i,k,iz,ID] * r3[k,2]
      end
    end

    @unroll for i = 1 : M
      @unroll for j = 1 : M - 2
        r1[i,1] -= A13[i,j,iz,ID] * r3[j,1]
        r1[i,2] -= A13[i,j,iz,ID] * r3[j,2]
      end
      r1[i,1] *= invfac
      r1[i,2] *= invfac
    end  

    @unroll for i = 1 : M - 2
      r2[i,1] = FT(0)
      r2[i,2] = FT(0)
      @unroll for j = 1 : M - 2
        r2[i,1] -= A23[i,j,iz,ID] * r3[j,1]
        r2[i,2] -= A23[i,j,iz,ID] * r3[j,2]
      end
      r2[i,1] *= invfac
      r2[i,2] *= invfac
    end

    s[1,1] = FT(0)
    s[2,1] = FT(0)
    s[3,1] = FT(0)
    s[4,1] = FT(0)
    s[1,2] = FT(0)
    s[2,2] = FT(0)
    s[3,2] = FT(0)
    s[4,2] = FT(0)
    @unroll for j = 1 : M - 2
      s[1,1] += -C14_3[1,j,iz,ID] * r3[j,1]
      s[1,2] += -C14_3[1,j,iz,ID] * r3[j,2]
      s[4,1] += -C14_3[2,j,iz,ID] * r3[j,1]
      s[4,2] += -C14_3[2,j,iz,ID] * r3[j,2]
      s[2,1] += -C23_2[1,j,iz,ID] * r2[j,1]  
      s[2,2] += -C23_2[1,j,iz,ID] * r2[j,2]
      s[3,1] += -C23_2[2,j,iz,ID] * r2[j,1]
      s[3,2] += -C23_2[2,j,iz,ID] * r2[j,2]
    end
    @unroll for j = 1 : M
      s[2,1] += -C23_1[1,j,iz,ID] * r1[j,1]  
      s[2,2] += -C23_1[1,j,iz,ID] * r1[j,2]
      s[3,1] += -C23_1[2,j,iz,ID] * r1[j,1]
      s[3,2] += -C23_1[2,j,iz,ID] * r1[j,2]
    end
    DD[ID,1,1,iz] += s[1,1]
    DD[ID,2,1,iz] += s[2,1]
    DD[ID,3,1,iz] += s[3,1]
    DD[ID,4,1,iz] += s[4,1]
    DD[ID,1,4,iz] += s[1,2]
    DD[ID,2,4,iz] += s[2,2]
    DD[ID,3,4,iz] += s[3,2]
    DD[ID,4,4,iz] += s[4,2]
  
    @unroll for i = 1 : M
      r1[i,1] = B1_23[i,1,iz,ID]
      r1[i,2] = B1_23[i,2,iz,ID]
    end  
    @unroll for i = 1 : M - 2
      r2[i,1] = B2_23[i,1,iz,ID]
      r2[i,2] = B2_23[i,2,iz,ID]
    end  
    @unroll for i = 1 : M - 2
      r3[i,1] = FT(0)
      r3[i,2] = FT(0)
      @unroll for j = 1 : M 
        r3[i,1] -= A31[i,j,iz,ID] * r1[j,1]
        r3[i,2] -= A31[i,j,iz,ID] * r1[j,2]
      end
      @unroll for j = 1 : M - 2
        r3[i,1] -= A32[i,j,iz,ID] * r2[j,1]
        r3[i,2] -= A32[i,j,iz,ID] * r2[j,2]
      end
      r3[i,1] *= invfac
      r3[i,2] *= invfac
    end

    @unroll for k = 1 : M - 3
      @unroll for i = k + 1 : M - 2
       r3[i,1] -= SA[i,k,iz,ID] * r3[k,1]
       r3[i,2] -= SA[i,k,iz,ID] * r3[k,2]
      end
    end
#   Backward loop
    @unroll for k = M - 2 : -1 : 1
      r3[k,1] /= SA[k,k,iz,ID]
      r3[k,2] /= SA[k,k,iz,ID]
      @unroll for i = 1 : k - 1
        r3[i,1] -= SA[i,k,iz,ID] * r3[k,1]
        r3[i,2] -= SA[i,k,iz,ID] * r3[k,2]
      end
    end

    @unroll for i = 1 : M 
      @unroll for j = 1 : M - 2
        r1[i,1] += -A13[i,j,iz,ID] * r3[j,1]
        r1[i,2] += -A13[i,j,iz,ID] * r3[j,2]
      end
      r1[i,1] *= invfac
      r1[i,2] *= invfac
    end
    @unroll for i = 1 : M - 2
      @unroll for j = 1 : M - 2
        r2[i,1] += -A23[i,j,iz,ID] * r3[j,1]
        r2[i,2] += -A23[i,j,iz,ID] * r3[j,2]
      end
      r2[i,1] *= invfac
      r2[i,2] *= invfac
    end

    s[1,1] = FT(0)
    s[2,1] = FT(0)
    s[3,1] = FT(0)
    s[4,1] = FT(0)
    s[1,2] = FT(0)
    s[2,2] = FT(0)
    s[3,2] = FT(0)
    s[4,2] = FT(0)
    @unroll for j = 1 : M - 2
      s[1,1] += - C14_3[1,j,iz,ID] * r3[j,1]
      s[1,2] += - C14_3[1,j,iz,ID] * r3[j,2]
      s[4,1] += - C14_3[2,j,iz,ID] * r3[j,1]
      s[4,2] += - C14_3[2,j,iz,ID] * r3[j,2]
      s[2,1] += - C23_2[1,j,iz,ID] * r2[j,1]
      s[2,2] += - C23_2[1,j,iz,ID] * r2[j,2]
      s[3,1] += - C23_2[2,j,iz,ID] * r2[j,1]
      s[3,2] += - C23_2[2,j,iz,ID] * r2[j,2]
    end
    @unroll for j = 1 : M
      s[2,1] += -C23_1[1,j,iz,ID] * r1[j,1]  
      s[2,2] += -C23_1[1,j,iz,ID] * r1[j,2]
      s[3,1] += -C23_1[2,j,iz,ID] * r1[j,1]
      s[3,2] += -C23_1[2,j,iz,ID] * r1[j,2]
    end
    DD[ID,1,2,iz] += s[1,1]
    DD[ID,2,2,iz] += s[2,1]
    DD[ID,3,2,iz] += s[3,1]
    DD[ID,4,2,iz] += s[4,1]
    DD[ID,1,3,iz] += s[1,2]
    DD[ID,2,3,iz] += s[2,2]
    DD[ID,3,3,iz] += s[3,2]
    DD[ID,4,3,iz] += s[4,2]

    if iz > 1
#     @atomic :monotonic SchurBand[7,sh-1,ID] -= FacGrav * invfac * B1m_34[1,iz,ID]
#     @atomic :monotonic SchurBand[6,sh,ID] -= FacGrav * invfac * B1m_34[2,iz,ID]
      @unroll for i = 1 : M
        r1[i,1] = FT(0)
        r1[i,2] = FT(0)
      end
      r1[1,1] = B1m_34[1,iz,ID]
      r1[1,2] = B1m_34[2,iz,ID]
      @unroll for i = 1 : M - 2
        r3[i,1] = - invfac * A31[i,1,iz,ID] * r1[1,1]
        r3[i,2] = - invfac * A31[i,1,iz,ID] * r1[1,2]
      end

      @unroll for k = 1 : M - 3
        @unroll for i = k + 1 : M - 2
          r3[i,1] -= SA[i,k,iz,ID] * r3[k,1]
          r3[i,2] -= SA[i,k,iz,ID] * r3[k,2]
        end
      end
#     Backward loop
      @unroll for k = M - 2 : -1 : 1
        r3[k,1] /= SA[k,k,iz,ID]
        r3[k,2] /= SA[k,k,iz,ID]
        @unroll for i = 1 : k - 1
          r3[i,1] -= SA[i,k,iz,ID] * r3[k,1]
          r3[i,2] -= SA[i,k,iz,ID] * r3[k,2]
        end
      end
      @unroll for i = 1 : M
        @unroll for j = 1 : M - 2
          r1[i,1] += -A13[i,j,iz,ID] * r3[j,1]
          r1[i,2] += -A13[i,j,iz,ID] * r3[j,2]
        end
        r1[i,1] *= invfac
        r1[i,2] *= invfac
      end
      @unroll for i = 1 : M - 2
        r2[i,1] = FT(0) 
        r2[i,2] = FT(0) 
        @unroll for j = 1 : M - 2
          r2[i,1] += -A23[i,j,iz,ID] * r3[j,1]
          r2[i,2] += -A23[i,j,iz,ID] * r3[j,2]
        end
        r2[i,1] *= invfac
        r2[i,2] *= invfac
      end

      s[1,1] = FT(0)
      s[2,1] = FT(0)
      s[3,1] = FT(0)
      s[4,1] = FT(0)
      s[1,2] = FT(0)
      s[2,2] = FT(0)
      s[3,2] = FT(0)
      s[4,2] = FT(0)
      @unroll for j = 1 : M - 2
        s[1,1] += - C14_3[1,j,iz,ID] * r3[j,1]
        s[1,2] += - C14_3[1,j,iz,ID] * r3[j,2]
        s[4,1] += - C14_3[2,j,iz,ID] * r3[j,1]
        s[4,2] += - C14_3[2,j,iz,ID] * r3[j,2]
        s[2,1] += - C23_2[1,j,iz,ID] * r2[j,1]
        s[2,2] += - C23_2[1,j,iz,ID] * r2[j,2]
        s[3,1] += - C23_2[2,j,iz,ID] * r2[j,1]
        s[3,2] += - C23_2[2,j,iz,ID] * r2[j,2]
      end
      @unroll for j = 1 : M
        s[2,1] += -C23_1[1,j,iz,ID] * r1[j,1]
        s[2,2] += -C23_1[1,j,iz,ID] * r1[j,2]
        s[3,1] += -C23_1[2,j,iz,ID] * r1[j,1]
        s[3,2] += -C23_1[2,j,iz,ID] * r1[j,2]
      end
      DL[ID,1,1,iz-1] += s[1,1]
      DL[ID,2,1,iz-1] += s[2,1]
      DL[ID,3,1,iz-1] += s[3,1]
      DL[ID,4,1,iz-1] += s[4,1]
      DL[ID,1,2,iz-1] += s[1,2]
      DL[ID,2,2,iz-1] += s[2,2]
      DL[ID,3,2,iz-1] += s[3,2]
      DL[ID,4,2,iz-1] += s[4,2]
    end
    if iz < nz
      #@atomic :monotonic SchurBand[2,sh+5,ID] -= FacGrav * invfac * B1p_12[1,iz,ID]
      #@atomic :monotonic SchurBand[1,sh+6,ID] -= FacGrav * invfac * B1p_12[2,iz,ID]
      @unroll for i = 1 : M
        r1[i,1] = FT(0)
        r1[i,2] = FT(0)
      end
      r1[M,1] = B1p_12[1,iz,ID]
      r1[M,2] = B1p_12[2,iz,ID]
      @unroll for i = 1 : M - 2
        r3[i,1] = - invfac * A31[i,M,iz,ID] * r1[M,1]
        r3[i,2] = - invfac * A31[i,M,iz,ID] * r1[M,2]
      end

      @unroll for k = 1 : M - 3
        @unroll for i = k + 1 : M - 2
          r3[i,1] -= SA[i,k,iz,ID] * r3[k,1]
          r3[i,2] -= SA[i,k,iz,ID] * r3[k,2]
        end
      end
#     Backward loop
      @unroll for k = M - 2 : -1 : 1
        r3[k,1] /= SA[k,k,iz,ID]
        r3[k,2] /= SA[k,k,iz,ID]
        @unroll for i = 1 : k - 1
          r3[i,1] -= SA[i,k,iz,ID] * r3[k,1]
          r3[i,2] -= SA[i,k,iz,ID] * r3[k,2]
        end
      end
      @unroll for i = 1 : M
        @unroll for j = 1 : M - 2
          r1[i,1] += -A13[i,j,iz,ID] * r3[j,1]
          r1[i,2] += -A13[i,j,iz,ID] * r3[j,2]
        end
        r1[i,1] *= invfac
        r1[i,2] *= invfac
      end
      @unroll for i = 1 : M - 2
        r2[i,1] = FT(0) 
        r2[i,2] = FT(0) 
        @unroll for j = 1 : M - 2
          r2[i,1] += -A23[i,j,iz,ID] * r3[j,1]
          r2[i,2] += -A23[i,j,iz,ID] * r3[j,2]
        end
        r2[i,1] *= invfac
        r2[i,2] *= invfac
      end

      s[1,1] = FT(0)
      s[2,1] = FT(0)
      s[3,1] = FT(0)
      s[4,1] = FT(0)
      s[1,2] = FT(0)
      s[2,2] = FT(0)
      s[3,2] = FT(0)
      s[4,2] = FT(0)
      @unroll for j = 1 : M - 2
        s[1,1] += - C14_3[1,j,iz,ID] * r3[j,1]
        s[1,2] += - C14_3[1,j,iz,ID] * r3[j,2]
        s[4,1] += - C14_3[2,j,iz,ID] * r3[j,1]
        s[4,2] += - C14_3[2,j,iz,ID] * r3[j,2]
        s[2,1] += - C23_2[1,j,iz,ID] * r2[j,1]
        s[2,2] += - C23_2[1,j,iz,ID] * r2[j,2]
        s[3,1] += - C23_2[2,j,iz,ID] * r2[j,1]
        s[3,2] += - C23_2[2,j,iz,ID] * r2[j,2]
      end
      @unroll for j = 1 : M
        s[2,1] += -C23_1[1,j,iz,ID] * r1[j,1]
        s[2,2] += -C23_1[1,j,iz,ID] * r1[j,2]
        s[3,1] += -C23_1[2,j,iz,ID] * r1[j,1]
        s[3,2] += -C23_1[2,j,iz,ID] * r1[j,2]
      end
      DU[ID,1,1,iz] += s[1,1]
      DU[ID,2,1,iz] += s[2,1]
      DU[ID,3,1,iz] += s[3,1]
      DU[ID,4,1,iz] += s[4,1]
      DU[ID,1,2,iz] += s[1,2]
      DU[ID,2,2,iz] += s[2,2]
      DU[ID,3,2,iz] += s[3,2]
      DU[ID,4,2,iz] += s[4,2]
    end
  end
end      
