@kernel inbounds = true function RiemanNonLinVertKernel!(RiemannSolver!,F,
  @Const(U),@Const(Aux),@Const(w), ::Val{M}, ::Val{NUMV}, ::Val{NAUX}) where {M, NUMV, NAUX}

  Iz, = @index(Global, NTuple)

  Nz = @uniform @ndrange()[1]

  VLL = @private eltype(F) (NUMV,)
  VRR = @private eltype(F) (NUMV,)
  AuxL = @private eltype(F) (NAUX,)
  AuxR = @private eltype(F) (NAUX,)
  FLoc = @private eltype(F) (NUMV,)

  RhoPos = @uniform 1
  wPos = @uniform 2
  ThPos = @uniform 3

  if Iz <= Nz
    if Iz > 1    
      @unroll for iAux = 1 : NAUX  
        AuxL[iAux] = Aux[M,Iz-1,iAux]
      end  
      @unroll for iv = 1 : NUMV  
        VLL[iv] = U[M,Iz-1,iv]
      end  
    else
      @unroll for iAux = 1 : NAUX  
        AuxL[iAux] = Aux[1,Iz,iAux]
      end  
      @unroll for iv = 1 : NUMV  
        VLL[iv] = U[1,Iz,iv]
      end  
      VLL[wPos] = -VLL[wPos]
    end  
    if Iz < Nz
      @unroll for iAux = 1 : NAUX  
        AuxR[iAux] = Aux[1,Iz,iAux]
      end  
      @unroll for iv = 1 : NUMV  
        VRR[iv] = U[1,Iz,iv]
      end  
    else  
      @unroll for iAux = 1 : NAUX  
        AuxR[iAux] = Aux[M,Iz-1,iAux]
      end  
      @unroll for iv = 1 : NUMV  
        VRR[iv] = U[M,Iz-1,iv]
      end  
      VRR[wPos] = -VRR[wPos]
    end

    @views RiemannSolver!(FLoc,VLL,VRR,AuxL,AuxR)

    Surf = eltype(U)(1) / w[1]  
    @unroll for iv = 1 : NUMV  
      FLoc[iv] *= Surf
    end  
    if Iz > 1 
      @unroll for iv = 1 : NUMV  
        @atomic :monotonic F[M,Iz-1,iv] -= FLoc[iv]
      end  
    end  
    if Iz < Nz
      @unroll for iv = 1 : NUMV  
        @atomic :monotonic F[1,Iz,iv] += FLoc[iv]
      end  
    end  
  end  
end

