
@kernel inbounds = true function LimitKernel!(DoF,qMin,qMax,@Const(Rhoq),@Const(Rho),@Const(Glob))

  iz = @index(Local, NTuple)
  Iz,IF,IT = @index(Global, NTuple)

  Nz = @uniform @ndrange()[1]
  NF = @uniform @ndrange()[2]
  NT = @uniform @ndrange()[3]


  qMin[1,Iz,IF,IT] = eltype(Rhoq)(1/0)
  qMax[1,Iz,IF,IT] = eltype(Rhoq)(-1/0)

  if Iz <= Nz && IF <= NF && IT <= NT
    for ID = 1 : DoF
      ind = Glob[ID,IF]
      qMin[1,Iz,IF,IT] = min(qMin[1,Iz,IF,IT],Rhoq[1,Iz,ind,IT] / Rho[1,Iz,ind])
      qMax[1,Iz,IF,IT] = max(qMax[1,Iz,IF,IT],Rhoq[1,Iz,ind,IT] / Rho[1,Iz,ind])
    end
  end
end  

@kernel inbounds = true function DivRhoTrUpwind3LimKernel!(FTr,@Const(Tr),@Const(U),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob),dt,@Const(w),@Const(qMin),@Const(qMax),@Const(Stencil))

# gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  @uniform l0  = eltype(FTr)(0)
  @uniform eta = eltype(FTr)(1.e-12)
  @uniform dlFD = eltype(FTr)(1.e-8)


  cCol = @localmem eltype(FTr) (N,N, ColumnTilesDim+3)
  uConCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  vConCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  DivRhoTr = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  DivRho = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  RhoTrColS = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  RhoColS = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  q = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  resp = @localmem eltype(FTr) (ColumnTilesDim)
  resc = @localmem eltype(FTr) (ColumnTilesDim)
  alpha = @localmem eltype(FTr) (ColumnTilesDim)
  lp = @localmem eltype(FTr) (ColumnTilesDim)
  lc = @localmem eltype(FTr) (ColumnTilesDim)
  sumJ = @localmem eltype(FTr) (ColumnTilesDim)
  qMinS = @localmem eltype(FTr) (ColumnTilesDim)
  qMaxS = @localmem eltype(FTr) (ColumnTilesDim)
  conv = @localmem (Bool) (ColumnTilesDim)
  if Iz <= Nz
    ID = I + (J - 1) * N  
    ind = Glob[ID,IF]
    cCol[I,J,iz+1] = Tr[1,Iz,ind] / U[1,Iz,ind,1]
    @views (uCon, vCon) = Contra12(-U[1,Iz,ind,1],U[1,Iz,ind,2],U[1,Iz,ind,3],dXdxI[1:2,1:2,:,ID,Iz,IF])
    uConCol[I,J,iz] = uCon
    vConCol[I,J,iz] = vCon
    if ID == 1
      resp[iz] = eltype(FTr)(0)  
      resc[iz] = eltype(FTr)(0)  
      sumJ[iz] = eltype(FTr)(0)  
      conv[iz] = true
      qMinS[iz] = qMin[1,Iz,Stencil[IF,1]]
      qMaxS[iz] = qMax[1,Iz,Stencil[IF,1]]
      for iS = 2 : 13
        qMinS[iz] = min(qMin[1,Iz,Stencil[IF,iS]],qMinS[iz])
        qMaxS[iz] = max(qMax[1,Iz,Stencil[IF,iS]],qMaxS[iz])
      end
    end  
  end
  if iz == 1
    Izm1 = max(Iz - 1,1)
    cCol[I,J,iz] = Tr[1,Izm1,ind] / U[1,Izm1,ind,1]
  end
  if iz == ColumnTilesDim || Iz == Nz
    Izp1 = min(Iz + 1,Nz)
    cCol[I,J,iz+2] = Tr[1,Izp1,ind] / U[1,Izp1,ind,1]
    Izp2 = min(Iz + 2,Nz)
    cCol[I,J,iz+3] = Tr[1,Izp2,ind] / U[1,Izp2,ind,1]
  end
  @synchronize

  if Iz <= Nz 
    ID = I + (J - 1) * N  
    @atomic :monotonic sumJ[iz] += JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]  
  end
  @synchronize

  if Iz < Nz 
    ID = I + (J - 1) * N  
    ind = Glob[ID,IF]
    cLL = cCol[I,J,iz]
    cL = cCol[I,J,iz+1]
    cR = cCol[I,J,iz+2]
    cRR = cCol[I,J,iz+3]

    @views wCon = Contra3(U[1,Iz:Iz+1,ind,1],U[1,Iz:Iz+1,ind,2],U[1,Iz:Iz+1,ind,3],
      U[1,Iz,ind,4],dXdxI[3,:,:,ID,Iz:Iz+1,IF])

    Izm1 = max(Iz - 1,1)
    Izp2 = min(Iz + 2, Nz)
    JLL = JJ[ID,1,Izm1,IF] + JJ[ID,2,Izm1,IF]
    JL = JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]
    JR = JJ[ID,1,Iz+1,IF] + JJ[ID,2,Iz+1,IF]
    JRR = JJ[ID,1,Izp2,IF] + JJ[ID,2,Izp2,IF]
    cFL, cFR = RecU4(cLL,cL,cR,cRR,JLL,JL,JR,JRR) 
    Flux = eltype(FTr)(0.25) * ((abs(wCon) + wCon) * cFL + (-abs(wCon) + wCon) * cFR)
    @atomic :monotonic FTr[1,Iz,ind] += -Flux / (M[Iz,ind,1] + M[Iz,ind,2])
    @atomic :monotonic FTr[1,Iz+1,ind] += Flux / (M[Iz+1,ind,1] + M[Iz+1,ind,2])
  end 

  if Iz <= Nz
    ID = I + (J - 1) * N  
    DivRhoTr[I,J,iz] = D[I,1] * uConCol[1,J,iz] * cCol[1,J,iz+1] 
    DivRhoTr[I,J,iz] += D[J,1] * vConCol[I,1,iz] * cCol[I,1,iz+1]
    DivRho[I,J,iz] = D[I,1] * uConCol[1,J,iz]
    DivRho[I,J,iz] += D[J,1] * vConCol[I,1,iz]
    for k = 2 : N
      DivRhoTr[I,J,iz] += D[I,k] * uConCol[k,J,iz] * cCol[k,J,iz+1] 
      DivRhoTr[I,J,iz] += D[J,k] * vConCol[I,k,iz] * cCol[I,k,iz+1]
      DivRho[I,J,iz] += D[I,k] * uConCol[k,J,iz]
      DivRho[I,J,iz] += D[J,k] * vConCol[I,k,iz]
    end
    ind = Glob[ID,IF]
    RhoTrColS[I,J,iz] = Tr[1,Iz,ind] + dt * DivRhoTr[I,J,iz] / (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF])
    RhoColS[I,J,iz] = U[1,Iz,ind,1] + dt * DivRho[I,J,iz] / (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF])
    #   Finite difference step
    q[I,J,iz] = medianGPU(qMinS[iz], RhoTrColS[I,J,iz] / RhoColS[I,J,iz] +
      l0,  qMaxS[iz])
    @atomic :monotonic resp[iz] += (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]) * w[I] * w[J] / sumJ[iz] * 
      (q[I,J,iz] * RhoColS[I,J,iz] - RhoTrColS[I,J,iz])
  end
  @synchronize
  if Iz <= Nz
    ID = I + (J - 1) * N  
    if abs(resp[iz]) <= eta 
      if ID == 1
        conv[iz] = false
      end  
    else
      qLoc = medianGPU(qMinS[iz],  RhoTrColS[I,J,iz] / RhoColS[I,J,iz] + 
        (l0 + dlFD),  qMaxS[iz])
      @atomic :monotonic resc[iz] += (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]) * w[I] * w[J] / sumJ[iz] * 
        (qLoc * RhoColS[I,J,iz]  - RhoTrColS[I,J,iz])  
    end
  end
  @synchronize

  if Iz <= Nz && I == 1 && J == 1 && conv[iz]
    if abs(resc[iz] - resp[iz]) <= eltype(FTr)(1.e-13)
      conv[iz] = false
    else
      alpha[iz] = dlFD / (resc[iz] - resp[iz])
      lp[iz] = l0
      lc[iz] = lp[iz] - alpha[iz] * resp[iz]
      resp[iz] = eltype(FTr)(0)
      resc[iz] = eltype(FTr)(0)
    end
  end  
  @synchronize
  for iTer = 1 : 8
    if Iz <= Nz && conv[iz]
      ID = I + (J - 1) * N  
      q[I,J,iz] = medianGPU(qMinS[iz], RhoTrColS[I,J,iz] / RhoColS[I,J,iz] +
        lc[iz],  qMaxS[iz])
      @atomic :monotonic resc[iz] += (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]) * w[I] * w[J] / sumJ[iz] *
        (q[I,J,iz] * RhoColS[I,J,iz] - RhoTrColS[I,J,iz])
    end  
    @synchronize
    if Iz <= Nz && I == 1 && J == 1 && conv[iz]
      if abs(resc[iz] - resp[iz]) <= eltype(FTr)(1.e-13) 
        conv[iz] = false
      else  
        alpha[iz] = (lp[iz] - lc[iz]) / (resp[iz] - resc[iz])
        resp[iz] = resc[iz]
        lp[iz] = lc[iz]
        lc[iz] = lc[iz] - alpha[iz] * resc[iz]
        resc[iz] = eltype(FTr)(0)  
      end  
    end  
    @synchronize
  end
  if Iz <= Nz
    ID = I + (J - 1) * N  
    ind = Glob[ID,IF]
    @atomic :monotonic FTr[1,Iz,ind] += (q[I,J,iz] * RhoColS[I,J,iz] - Tr[1,Iz,ind]) *
      (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]) / dt / (M[Iz,ind,1] + M[Iz,ind,2]) 
  end  
end

# Clip and assured sum (CAAS, limiter_option 9 in HOMME): clip the predictor to the
# stencil bounds and give the clipped mass back in proportion to the remaining capacity.
# The element sums are formed by every node from shared memory, no shared atomics.
@kernel inbounds = true function DivRhoTrUpwind3CAASKernel!(FTr,@Const(Tr),@Const(U),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob),dt,@Const(w),@Const(qMin),@Const(qMax),@Const(Stencil))

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NS = @uniform size(Stencil,2)

  cCol = @localmem eltype(FTr) (N,N, ColumnTilesDim+3)
  uConCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  vConCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  BoundCol = @localmem eltype(FTr) (N,N,2, ColumnTilesDim)
  MassCol = @localmem eltype(FTr) (N,N,3, ColumnTilesDim)
  Loc = @private eltype(FTr) (4)

  if Iz <= Nz
    ID = I + (J - 1) * N  
    ind = Glob[ID,IF]
    cCol[I,J,iz+1] = Tr[1,Iz,ind] / U[1,Iz,ind,1]
    @views (uCon, vCon) = Contra12(-U[1,Iz,ind,1],U[1,Iz,ind,2],U[1,Iz,ind,3],dXdxI[1:2,1:2,:,ID,Iz,IF])
    uConCol[I,J,iz] = uCon
    vConCol[I,J,iz] = vCon
#   Bounds over the stencil elements, distributed over the nodes
    qMinLoc = eltype(FTr)(1/0)
    qMaxLoc = eltype(FTr)(-1/0)
    for iS = ID : N * N : NS
      qMinLoc = min(qMinLoc,qMin[1,Iz,Stencil[IF,iS]])
      qMaxLoc = max(qMaxLoc,qMax[1,Iz,Stencil[IF,iS]])
    end
    BoundCol[I,J,1,iz] = qMinLoc
    BoundCol[I,J,2,iz] = qMaxLoc
    if iz == 1
      Izm1 = max(Iz - 1,1)
      cCol[I,J,iz] = Tr[1,Izm1,ind] / U[1,Izm1,ind,1]
    end
    if iz == ColumnTilesDim || Iz == Nz
      Izp1 = min(Iz + 1,Nz)
      cCol[I,J,iz+2] = Tr[1,Izp1,ind] / U[1,Izp1,ind,1]
      Izp2 = min(Iz + 2,Nz)
      cCol[I,J,iz+3] = Tr[1,Izp2,ind] / U[1,Izp2,ind,1]
    end
  end
  @synchronize

  if Iz < Nz 
    ID = I + (J - 1) * N  
    ind = Glob[ID,IF]
    cLL = cCol[I,J,iz]
    cL = cCol[I,J,iz+1]
    cR = cCol[I,J,iz+2]
    cRR = cCol[I,J,iz+3]

    @views wCon = Contra3(U[1,Iz:Iz+1,ind,1],U[1,Iz:Iz+1,ind,2],U[1,Iz:Iz+1,ind,3],
      U[1,Iz,ind,4],dXdxI[3,:,:,ID,Iz:Iz+1,IF])

    Izm1 = max(Iz - 1,1)
    Izp2 = min(Iz + 2, Nz)
    JLL = JJ[ID,1,Izm1,IF] + JJ[ID,2,Izm1,IF]
    JL = JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]
    JR = JJ[ID,1,Iz+1,IF] + JJ[ID,2,Iz+1,IF]
    JRR = JJ[ID,1,Izp2,IF] + JJ[ID,2,Izp2,IF]
    cFL, cFR = RecU4(cLL,cL,cR,cRR,JLL,JL,JR,JRR) 
    Flux = eltype(FTr)(0.25) * ((abs(wCon) + wCon) * cFL + (-abs(wCon) + wCon) * cFR)
    @atomic :monotonic FTr[1,Iz,ind] += -Flux / (M[Iz,ind,1] + M[Iz,ind,2])
    @atomic :monotonic FTr[1,Iz+1,ind] += Flux / (M[Iz+1,ind,1] + M[Iz+1,ind,2])
  end 

  if Iz <= Nz
    ID = I + (J - 1) * N  
    ind = Glob[ID,IF]
    DivRhoTr = D[I,1] * uConCol[1,J,iz] * cCol[1,J,iz+1] + D[J,1] * vConCol[I,1,iz] * cCol[I,1,iz+1]
    DivRho = D[I,1] * uConCol[1,J,iz] + D[J,1] * vConCol[I,1,iz]
    for k = 2 : N
      DivRhoTr += D[I,k] * uConCol[k,J,iz] * cCol[k,J,iz+1] + D[J,k] * vConCol[I,k,iz] * cCol[I,k,iz+1]
      DivRho += D[I,k] * uConCol[k,J,iz] + D[J,k] * vConCol[I,k,iz]
    end
    JLoc = JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]
    RhoTrS = Tr[1,Iz,ind] + dt * DivRhoTr / JLoc
    RhoS = U[1,Iz,ind,1] + dt * DivRho / JLoc
    qMinS = BoundCol[1,1,1,iz]
    qMaxS = BoundCol[1,1,2,iz]
    for j = 1 : N, i = 1 : N
      qMinS = min(qMinS,BoundCol[i,j,1,iz])
      qMaxS = max(qMaxS,BoundCol[i,j,2,iz])
    end
    qS = RhoTrS / RhoS
    qC = medianGPU(qMinS, qS, qMaxS)
    MassW = JLoc * w[I] * w[J] * RhoS
    MassCol[I,J,1,iz] = MassW * (qS - qC)
    MassCol[I,J,2,iz] = MassW * (qMaxS - qC)
    MassCol[I,J,3,iz] = MassW * (qC - qMinS)
    Loc[1] = qC
    Loc[2] = RhoS
    Loc[3] = qMinS
    Loc[4] = qMaxS
  end
  @synchronize

  if Iz <= Nz
    ID = I + (J - 1) * N  
    ind = Glob[ID,IF]
    dm = eltype(FTr)(0)
    CapU = eltype(FTr)(0)
    CapL = eltype(FTr)(0)
    for j = 1 : N, i = 1 : N
      dm += MassCol[i,j,1,iz]
      CapU += MassCol[i,j,2,iz]
      CapL += MassCol[i,j,3,iz]
    end
    qC = Loc[1]
    if dm > eltype(FTr)(0) && CapU > eltype(FTr)(0)
      qC += dm / CapU * (Loc[4] - qC)
    elseif dm < eltype(FTr)(0) && CapL > eltype(FTr)(0)
      qC += dm / CapL * (qC - Loc[3])
    end
    @atomic :monotonic FTr[1,Iz,ind] += (qC * Loc[2] - Tr[1,Iz,ind]) *
      (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]) / dt / (M[Iz,ind,1] + M[Iz,ind,2]) 
  end  
end

# ρ and ρθ as in DivRhoTrUpwind3New2Kernel!, all tracers with the CAAS limiter
@kernel inbounds = true function DivRhoTrUpwind3New2CAASKernel!(F,NumV,NumTr,@Const(U),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob),dt,@Const(w),@Const(qB))

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  cCol = @localmem eltype(F) (N,N, ColumnTilesDim+3)
  RhoCol = @localmem eltype(F) (N,N, ColumnTilesDim+1)
  uCol = @localmem eltype(F) (N,N, ColumnTilesDim+1)
  vCol = @localmem eltype(F) (N,N, ColumnTilesDim+1)
  wCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  JCol = @localmem eltype(F) (N,N, ColumnTilesDim+3)
  MCCol = @localmem eltype(F) (N,N, ColumnTilesDim+1)
  MassCol = @localmem eltype(F) (N,N,3, ColumnTilesDim)
  Loc = @private eltype(F) (5)
  if Iz <= Nz
    cCol[I,J,iz+1] = U[1,Iz,ind,5] / U[1,Iz,ind,1]
    JCol[I,J,iz+1] = JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]
    RhoCol[I,J,iz] = U[1,Iz,ind,1]
    uCol[I,J,iz] = U[1,Iz,ind,2]
    vCol[I,J,iz] = U[1,Iz,ind,3]
    wCol[I,J,iz] = U[1,Iz,ind,4]
    MCCol[I,J,iz] = M[Iz,ind,1] + M[Iz,ind,2]
  end
  if iz == 1
    Izm1 = max(Iz - 1,1)
    cCol[I,J,iz] = U[1,Izm1,ind,5] / U[1,Izm1,ind,1]
    JCol[I,J,iz] = JJ[ID,1,Izm1,IF] + JJ[ID,2,Izm1,IF] 
  end
  if iz == ColumnTilesDim || Iz == Nz
    Izp1 = min(Iz + 1,Nz)
    cCol[I,J,iz+2] = U[1,Izp1,ind,5] / U[1,Izp1,ind,1]
    JCol[I,J,iz+2] = JJ[ID,1,Izp1,IF] + JJ[ID,2,Izp1,IF] 
    Izp2 = min(Iz + 2,Nz)
    cCol[I,J,iz+3] = U[1,Izp2,ind,5] / U[1,Izp2,ind,1]
    JCol[I,J,iz+3] = JJ[ID,1,Izp2,IF] + JJ[ID,2,Izp2,IF] 
    RhoCol[I,J,iz+1] = U[1,Izp1,ind,1]
    uCol[I,J,iz+1] = U[1,Izp1,ind,2]
    vCol[I,J,iz+1] = U[1,Izp1,ind,3]
    MCCol[I,J,iz+1] = M[Izp1,ind,1] + M[Izp1,ind,2]
  end
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz < Nz 
    @views wCon = Contra3(RhoCol[I,J,iz:iz+1],uCol[I,J,iz:iz+1],vCol[I,J,iz:iz+1],
      wCol[I,J,iz],dXdxI[3,:,:,ID,Iz:Iz+1,IF])
    wCol[I,J,iz] = wCon
    cFL, cFR = RecU4(cCol[I,J,iz],cCol[I,J,iz+1],cCol[I,J,iz+2],cCol[I,J,iz+3],
      JCol[I,J,iz],JCol[I,J,iz+1],JCol[I,J,iz+2],JCol[I,J,iz+3]) 
#   central second order at the interfaces next to the boundaries, no ghost cells
    if Iz == 1 || Iz == Nz - 1
      cFL = eltype(F)(0.5) * (cCol[I,J,iz+1] + cCol[I,J,iz+2])
      cFR = cFL
    end
    Flux = eltype(F)(0.25) * ((abs(wCon) + wCon) * cFL + 
      (-abs(wCon) + wCon) * cFR)
    @atomic :monotonic F[1,Iz,ind,5] += -Flux / MCCol[I,J,iz]
    @atomic :monotonic F[1,Iz+1,ind,5] += Flux / MCCol[I,J,iz+1]
    Flux = eltype(F)(0.5)*wCon
    @atomic :monotonic F[1,Iz,ind,1] += -Flux / MCCol[I,J,iz]
    @atomic :monotonic F[1,Iz+1,ind,1] += Flux / MCCol[I,J,iz+1]
  end 

  if Iz <= Nz
    uCol[I,J,iz], vCol[I,J,iz] = Contra12(-RhoCol[I,J,iz],uCol[I,J,iz],vCol[I,J,iz],view(dXdxI,1:2,1:2,:,ID,Iz,IF))
  end  
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]
  if Iz <= Nz
    DivRhoTr = D[I,1] * uCol[1,J,iz] * cCol[1,J,iz+1] + D[J,1] * vCol[I,1,iz] * cCol[I,1,iz+1]
    DivTr1 = D[I,1] * cCol[1,J,iz+1] 
    DivTr2 = D[J,1] * cCol[I,1,iz+1]
    DivRho = D[I,1] * uCol[1,J,iz] + D[J,1] * vCol[I,1,iz] 
    for k = 2 : N
      DivRhoTr += D[I,k] * uCol[k,J,iz] * cCol[k,J,iz+1] + D[J,k] * vCol[I,k,iz] * cCol[I,k,iz+1]
      DivTr1 += D[I,k] * cCol[k,J,iz+1]
      DivTr2 += D[J,k] * cCol[I,k,iz+1]
      DivRho += D[I,k] * uCol[k,J,iz] + D[J,k] * vCol[I,k,iz]
    end
    @atomic :monotonic F[1,Iz,ind,5] += eltype(F)(0.5) * (DivRhoTr + uCol[I,J,iz] * DivTr1 + vCol[I,J,iz] * DivTr2 +
      cCol[I,J,iz+1] * DivRho)/ MCCol[I,J,iz]
    @atomic :monotonic F[1,Iz,ind,1] += DivRho / MCCol[I,J,iz]
#   density of the predictor for the tracer limiter
    Loc[5] = RhoCol[I,J,iz] + dt * DivRho / JCol[I,J,iz+1]
  end

  for iT = NumV + 1 : NumV + NumTr
#   all threads have finished reading cCol and MassCol of the previous tracer
    @synchronize
    ID = I + (J - 1) * N  
    ind = Glob[ID,IF]
    if Iz <= Nz
      cCol[I,J,iz+1] = U[1,Iz,ind,iT] / U[1,Iz,ind,1]
    end
    if iz == 1
      Izm1 = max(Iz - 1,1)
      cCol[I,J,iz] = U[1,Izm1,ind,iT] / U[1,Izm1,ind,1]
    end
    if iz == ColumnTilesDim || Iz == Nz
      Izp1 = min(Iz + 1,Nz)
      cCol[I,J,iz+2] = U[1,Izp1,ind,iT] / U[1,Izp1,ind,1]
      Izp2 = min(Iz + 2,Nz)
      cCol[I,J,iz+3] = U[1,Izp2,ind,iT] / U[1,Izp2,ind,1]
    end
    @synchronize

    ID = I + (J - 1) * N
    ind = Glob[ID,IF]
    if Iz < Nz
      wCon = wCol[I,J,iz]
      cFL, cFR = RecU4(cCol[I,J,iz],cCol[I,J,iz+1],cCol[I,J,iz+2],cCol[I,J,iz+3],
        JCol[I,J,iz],JCol[I,J,iz+1],JCol[I,J,iz+2],JCol[I,J,iz+3]) 
      Flux = eltype(F)(0.25) * ((abs(wCon) + wCon) * cFL + 
        (-abs(wCon) + wCon) * cFR)
      @atomic :monotonic F[1,Iz,ind,iT] += -Flux / MCCol[I,J,iz]
      @atomic :monotonic F[1,Iz+1,ind,iT] += Flux / MCCol[I,J,iz+1]
    end
#   horizontal predictor, clipped to the stencil bounds (CAAS)
    if Iz <= Nz
      DivRhoTr = D[I,1] * uCol[1,J,iz] * cCol[1,J,iz+1] + D[J,1] * vCol[I,1,iz] * cCol[I,1,iz+1]
      for k = 2 : N
        DivRhoTr += D[I,k] * uCol[k,J,iz] * cCol[k,J,iz+1] + D[J,k] * vCol[I,k,iz] * cCol[I,k,iz+1]
      end
      RhoS = Loc[5]
      qS = (U[1,Iz,ind,iT] + dt * DivRhoTr / JCol[I,J,iz+1]) / RhoS
      qMinS = qB[1,Iz,IF,iT-NumV]
      qMaxS = qB[1,Iz,IF,iT-NumV+NumTr]
      qC = medianGPU(qMinS, qS, qMaxS)
      MassW = JCol[I,J,iz+1] * w[I] * w[J] * RhoS
      MassCol[I,J,1,iz] = MassW * (qS - qC)
      MassCol[I,J,2,iz] = MassW * (qMaxS - qC)
      MassCol[I,J,3,iz] = MassW * (qC - qMinS)
      Loc[1] = qC
      Loc[3] = qMinS
      Loc[4] = qMaxS
    end
    @synchronize

    ID = I + (J - 1) * N
    ind = Glob[ID,IF]
    if Iz <= Nz
      dm = eltype(F)(0)
      CapU = eltype(F)(0)
      CapL = eltype(F)(0)
      for j = 1 : N, i = 1 : N
        dm += MassCol[i,j,1,iz]
        CapU += MassCol[i,j,2,iz]
        CapL += MassCol[i,j,3,iz]
      end
      qC = Loc[1]
      if dm > eltype(F)(0) && CapU > eltype(F)(0)
        qC += dm / CapU * (Loc[4] - qC)
      elseif dm < eltype(F)(0) && CapL > eltype(F)(0)
        qC += dm / CapL * (qC - Loc[3])
      end
      @atomic :monotonic F[1,Iz,ind,iT] += (qC * Loc[5] - U[1,Iz,ind,iT]) * JCol[I,J,iz+1] / dt / MCCol[I,J,iz]
    end
  end
end

# Bounds of the limiter over the stencil elements, from the element minima and maxima in q
@kernel inbounds = true function StencilBoundKernel!(qB,@Const(q),@Const(Stencil),NumTr)

  Iz,IF,IT = @index(Global, NTuple)

  Nz = @uniform @ndrange()[1]
  NF = @uniform @ndrange()[2]
  NT = @uniform @ndrange()[3]

  if Iz <= Nz && IF <= NF && IT <= NT
    qMinLoc = q[1,Iz,Stencil[IF,1],IT]
    qMaxLoc = q[1,Iz,Stencil[IF,1],NumTr+IT]
    for iS = 2 : size(Stencil,2)
      qMinLoc = min(qMinLoc,q[1,Iz,Stencil[IF,iS],IT])
      qMaxLoc = max(qMaxLoc,q[1,Iz,Stencil[IF,iS],NumTr+IT])
    end
    qB[1,Iz,IF,IT] = qMinLoc
    qB[1,Iz,IF,NumTr+IT] = qMaxLoc
  end
end

@kernel inbounds = true function DivRhoTrViscUpwind3LimKernel!(FTr,@Const(Tr),@Const(U),@Const(Cache),@Const(D),@Const(DW),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob),Koeff,dt,@Const(w),@Const(qMin),@Const(qMax),@Const(Stencil))

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  @uniform l0  = eltype(FTr)(0)
  @uniform eta = eltype(FTr)(1.e-12)
  @uniform dlFD = eltype(FTr)(1.e-8)

  cCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  CacheCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  RhoCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  uCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  vCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  wCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  DivRhoTr = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  DivRho = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  RhoTrColS = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  RhoColS = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  q = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  resp = @localmem eltype(FTr) (ColumnTilesDim)
  resc = @localmem eltype(FTr) (ColumnTilesDim)
  alpha = @localmem eltype(FTr) (ColumnTilesDim)
  lp = @localmem eltype(FTr) (ColumnTilesDim)
  lc = @localmem eltype(FTr) (ColumnTilesDim)
  sumJ = @localmem eltype(FTr) (ColumnTilesDim)
  qMinS = @localmem eltype(FTr) (ColumnTilesDim)
  qMaxS = @localmem eltype(FTr) (ColumnTilesDim)
  conv = @localmem (Bool) (ColumnTilesDim)
  if Iz <= Nz
    ID = I + (J - 1) * N  
    ind = Glob[ID,IF]
    CacheCol[I,J,iz] = Cache[Iz,ind]
    wCol[I,J,iz] = U[1,Iz,ind,4]
    RhoCol[I,J,iz] = U[1,Iz,ind,1]
    cCol[I,J,iz] = Tr[Iz,ind] / RhoCol[I,J,iz]
    uCol[I,J,iz] = U[1,Iz,ind,2]
    vCol[I,J,iz] = U[1,Iz,ind,3]
    DivRho[I,J,iz] = eltype(FTr)(0)
    DivRhoTr[I,J,iz] = eltype(FTr)(0)
    if ID == 1
      resp[iz] = eltype(FTr)(0)
      resc[iz] = eltype(FTr)(0)
      sumJ[iz] = eltype(FTr)(0)
      conv[iz] = true
      qMinS[iz] = minimum(qMin[Iz,Stencil[IF,:]])
      qMaxS[iz] = maximum(qMax[Iz,Stencil[IF,:]])
    end
  end
  @synchronize
  if Iz < Nz 
    ID = I + (J - 1) * N  
    ind = Glob[ID,IF]
    ind = Glob[ID,IF]
    cL = cCol[I,J,iz]
    cR = cCol[I,J,iz+1]
    if iz > 1
      cLL = cCol[I,J,iz-1]
    else
      Izm1 = max(Iz - 1,1)
      cLL = U[1,Izm1,ind,5] / U[1,Izm1,ind,1]
    end
    if iz < ColumnTilesDim - 1
      cRR = cCol[I,J,iz+2]
    else
      Izp2 = min(Iz + 2, Nz)
      cRR = U[1,Izp2,ind,5] / U[1,Izp2,ind,1]
    end

    @views wCon = Contra3(U[1,Iz:Iz+1,ind,1],U[1,Iz:Iz+1,ind,2],U[1,Iz:Iz+1,ind,3],
      U[1,Iz,ind,4],dXdxI[3,:,:,ID,Iz:Iz+1,IF])

    Izm1 = max(Iz - 1,1)
    Izp2 = min(Iz + 2, Nz)
    JLL = JJ[ID,1,Izm1,IF] + JJ[ID,2,Izm1,IF]
    JL = JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]
    JR = JJ[ID,1,Iz+1,IF] + JJ[ID,2,Iz+1,IF]
    JRR = JJ[ID,1,Izp2,IF] + JJ[ID,2,Izp2,IF]
    cFL, cFR = RecU4(cLL,cL,cR,cRR,JLL,JL,JR,JRR) 
    Flux = 0.25 * ((abs(wCon) + wCon) * cFL + (-abs(wCon) + wCon) * cFR)
    @atomic :monotonic FTr[Iz,ind] += -Flux / M[Iz,ind]
    @atomic :monotonic FTr[Iz+1,ind] += Flux / M[Iz+1,ind]
  end

  if Iz <= Nz
    ID = I + (J - 1) * N  
    Dxc = 0
    Dyc = 0
    for k = 1 : N
      Dxc = Dxc + D[I,k] * CacheCol[k,J,iz]
      Dyc = Dyc + D[J,k] * CacheCol[I,k,iz]
    end
    
    @views (GradDx, GradDy) = Grad12(RhoCol[I,J,iz],Dxc,Dyc,dXdxI[1:2,1:2,:,ID,Iz,IF],JJ[ID,:,Iz,IF])
    @views (tempx, tempy) = Contra12(-Koeff,GradDx,GradDy,dXdxI[1:2,1:2,:,ID,Iz,IF])
    for k = 1 : N
      @atomic :monotonic DivRhoTr[k,J,iz] += DW[k,I] * tempx
      @atomic :monotonic DivRhoTr[I,k,iz] += DW[k,J] * tempy
    end

    @views (tempxRho, tempyRho) = Contra12(-RhoCol[I,J,iz],uCol[I,J,iz],vCol[I,J,iz],dXdxI[1:2,1:2,:,ID,Iz,IF])
    for k = 1 : N
      @atomic :monotonic DivRho[k,J,iz] += D[k,I] * tempxRho
      @atomic :monotonic DivRho[I,k,iz] += D[k,J] * tempyRho
    end
    tempxTr = tempxRho * cCol[I,J,iz]
    tempyTr = tempyRho * cCol[I,J,iz]
    for k = 1 : N
      @atomic :monotonic DivRhoTr[k,J,iz] += D[k,I] * tempxTr
      @atomic :monotonic DivRhoTr[I,k,iz] += D[k,J] * tempyTr
    end
    @atomic :monotonic sumJ[iz] += JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]  
  end
  @synchronize

  if Iz <=Nz
    ID = I + (J - 1) * N  
    ind = Glob[ID,IF]  
    RhoTrColS[I,J,iz] = Tr[Iz,ind] + dt * DivRhoTr[I,J,iz] / (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF])
    RhoColS[I,J,iz] = U[1,Iz,ind,1] + dt * DivRho[I,J,iz] / (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF])
    #   Finite difference step
    q[I,J,iz] = medianGPU(qMinS[iz], RhoTrColS[I,J,iz] / RhoColS[I,J,iz] +
      l0,  qMaxS[iz])
    @atomic :monotonic resp[iz] += (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]) * w[I] * w[J] / sumJ[iz] * 
      (q[I,J,iz] * RhoColS[I,J,iz] - RhoTrColS[I,J,iz])
  end
  @synchronize
  if Iz <= Nz
    ID = I + (J - 1) * N  
    if abs(resp[iz]) <= eta 
      if ID == 1
        conv[iz] = false
      end  
    else
      qLoc = medianGPU(qMinS[iz],  RhoTrColS[I,J,iz] / RhoColS[I,J,iz] + 
        (l0 + dlFD), qMaxS[iz])
      @atomic :monotonic resc[iz] += (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]) * w[I] * w[J] / sumJ[iz] * 
        (qLoc * RhoColS[I,J,iz]  - RhoTrColS[I,J,iz])  
    end
  end
  @synchronize

  if Iz <= Nz && I == 1 && J == 1 && conv[iz]
    if abs(resc[iz] - resp[iz]) <= eltype(FTr)(1.e-13)
      conv[iz] = false
    else
      alpha[iz] = dlFD / (resc[iz] - resp[iz])
      lp[iz] = l0
      lc[iz] = lp[iz] - alpha[iz] * resp[iz]
      resp[iz] = eltype(FTr)(0)
      resc[iz] = eltype(FTr)(0)
    end
  end  
  @synchronize
  for iTer = 1 : 5
    if Iz <= Nz && conv[iz]
      ID = I + (J - 1) * N  
      q[I,J,iz] = medianGPU(qMinS[iz], RhoTrColS[I,J,iz] / RhoColS[I,J,iz] +
        lc[iz],  qMaxS[iz])
      @atomic :monotonic resc[iz] += (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]) * w[I] * w[J] / sumJ[iz] *
        (q[I,J,iz] * RhoColS[I,J,iz] - RhoTrColS[I,J,iz])
    end  
    @synchronize
    if Iz <= Nz && I == 1 && J == 1 && conv[iz]
      if abs(resc[iz] - resp[iz]) <= eltype(FTr)(1.e-13) 
        conv[iz] = false
      else  
        alpha[iz] = (lp[iz] - lc[iz]) / (resp[iz] - resc[iz])
        resp[iz] = resc[iz]
        lp[iz] = lc[iz]
        lc[iz] = lc[iz] - alpha[iz] * resc[iz]
        resc[iz] = eltype(FTr)(0)  
      end  
    end  
    @synchronize
  end
  if Iz <= Nz
    ID = I + (J - 1) * N  
    ind = Glob[ID,IF]
    @show Iz,ind,q[I,J,iz],RhoColS[I,J,iz]
    @atomic :monotonic FTr[Iz,ind] += (q[I,J,iz] * RhoColS[I,J,iz] - Tr[Iz,ind]) *
      (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]) / dt / M[Iz,ind]
  end  
end

@inline function medianGPU(a1,a2,a3)
  if a1 <= a2
    if a2 <= a3
      m = a2
    elseif a1 <= a3
      m = a3
    else
      m = a1
    end
  else
    if a1 <= a3 
      m = a1
    elseif a2 <= a3
      m = a3
    else
      m = a2
    end
  end
end

