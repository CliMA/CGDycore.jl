@kernel inbounds = true function DivRhoGradKernel!(F,@Const(U),@Const(D),@Const(DW),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob))

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N
  ind = Glob[ID,IF]

  cCol = @localmem eltype(F) (N+BANK,N, ColumnTilesDim)
  FCol = @localmem eltype(F) (N,N, ColumnTilesDim)

  if Iz <= Nz
    cCol[I,J,iz] = U[Iz,ind,5] / U[Iz,ind,1]
    FCol[I,J,iz] = 0.0
  end
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz <= Nz
    Dxc = D[I,1] * cCol[1,J,iz]
    Dyc = D[J,1] * cCol[I,1,iz]
    for k = 2 : N
      Dxc = Dxc + D[I,k] * cCol[k,J,iz]
      Dyc = Dyc + D[J,k] * cCol[I,k,iz] 
    end
    @views (GradDx, GradDy) = Grad12(Dxc,Dyc,dXdxI[1:2,1:2,:,ID,Iz,IF],JJ[ID,:,Iz,IF])
    @views (tempx, tempy) = Contra12(GradDx,GradDy,dXdxI[1:2,1:2,:,ID,Iz,IF])
    for k = 1 : N
      @atomic :monotonic FCol[k,J,iz] += DW[k,I] * tempx
      @atomic :monotonic FCol[I,k,iz] += DW[k,J] * tempy
    end
  end

  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]
  if Iz <= Nz
    @atomic :monotonic F[Iz,ind,5] += FCol[I,J,iz] / (M[Iz,ind,1] + M[Iz,ind,2])
  end
end


@kernel inbounds = true function DivRhoTrCentralKernel!(F,@Const(c),@Const(uC),@Const(vC),@Const(w),
  @Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob))

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  cCol = @localmem eltype(F) (N+BANK,N, ColumnTilesDim+2)
  uCol = @localmem eltype(F) (N+BANK,N, ColumnTilesDim+2)
  vCol = @localmem eltype(F) (N+BANK,N, ColumnTilesDim+2)
  wCol = @localmem eltype(F) (N+BANK,N, ColumnTilesDim+1)
  FCol = @localmem eltype(F) (N+BANK,N, ColumnTilesDim+2)
  if Iz <= Nz
    wCol[I,J,iz+1] = w[Iz,ind]
    cCol[I,J,iz+1] = c[Iz,ind]
    uCol[I,J,iz+1] = uC[Iz,ind]
    vCol[I,J,iz+1] = vC[Iz,ind]
    FCol[I,J,iz+1] = 0
    if iz == 1 && Iz > 1
      cCol[I,J,1] = c[Iz-1,ind]
      uCol[I,J,1] = uC[Iz-1,ind]
      vCol[I,J,1] = vC[Iz-1,ind]
      wCol[I,J,1] = w[Iz,ind]
      FCol[I,J,1] = 0
    elseif iz == 1 && Iz == 1
      cCol[I,J,1] = c[1,ind]
      wCol[I,J,1] = 0
      FCol[I,J,1] = 0
    end
    if iz == ColumnTilesDim && Iz < Nz
      cCol[I,J,ColumnTilesDim+2] = c[Iz+1,ind]
      uCol[I,J,ColumnTilesDim+2] = uC[Iz+1,ind]
      vCol[I,J,ColumnTilesDim+2] = vC[Iz+1,ind]
      FCol[I,J,ColumnTilesDim+2] = 0
    elseif iz == ColumnTilesDim && Iz == Nz
      cCol[I,J,ColumnTilesDim+2] = c[Nz,ind]
      FCol[I,J,ColumnTilesDim+2] = 0
    end
  end
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz < Nz 
    wCon = dXdxI[I,J,2,Iz,3,1,IF] * uCol[I,J,iz+1] + 
      dXdxI[I,J,2,Iz,3,2,IF] * vCol[I,J,iz+1] + 
      dXdxI[I,J,1,Iz+1,3,1,IF] * uCol[I,J,iz+2] + 
      dXdxI[I,J,2,Iz+1,3,2,IF] * vCol[I,J,iz+2] + 
      (dXdxI[I,J,2,Iz,3,3,IF] + dXdxI[I,J,1,Iz+1,3,3,IF]) * wCol[I,J,iz+1]
    cF = (JJ[ID,2,Iz,IF] * cCol[I,J,iz+1] + JJ[ID,1,Iz+1,IF] * cCol[I,J,iz+2]) /
      (JJ[ID,2,Iz,IF] + JJ[ID,1,Iz+1,IF])
    Flux = eltype(F)(0.5) * wCon * cF
    @atomic :monotonic FCol[I,J,iz+1] += -Flux
    @atomic :monotonic FCol[I,J,iz+2] += Flux
  end 

  if Iz <= Nz
    tempx = -cCol[I,J,iz+1] * ((dXdxI[I,J,1,Iz,1,1,IF] + dXdxI[I,J,2,Iz,1,1,IF]) * uCol[I,J,iz+1] +
      (dXdxI[I,J,1,Iz,1,2,IF] + dXdxI[I,J,2,Iz,1,2,IF]) * vCol[I,J,iz+1] +
      dXdxI[I,J,1,Iz,1,3,IF] * wCol[I,J,iz] + dXdxI[I,J,2,Iz,1,3,IF] * wCol[I,J,iz+1])
    tempy = -cCol[I,J,iz+1] * ((dXdxI[I,J,1,Iz,2,1,IF] + dXdxI[I,J,2,Iz,2,1,IF]) * uCol[I,J,iz+1] +
      (dXdxI[I,J,1,Iz,2,2,IF] + dXdxI[I,J,2,Iz,2,2,IF]) * vCol[I,J,iz+1] +
      dXdxI[I,J,1,Iz,2,3,IF] * wCol[I,J,iz] + dXdxI[I,J,2,Iz,2,3,IF] * wCol[I,J,iz+1])
    for k = 1 : N
      @atomic :monotonic FCol[k,J,iz+1] += D[k,I] * tempx
      @atomic :monotonic FCol[I,k,iz+1] += D[k,J] * tempy
    end
  end
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz <= Nz 
    @atomic :monotonic F[Iz,ind] += FCol[I,J,iz+1] / (M[Iz,ind,1] + M[Iz,ind,2])
    if iz == 1 && Iz >  1
      @atomic :monotonic F[Iz-1,ind] += FCol[I,J,iz] / (M[Iz-1,ind,1] + M[Iz-1,ind,2])
    end
    if iz == ColumnTilesDim && Iz <  Nz
      @atomic :monotonic F[Iz+1,ind] += FCol[I,J,iz+2] / (M[Iz+1,ind,1] + M[Iz+1,ind,2])
    end
  end
end

@kernel inbounds = true function DivRhoTrUpwindKernel!(F,@Const(c),@Const(Rho),@Const(uC),@Const(vC),@Const(w),
  @Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob))

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  cCol = @localmem eltype(F) (N,N, ColumnTilesDim+2)
  uCol = @localmem eltype(F) (N,N, ColumnTilesDim+2)
  vCol = @localmem eltype(F) (N,N, ColumnTilesDim+2)
  RhoCol = @localmem eltype(F) (N,N, ColumnTilesDim+2)
  wCol = @localmem eltype(F) (N,N, ColumnTilesDim+1)
  FCol = @localmem eltype(F) (N,N, ColumnTilesDim+2)
  if Iz <= Nz
    wCol[I,J,iz+1] = w[Iz,ind]
    RhoCol[I,J,iz+1] = Rho[Iz,ind]
    cCol[I,J,iz+1] = c[Iz,ind] / RhoCol[I,J,iz+1]
    uCol[I,J,iz+1] = uC[Iz,ind]
    vCol[I,J,iz+1] = vC[Iz,ind]
    FCol[I,J,iz+1] = 0
    if iz == 1 && Iz > 1
      RhoCol[I,J,1] = Rho[Iz-1,ind]
      cCol[I,J,1] = c[Iz-1,ind] / RhoCol[I,J,1]
      uCol[I,J,1] = uC[Iz-1,ind]
      vCol[I,J,1] = vC[Iz-1,ind]
      wCol[I,J,1] = w[Iz,ind]
      FCol[I,J,1] = 0
    elseif iz == 1 && Iz == 1
      RhoCol[I,J,1] = Rho[1,ind]
      cCol[I,J,1] = c[1,ind] / RhoCol[I,J,1]
      wCol[I,J,1] = 0
      FCol[I,J,1] = 0
    end
    if iz == ColumnTilesDim && Iz < Nz
      RhoCol[I,J,ColumnTilesDim+2] = Rho[Iz+1,ind]
      cCol[I,J,ColumnTilesDim+2] = c[Iz+1,ind] / RhoCol[I,J,ColumnTilesDim+2]
      uCol[I,J,ColumnTilesDim+2] = uC[Iz+1,ind]
      vCol[I,J,ColumnTilesDim+2] = vC[Iz+1,ind]
      FCol[I,J,ColumnTilesDim+2] = 0
    elseif iz == ColumnTilesDim && Iz == Nz
      RhoCol[I,J,ColumnTilesDim+2] = Rho[Nz,ind]
      cCol[I,J,ColumnTilesDim+2] = c[Nz,ind] / RhoCol[I,J,ColumnTilesDim+2]
      FCol[I,J,ColumnTilesDim+2] = 0
    end
  end
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz < Nz 
    wCon = RhoCol[I,J,iz+1] * (dXdxI[I,J,2,Iz,3,1,IF] * uCol[I,J,iz+1] + 
      dXdxI[I,J,2,Iz,3,2,IF] * vCol[I,J,iz+1] + dXdxI[I,J,2,Iz,3,3,IF] * wCol[I,J,iz+1]) +
      RhoCol[I,J,iz+2] * (dXdxI[I,J,1,Iz+1,3,1,IF] * uCol[I,J,iz+2] + 
      dXdxI[I,J,2,Iz+1,3,2,IF] * vCol[I,J,iz+2] + dXdxI[I,J,1,Iz+1,3,3,IF] * wCol[I,J,iz+1])
    cL = cCol[I,J,iz+1]
    cR = cCol[I,J,iz+2]
    Flux = 0.25 * ((abs(wCon) + wCon) * cL + (-abs(wCon) + wCon) * cR)
    @atomic :monotonic FCol[I,J,iz+1] += -Flux
    @atomic :monotonic FCol[I,J,iz+2] += Flux
  end 

  if Iz <= Nz
    tempx = -RhoCol[I,J,iz+1] * cCol[I,J,iz+1] * ((dXdxI[I,J,1,Iz,1,1,IF] + dXdxI[I,J,2,Iz,1,1,IF]) * uCol[I,J,iz+1] +
      (dXdxI[I,J,1,Iz,1,2,IF] + dXdxI[I,J,2,Iz,1,2,IF]) * vCol[I,J,iz+1] +
      dXdxI[I,J,1,Iz,1,3,IF] * wCol[I,J,iz] + dXdxI[I,J,2,Iz,1,3,IF] * wCol[I,J,iz+1])
    tempy = -RhoCol[I,J,iz+1] * cCol[I,J,iz+1] * ((dXdxI[I,J,1,Iz,2,1,IF] + dXdxI[I,J,2,Iz,2,1,IF]) * uCol[I,J,iz+1] +
      (dXdxI[I,J,1,Iz,2,2,IF] + dXdxI[I,J,2,Iz,2,2,IF]) * vCol[I,J,iz+1] +
      dXdxI[I,J,1,Iz,2,3,IF] * wCol[I,J,iz] + dXdxI[I,J,2,Iz,2,3,IF] * wCol[I,J,iz+1])
    for k = 1 : N
      @atomic :monotonic FCol[k,J,iz+1] += D[k,I] * tempx
      @atomic :monotonic FCol[I,k,iz+1] += D[k,J] * tempy
    end
  end
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]
  if Iz <= Nz 
    @atomic :monotonic F[Iz,ind] += FCol[I,J,iz+1] / (M[Iz,ind,1] + M[Iz,ind,2])
    if iz == 1 && Iz >  1
      @atomic :monotonic F[Iz-1,ind] += FCol[I,J,iz] / (M[Iz-1,ind,1] + M[Iz-1,ind,2])
    end
    if iz == ColumnTilesDim && Iz <  Nz
      @atomic :monotonic F[Iz+1,ind] += FCol[I,J,iz+2] / (M[Iz+1,ind,1] + M[Iz+1,ind,2])
    end
  end
end

@kernel inbounds = true function DivRhoKEUpwind3Kernel!(F,@Const(U),@Const(p),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob))

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  cCol = @localmem eltype(F) (N,N, ColumnTilesDim+3)
  uConCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  vConCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  if Iz <= Nz
    cCol[I,J,iz+1] = (U[Iz,ind,5] + p[Iz,ind]) / U[Iz,ind,1] 
    @views (uCon, vCon) = Contra12(-U[Iz,ind,1],U[Iz,ind,2],U[Iz,ind,3],dXdxI[1:2,1:2,:,ID,Iz,IF])
    uConCol[I,J,iz] = uCon
    vConCol[I,J,iz] = vCon
  end
  if iz == 1
    Izm1 = max(Iz - 1,1)
    cCol[I,J,iz] = (U[Izm1,ind,5] + p[Izm1,ind]) / U[Izm1,ind,1] 
#   if Iz == 1
#     cCol[I,J,iz] = eltype(F)(2) * cCol[I,J,iz] - U[2,ind,5] / U[2,ind,1]  
#   end  
  end
  if iz == ColumnTilesDim || Iz == Nz
    Izp1 = min(Iz + 1,Nz)
    cCol[I,J,iz+2] = (U[Izp1,ind,5] + p[Izp1,ind]) / U[Izp1,ind,1] 
    Izp2 = min(Iz + 2,Nz)
    cCol[I,J,iz+3] = (U[Izp2,ind,5] + p[Izp2,ind]) / U[Izp2,ind,1] 
  end
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz < Nz 
    cLL = cCol[I,J,iz]
    cL = cCol[I,J,iz+1]
    cR = cCol[I,J,iz+2]
    cRR = cCol[I,J,iz+3]

    @views wCon = Contra3(U[Iz:Iz+1,ind,1],U[Iz:Iz+1,ind,2],U[Iz:Iz+1,ind,3],
      U[Iz,ind,4],dXdxI[3,:,:,ID,Iz:Iz+1,IF])

    Izm1 = max(Iz - 1,1)
    Izp2 = min(Iz + 2, Nz)
    JLL = JJ[ID,1,Izm1,IF] + JJ[ID,2,Izm1,IF]
    JL = JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]
    JR = JJ[ID,1,Iz+1,IF] + JJ[ID,2,Iz+1,IF]
    JRR = JJ[ID,1,Izp2,IF] + JJ[ID,2,Izp2,IF]
    cFL, cFR = RecU4(cLL,cL,cR,cRR,JLL,JL,JR,JRR) 
    Flux = eltype(F)(0.25) * ((abs(wCon) + wCon) * cFL + (-abs(wCon) + wCon) * cFR)
    @atomic :monotonic F[Iz,ind,5] += -Flux / (M[Iz,ind,1] + M[Iz,ind,2])
    @atomic :monotonic F[Iz+1,ind,5] += Flux / (M[Iz+1,ind,1] + M[Iz+1,ind,2])
    Flux = eltype(F)(0.5)*wCon
    @atomic :monotonic F[Iz,ind,1] += -Flux / (M[Iz,ind,1] + M[Iz,ind,2])
    @atomic :monotonic F[Iz+1,ind,1] += Flux / (M[Iz+1,ind,1] + M[Iz+1,ind,2])
  end 

  if Iz <= Nz
    DivRho = D[I,1] * uConCol[1,J,iz] 
    DivRho += D[J,1] * vConCol[I,1,iz] 
    DivRhoTr = D[I,1] * uConCol[1,J,iz] * cCol[1,J,iz+1] 
    DivRhoTr += D[J,1] * vConCol[I,1,iz] * cCol[I,1,iz+1]
    for k = 2 : N
      DivRho += D[I,k] * uConCol[k,J,iz] 
      DivRho += D[J,k] * vConCol[I,k,iz] 
      DivRhoTr += D[I,k] * uConCol[k,J,iz] * cCol[k,J,iz+1] 
      DivRhoTr += D[J,k] * vConCol[I,k,iz] * cCol[I,k,iz+1]
    end
    @atomic :monotonic F[Iz,ind,1] += DivRho / (M[Iz,ind,1] + M[Iz,ind,2])
    @atomic :monotonic F[Iz,ind,5] += DivRhoTr / (M[Iz,ind,1] + M[Iz,ind,2])
  end
end
@kernel inbounds = true function DivRhoThUpwind3Kernel!(F,@Const(U),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob))

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  cCol = @localmem eltype(F) (N,N, ColumnTilesDim+3)
  uConCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  vConCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  if Iz <= Nz
    cCol[I,J,iz+1] = U[Iz,ind,5] / U[Iz,ind,1]
    @views (uCon, vCon) = Contra12(-U[Iz,ind,1],U[Iz,ind,2],U[Iz,ind,3],dXdxI[1:2,1:2,:,ID,Iz,IF])
    uConCol[I,J,iz] = uCon
    vConCol[I,J,iz] = vCon
  end
  if iz == 1
    Izm1 = max(Iz - 1,1)
    cCol[I,J,iz] = U[Izm1,ind,5] / U[Izm1,ind,1]
#   if Iz == 1
#     cCol[I,J,iz] = eltype(F)(2) * cCol[I,J,iz] - U[2,ind,5] / U[2,ind,1]  
#   end  
  end
  if iz == ColumnTilesDim || Iz == Nz
    Izp1 = min(Iz + 1,Nz)
    cCol[I,J,iz+2] = U[Izp1,ind,5] / U[Izp1,ind,1]
    Izp2 = min(Iz + 2,Nz)
    cCol[I,J,iz+3] = U[Izp2,ind,5] / U[Izp2,ind,1]
  end
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz < Nz 
    cLL = cCol[I,J,iz]
    cL = cCol[I,J,iz+1]
    cR = cCol[I,J,iz+2]
    cRR = cCol[I,J,iz+3]

    @views wCon = Contra3(U[Iz:Iz+1,ind,1],U[Iz:Iz+1,ind,2],U[Iz:Iz+1,ind,3],
      U[Iz,ind,4],dXdxI[3,:,:,ID,Iz:Iz+1,IF])

    Izm1 = max(Iz - 1,1)
    Izp2 = min(Iz + 2, Nz)
    JLL = JJ[ID,1,Izm1,IF] + JJ[ID,2,Izm1,IF]
    JL = JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]
    JR = JJ[ID,1,Iz+1,IF] + JJ[ID,2,Iz+1,IF]
    JRR = JJ[ID,1,Izp2,IF] + JJ[ID,2,Izp2,IF]
    cFL, cFR = RecU4(cLL,cL,cR,cRR,JLL,JL,JR,JRR) 
    Flux = eltype(F)(0.25) * ((abs(wCon) + wCon) * cFL + (-abs(wCon) + wCon) * cFR)
    @atomic :monotonic F[Iz,ind,5] += -Flux / (M[Iz,ind,1] + M[Iz,ind,2])
    @atomic :monotonic F[Iz+1,ind,5] += Flux / (M[Iz+1,ind,1] + M[Iz+1,ind,2])
    Flux = eltype(F)(0.5)*wCon
    @atomic :monotonic F[Iz,ind,1] += -Flux / (M[Iz,ind,1] + M[Iz,ind,2])
    @atomic :monotonic F[Iz+1,ind,1] += Flux / (M[Iz+1,ind,1] + M[Iz+1,ind,2])
  end 

  if Iz <= Nz
    DivRho = D[I,1] * uConCol[1,J,iz] 
    DivRho += D[J,1] * vConCol[I,1,iz] 
    DivRhoTr = D[I,1] * uConCol[1,J,iz] * cCol[1,J,iz+1] 
    DivRhoTr += D[J,1] * vConCol[I,1,iz] * cCol[I,1,iz+1]
    for k = 2 : N
      DivRho += D[I,k] * uConCol[k,J,iz] 
      DivRho += D[J,k] * vConCol[I,k,iz] 
      DivRhoTr += D[I,k] * uConCol[k,J,iz] * cCol[k,J,iz+1] 
      DivRhoTr += D[J,k] * vConCol[I,k,iz] * cCol[I,k,iz+1]
    end
    @atomic :monotonic F[Iz,ind,1] += DivRho / (M[Iz,ind,1] + M[Iz,ind,2])
    @atomic :monotonic F[Iz,ind,5] += DivRhoTr / (M[Iz,ind,1] + M[Iz,ind,2])
  end
end

@kernel inbounds = true function DivRhoKernel!(F,@Const(U),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob))

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  uConCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  vConCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  if Iz <= Nz
    @views (uCon, vCon) = Contra12(-U[Iz,ind,1],U[Iz,ind,2],U[Iz,ind,3],dXdxI[1:2,1:2,:,ID,Iz,IF])
    uConCol[I,J,iz] = uCon
    vConCol[I,J,iz] = vCon
  end
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz < Nz 
    @views wCon = Contra3(U[Iz:Iz+1,ind,1],U[Iz:Iz+1,ind,2],U[Iz:Iz+1,ind,3],
      U[Iz,ind,4],dXdxI[3,:,:,ID,Iz:Iz+1,IF])

    Flux = eltype(F)(0.5) * wCon
    @atomic :monotonic F[Iz,ind,1] += -Flux / (M[Iz,ind,1] + M[Iz,ind,2])
    @atomic :monotonic F[Iz+1,ind,1] += Flux / (M[Iz+1,ind,1] + M[Iz+1,ind,2])
  end 

  if Iz <= Nz
    DivRho = D[I,1] * uConCol[1,J,iz] 
    DivRho += D[J,1] * vConCol[I,1,iz] 
    for k = 2 : N
      DivRho += D[I,k] * uConCol[k,J,iz] 
      DivRho += D[J,k] * vConCol[I,k,iz] 
    end
    @atomic :monotonic F[Iz,ind,1] += DivRho / (M[Iz,ind,1] + M[Iz,ind,2])
  end
end

@kernel inbounds = true function DivRhoTrUpwind3Kernel!(FTr,@Const(Tr),@Const(U),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob))

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  cCol = @localmem eltype(FTr) (N,N, ColumnTilesDim+3)
  uConCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  vConCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  if Iz <= Nz
    cCol[I,J,iz+1] = Tr[Iz,ind] / U[Iz,ind,1]
    @views (uCon, vCon) = Contra12(-U[Iz,ind,1],U[Iz,ind,2],U[Iz,ind,3],dXdxI[1:2,1:2,:,ID,Iz,IF])
    uConCol[I,J,iz] = uCon
    vConCol[I,J,iz] = vCon
  end
  if iz == 1
    Izm1 = max(Iz - 1,1)
    cCol[I,J,iz] = Tr[Izm1,ind] / U[Izm1,ind,1]
  end
  if iz == ColumnTilesDim || Iz == Nz
    Izp1 = min(Iz + 1,Nz)
    cCol[I,J,iz+2] = Tr[Izp1,ind] / U[Izp1,ind,1]
    Izp2 = min(Iz + 2,Nz)
    cCol[I,J,iz+3] = Tr[Izp2,ind] / U[Izp2,ind,1]
  end
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz < Nz 
    cLL = cCol[I,J,iz]
    cL = cCol[I,J,iz+1]
    cR = cCol[I,J,iz+2]
    cRR = cCol[I,J,iz+3]

    @views wCon = Contra3(U[Iz:Iz+1,ind,1],U[Iz:Iz+1,ind,2],U[Iz:Iz+1,ind,3],
      U[Iz,ind,4],dXdxI[3,:,:,ID,Iz:Iz+1,IF])

    Izm1 = max(Iz - 1,1)
    Izp2 = min(Iz + 2, Nz)
    JLL = JJ[ID,1,Izm1,IF] + JJ[ID,2,Izm1,IF]
    JL = JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]
    JR = JJ[ID,1,Iz+1,IF] + JJ[ID,2,Iz+1,IF]
    JRR = JJ[ID,1,Izp2,IF] + JJ[ID,2,Izp2,IF]
    cFL, cFR = RecU4(cLL,cL,cR,cRR,JLL,JL,JR,JRR) 
    Flux = eltype(FTr)(0.25) * ((abs(wCon) + wCon) * cFL + (-abs(wCon) + wCon) * cFR)
    @atomic :monotonic FTr[Iz,ind] += -Flux / (M[Iz,ind,1] + M[Iz,ind,2])
    @atomic :monotonic FTr[Iz+1,ind] += Flux / (M[Iz+1,ind,1] + M[Iz+1,ind,2])
  end 

  if Iz <= Nz
    DivRhoTr = D[I,1] * uConCol[1,J,iz] * cCol[1,J,iz+1] 
    DivRhoTr += D[J,1] * vConCol[I,1,iz] * cCol[I,1,iz+1]
    for k = 2 : N
      DivRhoTr += D[I,k] * uConCol[k,J,iz] * cCol[k,J,iz+1] 
      DivRhoTr += D[J,k] * vConCol[I,k,iz] * cCol[I,k,iz+1]
    end
    @atomic :monotonic FTr[Iz,ind] += DivRhoTr / (M[Iz,ind,1] + M[Iz,ind,2])
  end
end

@kernel inbounds = true function DivRhoTrUpwind3New1Kernel!(FTr,@Const(Tr),@Const(U),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob))

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  cCol = @localmem eltype(FTr) (N,N, ColumnTilesDim+3)
  RhoCol = @localmem eltype(FTr) (N,N, ColumnTilesDim+1)
  uCol = @localmem eltype(FTr) (N,N, ColumnTilesDim+1)
  vCol = @localmem eltype(FTr) (N,N, ColumnTilesDim+1)
  wCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  JCol = @localmem eltype(FTr) (N,N, ColumnTilesDim+3)
  if Iz <= Nz
    cCol[I,J,iz+1] = Tr[Iz,ind,1] / U[Iz,ind,1]
    JCol[I,J,iz+1] = JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]
    RhoCol[I,J,iz] = U[Iz,ind,1]
    uCol[I,J,iz] = U[Iz,ind,2]
    vCol[I,J,iz] = U[Iz,ind,3]
    wCol[I,J,iz] = U[Iz,ind,4]
  end
  if iz == 1
    Izm1 = max(Iz - 1,1)
    cCol[I,J,iz] = Tr[Izm1,ind,1] / U[Izm1,ind,1]
    JCol[I,J,iz] = JJ[ID,1,Izm1,IF] + JJ[ID,2,Izm1,IF] 
  end
  if iz == ColumnTilesDim || Iz == Nz
    Izp1 = min(Iz + 1,Nz)
    cCol[I,J,iz+2] = Tr[Izp1,ind,1] / U[Izp1,ind,1]
    JCol[I,J,iz+2] = JJ[ID,1,Izp1,IF] + JJ[ID,2,Izp1,IF] 
    Izp2 = min(Iz + 2,Nz)
    cCol[I,J,iz+3] = Tr[Izp2,ind,1] / U[Izp2,ind,1]
    JCol[I,J,iz+3] = JJ[ID,1,Izp2,IF] + JJ[ID,2,Izp2,IF] 
    RhoCol[I,J,iz+1] = U[Izp1,ind,1]
    uCol[I,J,iz+1] = U[Izp1,ind,2]
    vCol[I,J,iz+1] = U[Izp1,ind,3]
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
    Flux = eltype(FTr)(0.25) * ((abs(wCon) + wCon) * cFL + 
      (-abs(wCon) + wCon) * cFR)
    @atomic :monotonic FTr[Iz,ind,1] += -Flux / (M[Iz,ind,1] + M[Iz,ind,2])
    @atomic :monotonic FTr[Iz+1,ind,1] += Flux / (M[Iz+1,ind,1] + M[Iz+1,ind,2])
  end 

  if Iz <= Nz
    uCon, vCon = Contra12(-RhoCol[I,J,iz],uCol[I,J,iz],vCol[I,J,iz],view(dXdxI,1:2,1:2,:,ID,Iz,IF))
    uCol[I,J,iz] = uCon
    vCol[I,J,iz] = vCon
  end  
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]
  if Iz <= Nz
    DivRhoTr = D[I,1] * uCol[1,J,iz] * cCol[1,J,iz+1] + D[J,1] * vCol[I,1,iz] * cCol[I,1,iz+1]
    for k = 2 : N
      DivRhoTr += D[I,k] * uCol[k,J,iz] * cCol[k,J,iz+1] + D[J,k] * vCol[I,k,iz] * cCol[I,k,iz+1]
    end
    @atomic :monotonic FTr[Iz,ind,1] += DivRhoTr / (M[Iz,ind,1] + M[Iz,ind,2])
  end
  
# Second tracer  
  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]
  if Iz <= Nz
    cCol[I,J,iz+1] = Tr[Iz,ind,2] / U[Iz,ind,1]
  end
  if iz == 1
    Izm1 = max(Iz - 1,1)
    cCol[I,J,iz] = Tr[Izm1,ind,2] / U[Izm1,ind,1]
  end
  if iz == ColumnTilesDim || Iz == Nz
    Izp1 = min(Iz + 1,Nz)
    cCol[I,J,iz+2] = Tr[Izp1,ind,2] / U[Izp1,ind,1]
    Izp2 = min(Iz + 2,Nz)
    cCol[I,J,iz+3] = Tr[Izp2,ind,2] / U[Izp2,ind,1]
  end
  @synchronize

  ID = I + (J - 1) * N
  ind = Glob[ID,IF]
  if Iz < Nz
    wCon = wCol[I,J,iz]
    cFL, cFR = RecU4(cCol[I,J,iz],cCol[I,J,iz+1],cCol[I,J,iz+2],cCol[I,J,iz+3],
      JCol[I,J,iz],JCol[I,J,iz+1],JCol[I,J,iz+2],JCol[I,J,iz+3]) 
    Flux = eltype(FTr)(0.25) * ((abs(wCon) + wCon) * cFL + 
      (-abs(wCon) + wCon) * cFR)
    @atomic :monotonic FTr[Iz,ind,2] += -Flux / (M[Iz,ind,1] + M[Iz,ind,2])
    @atomic :monotonic FTr[Iz+1,ind,2] += Flux / (M[Iz+1,ind,1] + M[Iz+1,ind,2])
  end
  if Iz <= Nz
    DivRhoTr = D[I,1] * uCol[1,J,iz] * cCol[1,J,iz+1] + D[J,1] * vCol[I,1,iz] * cCol[I,1,iz+1]
    for k = 2 : N
      DivRhoTr += D[I,k] * uCol[k,J,iz] * cCol[k,J,iz+1] + D[J,k] * vCol[I,k,iz] * cCol[I,k,iz+1]
    end
    @atomic :monotonic FTr[Iz,ind,2] += DivRhoTr / (M[Iz,ind,1] + M[Iz,ind,2])
  end
end

@kernel inbounds = true function DivRhoTrUpwind3NewKernel!(FTr,@Const(Tr),@Const(U),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob))

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  cCol = @localmem eltype(FTr) (N,N, ColumnTilesDim+3)
  RhoCol = @localmem eltype(FTr) (N,N, ColumnTilesDim+1)
  uCol = @localmem eltype(FTr) (N,N, ColumnTilesDim+1)
  vCol = @localmem eltype(FTr) (N,N, ColumnTilesDim+1)
  wCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  JCol = @localmem eltype(FTr) (N,N, ColumnTilesDim+3)
  if Iz <= Nz
    cCol[I,J,iz+1] = Tr[Iz,ind] / U[Iz,ind,1]
    JCol[I,J,iz+1] = JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]
    RhoCol[I,J,iz] = U[Iz,ind,1]
    uCol[I,J,iz] = U[Iz,ind,2]
    vCol[I,J,iz] = U[Iz,ind,3]
    wCol[I,J,iz] = U[Iz,ind,4]
  end
  if iz == 1
    Izm1 = max(Iz - 1,1)
    cCol[I,J,iz] = Tr[Izm1,ind] / U[Izm1,ind,1]
    JCol[I,J,iz] = JJ[ID,1,Izm1,IF] + JJ[ID,2,Izm1,IF] 
  end
  if iz == ColumnTilesDim || Iz == Nz
    Izp1 = min(Iz + 1,Nz)
    cCol[I,J,iz+2] = Tr[Izp1,ind] / U[Izp1,ind,1]
    JCol[I,J,iz+2] = JJ[ID,1,Izp1,IF] + JJ[ID,2,Izp1,IF] 
    Izp2 = min(Iz + 2,Nz)
    cCol[I,J,iz+3] = Tr[Izp2,ind] / U[Izp2,ind,1]
    JCol[I,J,iz+3] = JJ[ID,1,Izp2,IF] + JJ[ID,2,Izp2,IF] 
    RhoCol[I,J,iz+1] = U[Izp1,ind,1]
    uCol[I,J,iz+1] = U[Izp1,ind,2]
    vCol[I,J,iz+1] = U[Izp1,ind,3]
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
    Flux = eltype(FTr)(0.25) * ((abs(wCon) + wCon) * cFL + 
      (-abs(wCon) + wCon) * cFR)
    @atomic :monotonic FTr[Iz,ind] += -Flux / (M[Iz,ind,1] + M[Iz,ind,2])
    @atomic :monotonic FTr[Iz+1,ind] += Flux / (M[Iz+1,ind,1] + M[Iz+1,ind,2])
  end 
  if Iz <= Nz
    uCon, vCon = Contra12(-RhoCol[I,J,iz],uCol[I,J,iz],vCol[I,J,iz],view(dXdxI,1:2,1:2,:,ID,Iz,IF))
    uCol[I,J,iz] = uCon
    vCol[I,J,iz] = vCon
  end  
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]
  if Iz <= Nz
    DivRhoTr = D[I,1] * uCol[1,J,iz] * cCol[1,J,iz+1] + D[J,1] * vCol[I,1,iz] * cCol[I,1,iz+1]
    for k = 2 : N
      DivRhoTr += D[I,k] * uCol[k,J,iz] * cCol[k,J,iz+1] + D[J,k] * vCol[I,k,iz] * cCol[I,k,iz+1]
    end
    @atomic :monotonic FTr[Iz,ind] += DivRhoTr / (M[Iz,ind,1] + M[Iz,ind,2])
  end
end

@kernel inbounds = true function AdvectionTrUpwind3Kernel!(FTr,@Const(Tr),@Const(U),@Const(w),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob))

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF,IE = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]
  NE = @uniform @ndrange()[5]

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  cCol = @localmem eltype(FTr) (N,N, ColumnTilesDim+3)
  if Iz <= Nz
    cCol[I,J,iz+1] = Tr[Iz,ind,IE]
  end
  if iz == 1
    Izm1 = max(Iz - 1,1)
    cCol[I,J,iz] = Tr[Izm1,ind,IE] 
  end
  if iz == ColumnTilesDim || Iz == Nz
    Izp1 = min(Iz + 1,Nz)
    cCol[I,J,iz+2] = Tr[Izp1,ind,IE]
    Izp2 = min(Iz + 2,Nz)
    cCol[I,J,iz+3] = Tr[Izp2,ind,IE]
  end
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz < Nz 
    cLL = cCol[I,J,iz]
    cL = cCol[I,J,iz+1]
    cR = cCol[I,J,iz+2]
    cRR = cCol[I,J,iz+3]

    @views wCon = Contra3(U[Iz:Iz+1,ind,2],U[Iz:Iz+1,ind,3],
      w[Iz,ind],dXdxI[3,:,:,ID,Iz:Iz+1,IF])

    Izm1 = max(Iz - 1,1)
    Izp2 = min(Iz + 2, Nz)
    JLL = JJ[ID,1,Izm1,IF] + JJ[ID,2,Izm1,IF]
    JL = JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]
    JR = JJ[ID,1,Iz+1,IF] + JJ[ID,2,Iz+1,IF]
    JRR = JJ[ID,1,Izp2,IF] + JJ[ID,2,Izp2,IF]
    cFL, cFR = RecU4(cLL,cL,cR,cRR,JLL,JL,JR,JRR) 
    Flux = eltype(FTr)(0.25) * ((abs(wCon) + wCon) * cFL + (-abs(wCon) + wCon) * cFR)
    @atomic :monotonic FTr[Iz,ind,IE] += (-Flux + wCon * cL) / (M[Iz,ind,1] + M[Iz,ind,2])
    @atomic :monotonic FTr[Iz+1,ind,IE] += (Flux - wCon * cR) / (M[Iz+1,ind,1] + M[Iz+1,ind,2])
  end 

  if Iz <= Nz
    GradxTr = D[I,1] * cCol[1,J,iz+1] 
    GradyTr = D[J,1] * cCol[I,1,iz+1]
    for k = 2 : N
      GradxTr += D[I,k] * cCol[k,J,iz+1] 
      GradyTr += D[J,k] * cCol[I,k,iz+1]
    end
    @views (uCon, vCon) = Contra12(U[Iz,ind,2],U[Iz,ind,3],dXdxI[1:2,1:2,:,ID,Iz,IF])
    @atomic :monotonic FTr[Iz,ind,IE] += -(GradxTr * uCon + GradyTr * vCon) / (M[Iz,ind,1] + M[Iz,ind,2])
  end
end

@kernel inbounds = true function DivRhoEDMFKernel!(F,@Const(U),@Const(w),@Const(Rho),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob))

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF,IE = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]
  NE = @uniform @ndrange()[5]

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  uConCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  vConCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  if Iz <= Nz
    @views (uCon, vCon) = Contra12(-Rho[Iz,ind,IE],U[Iz,ind,2],U[Iz,ind,3],dXdxI[1:2,1:2,:,ID,Iz,IF])
    uConCol[I,J,iz] = uCon
    vConCol[I,J,iz] = vCon
  end
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz < Nz 
    @views wCon = Contra3(Rho[Iz:Iz+1,ind,IE],U[Iz:Iz+1,ind,2],U[Iz:Iz+1,ind,3],
      w[Iz,ind,IE],dXdxI[3,:,:,ID,Iz:Iz+1,IF])

    Flux = eltype(F)(0.5) * wCon
    @atomic :monotonic F[Iz,ind,IE] += -Flux / (M[Iz,ind,1] + M[Iz,ind,2])
    @atomic :monotonic F[Iz+1,ind,IE] += Flux / (M[Iz+1,ind,1] + M[Iz+1,ind,2])
  end 

  if Iz <= Nz
    DivRho = D[I,1] * uConCol[1,J,iz] 
    DivRho += D[J,1] * vConCol[I,1,iz] 
    for k = 2 : N
      DivRho += D[I,k] * uConCol[k,J,iz] 
      DivRho += D[J,k] * vConCol[I,k,iz] 
    end
    @atomic :monotonic F[Iz,ind,1] += DivRho / (M[Iz,ind,1] + M[Iz,ind,2])
  end
end

@inline function Contra12(Rho,u,v,dXdxI)
  uCon = Rho * ((dXdxI[1,1,1] + dXdxI[1,1,2]) * u +
  (dXdxI[1,2,1] + dXdxI[1,2,2]) * v)
  vCon = Rho * ((dXdxI[2,1,1] + dXdxI[2,1,2]) * u +
  (dXdxI[2,2,1] + dXdxI[2,2,2]) * v)
  return uCon, vCon
end

@inline function Contra12(u,v,dXdxI)
  uCon = (dXdxI[1,1,1] + dXdxI[1,1,2]) * u +
  (dXdxI[1,2,1] + dXdxI[1,2,2]) * v
  vCon = (dXdxI[2,1,1] + dXdxI[2,1,2]) * u +
  (dXdxI[2,2,1] + dXdxI[2,2,2]) * v
  return uCon, vCon
end

@inline function Grad12(Rho,u,v,dXdxI,J)
  uCon = Rho * ((dXdxI[1,1,1] + dXdxI[1,1,2]) * u +
  (dXdxI[2,1,1] + dXdxI[2,1,2]) * v) / (J[1] + J[2])
  vCon = Rho * ((dXdxI[1,2,1] + dXdxI[1,2,2]) * u +
  (dXdxI[2,2,1] + dXdxI[2,2,2]) * v) / (J[1] + J[2])
  return uCon, vCon
end

@inline function Grad12(u,v,dXdxI,J)
  uCon = ((dXdxI[1,1,1] + dXdxI[1,1,2]) * u +
  (dXdxI[2,1,1] + dXdxI[2,1,2]) * v) / (J[1] + J[2])
  vCon = ((dXdxI[1,2,1] + dXdxI[1,2,2]) * u +
  (dXdxI[2,2,1] + dXdxI[2,2,2]) * v) / (J[1] + J[2])
  return uCon, vCon
end

@inline function Grad12(u,v,dXdxI)
  uCon = (dXdxI[1,1,1] + dXdxI[1,1,2]) * u +
  (dXdxI[2,1,1] + dXdxI[2,1,2]) * v
  vCon = (dXdxI[1,2,1] + dXdxI[1,2,2]) * u +
  (dXdxI[2,2,1] + dXdxI[2,2,2]) * v
  return uCon, vCon
end

@inline function Grad3(u,v,dXdxI)
  wCon1 = dXdxI[1,3,1] * u + dXdxI[2,3,1] * v
  wCon2 = dXdxI[1,3,2] * u + dXdxI[2,3,2] * v
  return wCon1, wCon2
end

@inline function Curl12(u,v,dXdxI)
  uCon = (dXdxI[1,1,1] + dXdxI[1,1,2]) * v -
  (dXdxI[1,2,1] + dXdxI[1,2,2]) * u
  vCon = (dXdxI[2,1,1] + dXdxI[2,1,2]) * v -
  (dXdxI[2,2,1] + dXdxI[2,2,2]) * u
  return uCon, vCon
end

@inline function Rot12(u,v,dXdxI)
  uCon = (dXdxI[1,2,1] + dXdxI[1,2,2]) * u +
  (dXdxI[2,2,1] + dXdxI[2,2,2]) * v
  vCon = -(dXdxI[1,1,1] + dXdxI[1,1,2]) * u -
  (dXdxI[2,1,1] + dXdxI[2,1,2]) * v
  return uCon, vCon
end

@inline function Contra3(u,v,w,dXdxI)
  wCon = dXdxI[1,2,1] * u[1] + dXdxI[2,2,1] * v[1] + dXdxI[3,2,1] * w + 
    dXdxI[1,1,2] * u[2] + dXdxI[2,1,2] * v[2] + dXdxI[3,1,2] * w
end

@inline function Contra3(Rho,u,v,w,dXdxI)
  wCon = Rho[1] * (dXdxI[1,2,1] * u[1] + dXdxI[2,2,1] * v[1] + dXdxI[3,2,1] * w) + 
    Rho[2] * (dXdxI[1,1,2] * u[2] + dXdxI[2,1,2] * v[2] + dXdxI[3,1,2] * w)
end

@inline function Contra3W(Rho,u,v,w,dXdxI,J)
  wCon = eltype(Rho)(2) * (J[2,1] * Rho[1] * (dXdxI[1,2,1] * u[1] + dXdxI[2,2,1] * v[1] + dXdxI[3,2,1] * w) + 
    J[1,2] * Rho[2] * (dXdxI[1,1,2] * u[2] + dXdxI[2,1,2] * v[2] + dXdxI[3,1,2] * w)) /
    (J[2,1] + J[1,2])
end
  
@inline function RecU41(cLL,cL,cR,cRR,JLL,JL,JR,JRR)

  kR = (JL / (JL + JR)) * ((JLL + JL) / (JLL + JL + JR))
  kL = -(JL / (JLL + JL)) * (JR / (JLL + JL + JR))
  cFL = kL * cLL + (1 - kL - kR)*cL + kR * cR

  kL = (JR / (JR + JL)) * ((JRR + JR)/(JL + JR + JRR))
  kR = -(JR /(JRR + JR)) *(JL /(JL + JR + JRR))
  cFR = kL * cL + (1 - kL - kR) * cR + kR * cRR
 
  return (cFL,cFR)
end

@inline function RecU4(cLL,cL,cR,cRR,JLL,JL,JR,JRR)

  kR = (JL / (JL + JR)) * ((JLL + JL) / (JLL + JL + JR))
  kL = (JL / (JLL + JL)) * (JR / (JLL + JL + JR))
  #r = (cL - cLL) / (cR - cL + sign(cR- cL)*1.e-20 + 1.e-20)
  r = (cL - cLL) / (cR - cL + 1.e-20)
  cFL = cL + max(eltype(cLL)(0), min(r, min(kL * r + kR, eltype(cLL)(1)))) * (cR - cL) 

  kL = (JR / (JR + JL)) * ((JRR + JR)/(JL + JR + JRR))
  kR = (JR /(JRR + JR)) *(JL /(JL + JR + JRR))
# r = (cR - cRR) / (cL - cR  + sign(cL- cR)*1.e-20 + 1.e-20)
  r = (cR - cRR) / (cL - cR  + 1.e-20)
  cFR = cR + max(eltype(cLL)(0), min(r, min(kR * r + kL, eltype(cLL)(1)))) * (cL - cR) 
 
  return (cFL,cFR)
end

@kernel inbounds = true function DampKernel!(Damp,F,U,zP)
  Iz,IC = @index(Global, NTuple)
  NumG = @uniform @ndrange()[2]

  if IC <= NumG
    Fu,Fv,Fw = Damp(zP[Iz,IC],view(U,Iz,IC,1:5))
    F[Iz,IC,2] += Fu
    F[Iz,IC,3] += Fv
    F[Iz,IC,4] += Fw
  end
end  

@kernel inbounds = true function ForceKernel!(Force,F,U,p,xS)
  Iz,IC = @index(Global, NTuple)
  NumG = @uniform @ndrange()[2]

  if IC <= NumG
    Force(view(F,Iz,IC,:),view(U,Iz,IC,:),p[Iz,IC],xS[2,IC])
#   F[Iz,IC,1] += FRho
#   F[Iz,IC,2] += Fu
#   F[Iz,IC,3] += Fv
#   F[Iz,IC,4] += Fw
#   F[Iz,IC,5] += FRhoTh
  end
end  

@kernel inbounds = true function MicrophysicsKernel!(Source,F,U,p)
  Iz,IC = @index(Global, NTuple)
  NumG = @uniform @ndrange()[2]

  if IC <= NumG
    Source(view(F,Iz,IC,:),view(U,Iz,IC,:),p[Iz,IC])
  end
end  

@kernel inbounds = true function VerticalDiffusionScalarKernel!(FTr,@Const(Tr),@Const(Rho),@Const(K),
  @Const(dz))
  iz, = @index(Local, NTuple)
  Iz,IC = @index(Global, NTuple)

  Nz = @uniform @ndrange()[1]
  NumG = @uniform @ndrange()[2]

  if Iz < Nz && IC <= NumG
    dzT = dz[Iz+1,IC]   
    dzB = dz[Iz,IC]   
    grad = eltype(FTr)(2) * K[Iz,IC] * (Tr[Iz+1,IC] / Rho[Iz+1,IC] - 
      Tr[Iz,IC] / Rho[Iz,IC]) / (dzT + dzB)
    @atomic :monotonic FTr[Iz,IC] +=  grad / dzB
    @atomic :monotonic FTr[Iz+1,IC] +=  -grad / dzT
  end  
end  

@kernel inbounds = true function VerticalDiffusionScalarNewKernel!(FTr,@Const(Tr),@Const(Rho),@Const(K),
  @Const(dz))
  iz,iC = @index(Local, NTuple)
  Iz,IC = @index(Global, NTuple)

  nz = @uniform @groupsize()[1]
  NodeTiles = @uniform @groupsize()[2]
  Nz = @uniform @ndrange()[1]
  NumG = @uniform @ndrange()[2]

  qLoc = @localmem eltype(FTr) (nz,NodeTiles)
  dzLoc = @localmem eltype(FTr) (nz,NodeTiles)
  KLoc = @localmem eltype(FTr) (nz,NodeTiles)
  if Iz <= Nz && IC <= NumG
    qLoc[iz,iC] = Tr[Iz,IC] / Rho[Iz,IC]
    dzLoc[iz,iC] = dz[Iz,iC]
  end  
  if Iz < Nz && IC <= NumG
    KLoc[iz,iC] = K[Iz,IC]
  end  

  @synchronize

  if IC <= NumG
    if Iz < Nz  
      gradT = eltype(FTr)(2) * KLoc[iz,iC] * (qLoc[iz+1,iC] - 
        qLoc[iz,iC]) / (dzLoc[iz+1,iC] + dzLoc[iz,iC])
    else
      gradT = eltype(FTr)(0)  
    end    
    if Iz > 1  
      gradB = eltype(FTr)(2) * KLoc[iz-1,iC] * (qLoc[iz,iC] - 
        qLoc[iz-1,iC]) / (dzLoc[iz,iC] + dzLoc[iz-1,iC])
    else
      gradB = eltype(FTr)(0)  
    end    
    FTr[Iz,IC]  += (gradT - gradB) / dzLoc[iz,iC]
  end  
end  

@kernel inbounds = true function SurfaceFluxScalarsKernel!(SurfaceFluxRhs!,F,@Const(U),@Const(p),
  @Const(TSurf),@Const(RhoVSurf),@Const(uStar),@Const(CT),@Const(CH),@Const(dz))

  IC, = @index(Global, NTuple)

  NumG = @uniform @ndrange()[1]

  if IC <= NumG
    SurfaceFluxRhs!(view(F,1,IC,:),view(U,1,IC,:),p[1,IC],dz[1,IC],
      uStar[IC],CT[IC],CH[IC],TSurf[IC],RhoVSurf[IC])  
  end  
end

