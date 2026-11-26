@kernel inbounds = true function HyperViscKernel!(F,@Const(U),@Const(D),@Const(DW),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob)) 

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N
  ind = Glob[ID,IF]

  ThCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  uCCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  vCCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  uDCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  vDCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  Curl = @localmem eltype(F) (N,N, ColumnTilesDim)
  Div = @localmem eltype(F) (N,N, ColumnTilesDim)
  ThCxCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  ThCyCol = @localmem eltype(F) (N,N, ColumnTilesDim)

  if Iz <= Nz
    @views uC, vC = Curl12(U[1,Iz,ind,2],U[1,Iz,ind,3],dXdxI[1:2,1:2,:,ID,Iz,IF])
    uCCol[I,J,iz] = uC
    vCCol[I,J,iz] = vC
    @views uD, vD = Contra12(U[1,Iz,ind,2],U[1,Iz,ind,3],dXdxI[1:2,1:2,:,ID,Iz,IF])
    uDCol[I,J,iz] = uD
    vDCol[I,J,iz] = vD
    ThCol[I,J,iz] = U[1,Iz,ind,5] / U[1,Iz,ind,1]
  end
  @synchronize

  ID = I + (J - 1) * N
  ind = Glob[ID,IF]

  if Iz <= Nz
    Dxc = D[I,1] * ThCol[1,J,iz]
    Dyc = D[J,1] * ThCol[I,1,iz]
    Curl[I,J,iz] = D[I,1] * uCCol[1,J,iz] + D[J,1] * vCCol[I,1,iz] 
    Div[I,J,iz] = D[I,1] * uDCol[1,J,iz] + D[J,1] * vDCol[I,1,iz] 
    for k = 2 : N
      Dxc += D[I,k] * ThCol[k,J,iz]
      Dyc += D[J,k] * ThCol[I,k,iz] 
      Curl[I,J,iz] += D[I,k] * uCCol[k,J,iz] + D[J,k] * vCCol[I,k,iz] 
      Div[I,J,iz] += D[I,k] * uDCol[k,J,iz] + D[J,k] * vDCol[I,k,iz] 
    end
    Curl[I,J,iz] /= (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF])
    Div[I,J,iz] /= (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF])
    @views (GradDx, GradDy) = Grad12(Dxc,Dyc,dXdxI[1:2,1:2,:,ID,Iz,IF],JJ[ID,:,Iz,IF])
    @views (tempx, tempy) = Contra12(GradDx,GradDy,dXdxI[1:2,1:2,:,ID,Iz,IF])
    ThCxCol[I,J,iz] = tempx
    ThCyCol[I,J,iz] = tempy
  end

  @synchronize 

  ID = I + (J - 1) * N
  ind = Glob[ID,IF]

  if Iz <= Nz
    DxCurl = DW[I,1] * Curl[1,J,iz]
    DyCurl = DW[J,1] * Curl[I,1,iz]
    DxDiv = DW[I,1] * Div[1,J,iz]
    DyDiv = DW[J,1] * Div[I,1,iz]
    DivTh = DW[I,1] * ThCxCol[1,J,iz] + DW[J,1] * ThCyCol[I,1,iz]
    for k = 2 : N
      DxCurl += DW[I,k] * Curl[k,J,iz]
      DyCurl += DW[J,k] * Curl[I,k,iz]
      DxDiv += DW[I,k] * Div[k,J,iz]
      DyDiv += DW[J,k] * Div[I,k,iz]
      DivTh += DW[I,k] * ThCxCol[k,J,iz] + DW[J,k] * ThCyCol[I,k,iz]
    end
    @views FuC, FvC = Rot12(DxCurl,DyCurl,dXdxI[1:2,1:2,:,ID,Iz,IF])
    @views FuD, FvD = Grad12(DxDiv,DyDiv,dXdxI[1:2,1:2,:,ID,Iz,IF]) 
    @atomic :monotonic F[1,Iz,ind,1] += -FuC / (M[Iz,ind,1] + M[Iz,ind,2])
    @atomic :monotonic F[1,Iz,ind,2] += -FvC / (M[Iz,ind,1] + M[Iz,ind,2])
    @atomic :monotonic F[1,Iz,ind,3] += FuD / (M[Iz,ind,1] + M[Iz,ind,2])
    @atomic :monotonic F[1,Iz,ind,4] += FvD / (M[Iz,ind,1] + M[Iz,ind,2])
    @atomic :monotonic F[1,Iz,ind,5] += DivTh / (M[Iz,ind,1] + M[Iz,ind,2])
  end
end

@kernel inbounds = true function HyperViscKoeffKernel!(F,@Const(U),@Const(Cache),@Const(D),@Const(DW),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob),KoeffCurl,KoeffGrad,KoeffDiv)

# gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ThCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  uCCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  vCCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  uDCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  vDCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  Curl = @localmem eltype(F) (N,N, ColumnTilesDim)
  Div = @localmem eltype(F) (N,N, ColumnTilesDim)
  ThCxCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  ThCyCol = @localmem eltype(F) (N,N, ColumnTilesDim)

  ID = I + (J - 1) * N
  ind = Glob[ID,IF]

  if Iz <= Nz
    @views uC, vC = Curl12(Cache[1,Iz,ind,1],Cache[1,Iz,ind,2],dXdxI[1:2,1:2,:,ID,Iz,IF])
    uCCol[I,J,iz] = uC
    vCCol[I,J,iz] = vC
    @views uD, vD = Contra12(Cache[1,Iz,ind,3],Cache[1,Iz,ind,4],dXdxI[1:2,1:2,:,ID,Iz,IF])
    uDCol[I,J,iz] = uD
    vDCol[I,J,iz] = vD
    ThCol[I,J,iz] = Cache[1,Iz,ind,5] 
  end
  @synchronize

  ID = I + (J - 1) * N
  ind = Glob[ID,IF]

  if Iz <= Nz
    Dxc = D[I,1] * ThCol[1,J,iz]
    Dyc = D[J,1] * ThCol[I,1,iz]
    Curl[I,J,iz] = D[I,1] * uCCol[1,J,iz] + D[J,1] * vCCol[I,1,iz]
    Div[I,J,iz] = D[I,1] * uDCol[1,J,iz] + D[J,1] * vDCol[I,1,iz]
    for k = 2 : N
      Dxc += D[I,k] * ThCol[k,J,iz]
      Dyc += D[J,k] * ThCol[I,k,iz]
      Curl[I,J,iz] += D[I,k] * uCCol[k,J,iz] + D[J,k] * vCCol[I,k,iz]
      Div[I,J,iz] += D[I,k] * uDCol[k,J,iz] + D[J,k] * vDCol[I,k,iz]
    end
    Curl[I,J,iz] /= (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF])
    Div[I,J,iz] /= (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF])
    @views (GradDx, GradDy) = Grad12(Dxc,Dyc,dXdxI[1:2,1:2,:,ID,Iz,IF],JJ[ID,:,Iz,IF])
    @views (tempx, tempy) = Contra12(U[1,Iz,ind,1],GradDx,GradDy,dXdxI[1:2,1:2,:,ID,Iz,IF])
    ThCxCol[I,J,iz] = tempx
    ThCyCol[I,J,iz] = tempy
  end

  @synchronize

  ID = I + (J - 1) * N
  ind = Glob[ID,IF]
  if Iz <= Nz
    DxCurl = DW[I,1] * Curl[1,J,iz]
    DyCurl = DW[J,1] * Curl[I,1,iz]
    DxDiv = DW[I,1] * Div[1,J,iz]
    DyDiv = DW[J,1] * Div[I,1,iz]
    DivTh = DW[I,1] * ThCxCol[1,J,iz] + DW[J,1] * ThCyCol[I,1,iz]
    for k = 2 : N
      DxCurl += DW[I,k] * Curl[k,J,iz]
      DyCurl += DW[J,k] * Curl[I,k,iz]
      DxDiv += DW[I,k] * Div[k,J,iz]
      DyDiv += DW[J,k] * Div[I,k,iz]
      DivTh += DW[I,k] * ThCxCol[k,J,iz] + DW[J,k] * ThCyCol[I,k,iz]
    end
    @views FuC, FvC = Rot12(DxCurl,DyCurl,dXdxI[1:2,1:2,:,ID,Iz,IF])
    @views FuD, FvD = Grad12(DxDiv,DyDiv,dXdxI[1:2,1:2,:,ID,Iz,IF])
    @atomic :monotonic F[1,Iz,ind,2] += -(-KoeffCurl * FuC + KoeffGrad * FuD) / (M[Iz,ind,1] + M[Iz,ind,2])
    @atomic :monotonic F[1,Iz,ind,3] += -(-KoeffCurl * FvC + KoeffGrad * FvD) / (M[Iz,ind,1] + M[Iz,ind,2])
    @atomic :monotonic F[1,Iz,ind,5] += -KoeffDiv * DivTh / (M[Iz,ind,1] + M[Iz,ind,2])
  end
end

# The metric (a_ij = dXdxI[i,j,1] + dXdxI[i,j,2], J), the density and the nodal mass are loaded
# once into registers, the loop over the tracers only loads the tracer.
@kernel inbounds = true function HyperViscTracerKernel!(FTr,NumTr,@Const(Tr),@Const(Rho),@Const(D),@Const(DW),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob)) 

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  TrCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  TrCxCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  TrCyCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)

  FT = eltype(FTr)
  ID = I + (J - 1) * N  
  Inside = Iz <= Nz && IF <= NF
  ind = Inside ? Glob[ID,IF] : 1
  a11 = a12 = a21 = a22 = JInv = RhoInv = MInv = FT(0)
  if Inside
    a11 = dXdxI[1,1,1,ID,Iz,IF] + dXdxI[1,1,2,ID,Iz,IF]
    a12 = dXdxI[1,2,1,ID,Iz,IF] + dXdxI[1,2,2,ID,Iz,IF]
    a21 = dXdxI[2,1,1,ID,Iz,IF] + dXdxI[2,1,2,ID,Iz,IF]
    a22 = dXdxI[2,2,1,ID,Iz,IF] + dXdxI[2,2,2,ID,Iz,IF]
    JInv = FT(1) / (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF])
    RhoInv = FT(1) / Rho[1,Iz,ind]
    MInv = FT(1) / (M[Iz,ind,1] + M[Iz,ind,2])
  end

  for it = 1 : NumTr
    if Inside
      TrCol[I,J,iz] = Tr[1,Iz,ind,it] * RhoInv
    end
    @synchronize

    if Inside
      Dxc = D[I,1] * TrCol[1,J,iz]
      Dyc = D[J,1] * TrCol[I,1,iz]
      for k = 2 : N
        Dxc += D[I,k] * TrCol[k,J,iz]
        Dyc += D[J,k] * TrCol[I,k,iz] 
      end
      GradDx = (a11 * Dxc + a21 * Dyc) * JInv
      GradDy = (a12 * Dxc + a22 * Dyc) * JInv
      TrCxCol[I,J,iz] = a11 * GradDx + a12 * GradDy
      TrCyCol[I,J,iz] = a21 * GradDx + a22 * GradDy
    end

    @synchronize 

    if Inside
      DivTr = DW[I,1] * TrCxCol[1,J,iz] + DW[J,1] * TrCyCol[I,1,iz]
      for k = 2 : N
        DivTr += DW[I,k] * TrCxCol[k,J,iz] + DW[J,k] * TrCyCol[I,k,iz]
      end
      @atomic :monotonic FTr[1,Iz,ind,it] += DivTr * MInv
    end
  end
end

@kernel inbounds = true function HyperViscWKernel!(Fw,@Const(w),@Const(D),@Const(DW),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob)) 

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  wCol = @localmem eltype(Fw) (N,N, ColumnTilesDim)
  wCxCol = @localmem eltype(Fw) (N,N, ColumnTilesDim)
  wCyCol = @localmem eltype(Fw) (N,N, ColumnTilesDim)
  if Iz < Nz && IF <= NF
    wCol[I,J,iz] = w[Iz,ind] 
  end
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz < Nz && IF <= NF
    DxcF = D[I,1] * wCol[1,J,iz]
    DycF = D[J,1] * wCol[I,1,iz]
    for k = 2 : N
      DxcF += D[I,k] * wCol[k,J,iz]
      DycF += D[J,k] * wCol[I,k,iz] 
    end
    GradDx = ((dXdxI[1,1,2,ID,Iz,IF] + dXdxI[1,1,1,ID,Iz+1,IF]) * DxcF +
      (dXdxI[2,1,2,ID,Iz,IF] + dXdxI[2,1,1,ID,Iz+1,IF]) * DycF) / (JJ[ID,2,Iz,IF] + JJ[ID,1,Iz+1,IF])
    GradDy = ((dXdxI[1,2,2,ID,Iz,IF] + dXdxI[1,2,1,ID,Iz+1,IF]) * DxcF +
      (dXdxI[2,2,2,ID,Iz,IF] + dXdxI[2,2,1,ID,Iz+1,IF]) * DycF) / (JJ[ID,2,Iz,IF] + JJ[ID,1,Iz+1,IF])
    tempx = (dXdxI[1,1,2,ID,Iz,IF] + dXdxI[1,1,1,ID,Iz+1,IF]) * GradDx +
      (dXdxI[1,2,2,ID,Iz,IF] + dXdxI[1,2,1,ID,Iz+1,IF]) * GradDy
    tempy = (dXdxI[2,1,2,ID,Iz,IF] + dXdxI[2,1,1,ID,Iz+1,IF]) * GradDx +
      (dXdxI[2,2,2,ID,Iz,IF] + dXdxI[2,2,1,ID,Iz+1,IF]) * GradDy
    wCxCol[I,J,iz] = tempx
    wCyCol[I,J,iz] = tempy
  end

  @synchronize 

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]
  if Iz < Nz && IF <= NF
    Divw = DW[I,1] * wCxCol[1,J,iz] + DW[J,1] * wCyCol[I,1,iz]
    for k = 2 : N
      Divw += DW[I,k] * wCxCol[k,J,iz] + DW[J,k] * wCyCol[I,k,iz]
    end
    @atomic :monotonic Fw[Iz,ind] += Divw / (M[Iz,ind,2] + M[Iz+1,ind,1])
  end
end

# as HyperViscTracerKernel!, metric, density and nodal mass loaded once
@kernel inbounds = true function HyperViscTracerKoeffKernel!(FTr,NumTr,@Const(Cache),@Const(Rho),@Const(D),@Const(DW),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob),KoeffDiv) 

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  TrCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  TrCxCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  TrCyCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)

  FT = eltype(FTr)
  ID = I + (J - 1) * N  
  Inside = Iz <= Nz && IF <= NF
  ind = Inside ? Glob[ID,IF] : 1
  a11 = a12 = a21 = a22 = JInv = RhoL = MInv = FT(0)
  if Inside
    a11 = dXdxI[1,1,1,ID,Iz,IF] + dXdxI[1,1,2,ID,Iz,IF]
    a12 = dXdxI[1,2,1,ID,Iz,IF] + dXdxI[1,2,2,ID,Iz,IF]
    a21 = dXdxI[2,1,1,ID,Iz,IF] + dXdxI[2,1,2,ID,Iz,IF]
    a22 = dXdxI[2,2,1,ID,Iz,IF] + dXdxI[2,2,2,ID,Iz,IF]
    JInv = FT(1) / (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF])
    RhoL = Rho[1,Iz,ind]
    MInv = -KoeffDiv / (M[Iz,ind,1] + M[Iz,ind,2])
  end

  for it = 1 : NumTr
    if Inside
      TrCol[I,J,iz] = Cache[1,Iz,ind,it] 
    end
    @synchronize

    if Inside
      Dxc = D[I,1] * TrCol[1,J,iz]
      Dyc = D[J,1] * TrCol[I,1,iz]
      for k = 2 : N
        Dxc += D[I,k] * TrCol[k,J,iz]
        Dyc += D[J,k] * TrCol[I,k,iz] 
      end
      GradDx = (a11 * Dxc + a21 * Dyc) * JInv
      GradDy = (a12 * Dxc + a22 * Dyc) * JInv
      TrCxCol[I,J,iz] = RhoL * (a11 * GradDx + a12 * GradDy)
      TrCyCol[I,J,iz] = RhoL * (a21 * GradDx + a22 * GradDy)
    end

    @synchronize 

    if Inside
      DivTr = DW[I,1] * TrCxCol[1,J,iz] + DW[J,1] * TrCyCol[I,1,iz]
      for k = 2 : N
        DivTr += DW[I,k] * TrCxCol[k,J,iz] + DW[J,k] * TrCyCol[I,k,iz]
      end
      @atomic :monotonic FTr[1,Iz,ind,it] += DivTr * MInv
    end
  end
end

@kernel inbounds = true function HyperViscWKoeffKernel!(Fw,@Const(w),@Const(D),@Const(DW),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob),KoeffDivW) 

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  wCol = @localmem eltype(Fw) (N,N, ColumnTilesDim)
  wCxCol = @localmem eltype(Fw) (N,N, ColumnTilesDim)
  wCyCol = @localmem eltype(Fw) (N,N, ColumnTilesDim)
  if Iz < Nz && IF <= NF
    wCol[I,J,iz] = w[Iz,ind] 
  end
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz < Nz && IF <= NF
    DxcF = D[I,1] * wCol[1,J,iz]
    DycF = D[J,1] * wCol[I,1,iz]
    for k = 2 : N
      DxcF += D[I,k] * wCol[k,J,iz]
      DycF += D[J,k] * wCol[I,k,iz] 
    end
    GradDx = ((dXdxI[1,1,2,ID,Iz,IF] + dXdxI[1,1,1,ID,Iz+1,IF]) * DxcF +
      (dXdxI[2,1,2,ID,Iz,IF] + dXdxI[2,1,1,ID,Iz+1,IF]) * DycF) / (JJ[ID,2,Iz,IF] + JJ[ID,1,Iz+1,IF])
    GradDy = ((dXdxI[1,2,2,ID,Iz,IF] + dXdxI[1,2,1,ID,Iz+1,IF]) * DxcF +
      (dXdxI[2,2,2,ID,Iz,IF] + dXdxI[2,2,1,ID,Iz+1,IF]) * DycF) / (JJ[ID,2,Iz,IF] + JJ[ID,1,Iz+1,IF])
    tempx = (dXdxI[1,1,2,ID,Iz,IF] + dXdxI[1,1,1,ID,Iz+1,IF]) * GradDx +
      (dXdxI[1,2,2,ID,Iz,IF] + dXdxI[1,2,1,ID,Iz+1,IF]) * GradDy
    tempy = (dXdxI[2,1,2,ID,Iz,IF] + dXdxI[2,1,1,ID,Iz+1,IF]) * GradDx +
      (dXdxI[2,2,2,ID,Iz,IF] + dXdxI[2,2,1,ID,Iz+1,IF]) * GradDy
    wCxCol[I,J,iz] = tempx
    wCyCol[I,J,iz] = tempy
  end

  @synchronize 

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]
  if Iz < Nz && IF <= NF
    Divw = DW[I,1] * wCxCol[1,J,iz] + DW[J,1] * wCyCol[I,1,iz]
    for k = 2 : N
      Divw += DW[I,k] * wCxCol[k,J,iz] + DW[J,k] * wCyCol[I,k,iz]
    end
    @atomic :monotonic Fw[Iz,ind] += -KoeffDivW * Divw / (M[Iz,ind,2] + M[Iz+1,ind,1])
  end
end


@kernel inbounds = true function HyperViscWEDMFKernel!(Fw,@Const(w),@Const(Rho),@Const(D),@Const(DW),@Const(dXdxI),
  @Const(JJ),@Const(W),@Const(Glob)) 

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF,IE = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]
  NE = @uniform @ndrange()[5]

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  wCol = @localmem eltype(Fw) (N,N, ColumnTilesDim)
  wCxCol = @localmem eltype(Fw) (N,N, ColumnTilesDim)
  wCyCol = @localmem eltype(Fw) (N,N, ColumnTilesDim)
  if Iz < Nz && IF <= NF && IE <= NE
    wCol[I,J,iz] = w[Iz,ind,IE] 
  end
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz < Nz && IF <= NF && IE <= NE
    DxcF = D[I,1] * wCol[1,J,iz]
    DycF = D[J,1] * wCol[I,1,iz]
    for k = 2 : N
      DxcF += D[I,k] * wCol[k,J,iz]
      DycF += D[J,k] * wCol[I,k,iz] 
    end
    GradDx = ((dXdxI[1,1,2,ID,Iz,IF] + dXdxI[1,1,1,ID,Iz+1,IF]) * DxcF +
      (dXdxI[2,1,2,ID,Iz,IF] + dXdxI[2,1,1,ID,Iz+1,IF]) * DycF) / (JJ[ID,2,Iz,IF] + JJ[ID,1,Iz+1,IF])
    GradDy = ((dXdxI[1,2,2,ID,Iz,IF] + dXdxI[1,2,1,ID,Iz+1,IF]) * DxcF +
      (dXdxI[2,2,2,ID,Iz,IF] + dXdxI[2,2,1,ID,Iz+1,IF]) * DycF) / (JJ[ID,2,Iz,IF] + JJ[ID,1,Iz+1,IF])
    tempx = (dXdxI[1,1,2,ID,Iz,IF] + dXdxI[1,1,1,ID,Iz+1,IF]) * GradDx +
      (dXdxI[1,2,2,ID,Iz,IF] + dXdxI[1,2,1,ID,Iz+1,IF]) * GradDy
    tempy = (dXdxI[2,1,2,ID,Iz,IF] + dXdxI[2,1,1,ID,Iz+1,IF]) * GradDx +
      (dXdxI[2,2,2,ID,Iz,IF] + dXdxI[2,2,1,ID,Iz+1,IF]) * GradDy
    wCxCol[I,J,iz] = tempx
    wCyCol[I,J,iz] = tempy
  end

  @synchronize 

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]
  if Iz < Nz && IF <= NF && IE <= NE
    Divw = DW[I,1] * wCxCol[1,J,iz] + DW[J,1] * wCyCol[I,1,iz]
    for k = 2 : N
      Divw += DW[I,k] * wCxCol[k,J,iz] + DW[J,k] * wCyCol[I,k,iz]
    end
    @atomic :monotonic Fw[Iz,ind,IE] += Divw / (M[Iz,ind,2] + M[Iz+1,ind,1]) 
  end
end

@kernel inbounds = true function HyperViscWKoeffEDMFKernel!(Fw,@Const(w),@Const(D),@Const(DW),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob),KoeffDivW) 

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF,IE = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]
  NE = @uniform @ndrange()[5]

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  wCol = @localmem eltype(Fw) (N,N, ColumnTilesDim)
  wCxCol = @localmem eltype(Fw) (N,N, ColumnTilesDim)
  wCyCol = @localmem eltype(Fw) (N,N, ColumnTilesDim)
  if Iz < Nz && IF <= NF && I <= NE
    wCol[I,J,iz] = w[Iz,ind,IE] 
  end
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz < Nz && IF <= NF && IE <= NE
    DxcF = D[I,1] * wCol[1,J,iz]
    DycF = D[J,1] * wCol[I,1,iz]
    for k = 2 : N
      DxcF += D[I,k] * wCol[k,J,iz]
      DycF += D[J,k] * wCol[I,k,iz] 
    end
    GradDx = ((dXdxI[1,1,2,ID,Iz,IF] + dXdxI[1,1,1,ID,Iz+1,IF]) * DxcF +
      (dXdxI[2,1,2,ID,Iz,IF] + dXdxI[2,1,1,ID,Iz+1,IF]) * DycF) / (JJ[ID,2,Iz,IF] + JJ[ID,1,Iz+1,IF])
    GradDy = ((dXdxI[1,2,2,ID,Iz,IF] + dXdxI[1,2,1,ID,Iz+1,IF]) * DxcF +
      (dXdxI[2,2,2,ID,Iz,IF] + dXdxI[2,2,1,ID,Iz+1,IF]) * DycF) / (JJ[ID,2,Iz,IF] + JJ[ID,1,Iz+1,IF])
    tempx = (dXdxI[1,1,2,ID,Iz,IF] + dXdxI[1,1,1,ID,Iz+1,IF]) * GradDx +
      (dXdxI[1,2,2,ID,Iz,IF] + dXdxI[1,2,1,ID,Iz+1,IF]) * GradDy
    tempy = (dXdxI[2,1,2,ID,Iz,IF] + dXdxI[2,1,1,ID,Iz+1,IF]) * GradDx +
      (dXdxI[2,2,2,ID,Iz,IF] + dXdxI[2,2,1,ID,Iz+1,IF]) * GradDy
    wCxCol[I,J,iz] = tempx
    wCyCol[I,J,iz] = tempy
  end

  @synchronize 

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]
  if Iz < Nz && IF <= NF && IE <= NE
    Divw = DW[I,1] * wCxCol[1,J,iz] + DW[J,1] * wCyCol[I,1,iz]
    for k = 2 : N
      Divw += DW[I,k] * wCxCol[k,J,iz] + DW[J,k] * wCyCol[I,k,iz]
    end
    @atomic :monotonic Fw[Iz,ind,IE] += -KoeffDivW * Divw / (M[Iz,ind,2] + M[Iz+1,ind,1])
  end
end

@kernel inbounds = true function HyperViscTracerEDMFKernel!(FTr,@Const(Tr),@Const(D),@Const(DW),@Const(dXdxI),
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

  TrCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  TrCxCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  TrCyCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  if Iz <= Nz && IF <= NF && IE <= NE
    TrCol[I,J,iz] = Tr[Iz,ind,IE] 
  end
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz <= Nz && IF <= NF && IE <= NE
    Dxc = D[I,1] * TrCol[1,J,iz]
    Dyc = D[J,1] * TrCol[I,1,iz]
    for k = 2 : N
      Dxc += D[I,k] * TrCol[k,J,iz]
      Dyc += D[J,k] * TrCol[I,k,iz] 
    end
    @views (GradDx, GradDy) = Grad12(Dxc,Dyc,dXdxI[1:2,1:2,:,ID,Iz,IF],JJ[ID,:,Iz,IF])
    @views (tempx, tempy) = Contra12(GradDx,GradDy,dXdxI[1:2,1:2,:,ID,Iz,IF])
    TrCxCol[I,J,iz] = tempx
    TrCyCol[I,J,iz] = tempy
  end

  @synchronize 

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]
  if Iz <= Nz && IF <= NF && IE <= NE
    DivTr = DW[I,1] * TrCxCol[1,J,iz] + DW[J,1] * TrCyCol[I,1,iz]
    for k = 2 : N
      DivTr += DW[I,k] * TrCxCol[k,J,iz] + DW[J,k] * TrCyCol[I,k,iz]
    end
    @atomic :monotonic FTr[Iz,ind,IE] += DivTr / (M[Iz,ind,1] + M[Iz,ind,2])
  end
end


@kernel inbounds = true function HyperViscTracerKoeffEDMFKernel!(FTr,@Const(Cache),@Const(D),@Const(DW),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob),KoeffDiv) 

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF,IE = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]
  NE = @uniform @ndrange()[5]

  TrCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  TrCxCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  TrCyCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz <= Nz && IF <= NF && IE <= NE
    TrCol[I,J,iz] = Cache[Iz,ind,IE] 
  end
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz <= Nz && IF <= NF && IE <= NE
    Dxc = D[I,1] * TrCol[1,J,iz]
    Dyc = D[J,1] * TrCol[I,1,iz]
    for k = 2 : N
      Dxc += D[I,k] * TrCol[k,J,iz]
      Dyc += D[J,k] * TrCol[I,k,iz] 
    end
    @views (GradDx, GradDy) = Grad12(Dxc,Dyc,dXdxI[1:2,1:2,:,ID,Iz,IF],JJ[ID,:,Iz,IF])
    @views (tempx, tempy) = Contra12(GradDx,GradDy,dXdxI[1:2,1:2,:,ID,Iz,IF])
    TrCxCol[I,J,iz] = tempx
    TrCyCol[I,J,iz] = tempy
  end

  @synchronize 

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]
  if Iz <= Nz && IF <= NF && IE <= NE
    DivTr = DW[I,1] * TrCxCol[1,J,iz] + DW[J,1] * TrCyCol[I,1,iz]
    for k = 2 : N
      DivTr += DW[I,k] * TrCxCol[k,J,iz] + DW[J,k] * TrCyCol[I,k,iz]
    end
    @atomic :monotonic FTr[Iz,ind,IE] += -KoeffDiv * DivTr / (M[Iz,ind,1] + M[Iz,ind,2])
  end
end














